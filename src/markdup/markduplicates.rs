use std::{
    borrow::{Borrow, BorrowMut},
    collections::{HashMap, HashSet},
    fmt::Write,
    mem::size_of,
    path::Path,
};

use clap::ValueEnum;
use rust_htslib::bam::{
    ext::BamRecordExtensions, extra_ext::RecordExt, header::HeaderRecord, record::Aux, Header, HeaderView, IndexedReader, Read, Record, Writer
};

use macro_sup::set_mlog;

set_mlog!(MarkDuplicates::log);

use crate::{
    cmdline::cli::Cli,
    hts::{
        duplicate_scoring_strategy::{DuplicateScoringStrategy, ScoringStrategy},
        header::PgIdGenerator,
        metrics::{self, MetricsFile, MetricsHeader},
        utils::{
            histogram::{f64H, Histogram},
            sorting_collection::{CowForSC, SortingCollection},
            sorting_long_collection::{SortingLongCollection, SortingLongCollectionDrain},
        },
        SAMTag, SortOrder,
    },
    markdup::{
        estimate_library_complexity::EstimateLibraryComplexity,
        umi_utils::UmiUtil,
        utils::{
            library_id_generator::{DuplicationMetrics, LibraryIdGenerator},
            read_ends::{ReadEnds, ReadEndsExt},
            read_ends_md_map::{
                DiskBasedReadEndsForMarkDuplicatesMap, MemoryBasedReadEndsForMarkDuplicatesMap,
                ReadEndsForMarkDuplicatesMap,
            },
            read_name_parser::ReadNameParserExt,
            representative_read_indexer,
        },
    },
    utils::{
        hash_code,
        logging::{init_global_logger, ProgressLogger},
        CommonHasher,
    },
};

use super::{
    markduplicates_helper::MarkDuplicatesHelper,
    utils::{
        optical_duplicate_finder::OpticalDuplicateFinder, physical_location::PhysicalLocation,
        read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates,
        read_ends_for_mark_duplicates_with_barcodes::ReadEndsBarcodeData,
        representative_read_indexer::RepresentativeReadIndexer,
    },
};

// type Error = Box<dyn std::error::Error + Send + Sync>;
use anyhow::{anyhow, Error};

pub(crate) struct MarkDuplicates {
    pub(super) pg_ids_seen: HashSet<String>,
    pub(super) library_id_generator: LibraryIdGenerator,
    pub(super) optical_duplicate_finder: OpticalDuplicateFinder,
    pub(super) cli: Cli,

    calc_helper: MarkDuplicatesHelper,

    pub(super) frag_sort: SortingCollection<ReadEndsForMarkDuplicates>,
    pub(super) pair_sort: SortingCollection<ReadEndsForMarkDuplicates>,

    pub(super) duplicate_indexes: SortingLongCollection,
    pub(super) optical_duplicate_indexes: Option<SortingLongCollection>,

    pub(super) representative_read_indices_for_duplicates:
        SortingCollection<RepresentativeReadIndexer>,

    num_duplicate_indices: usize,

    default_headers:Option<Vec<MetricsHeader>>,
}

impl MarkDuplicates {
    const DUPLICATE_TYPE_TAG: &'static str = "DT";
    const DUPLICATE_TYPE_LIBRARY: &'static str = "LB";
    const DUPLICATE_TYPE_SEQUENCING: &'static str = "SQ";
    const DUPLICATE_SET_INDEX_TAG: &'static str = "DI";
    const DUPLICATE_SET_SIZE_TAG: &'static str = "DS";

    pub(crate) fn new(mut cli: Cli) -> Self {
        const size_in_bytes: usize = size_of::<ReadEndsForMarkDuplicates>();
        let max_in_memory = ((24 * 10_usize.pow(9)) as f64 * cli.SORTING_COLLECTION_SIZE_RATIO
            / size_in_bytes as f64) as usize;
        let tmp_dir = cli.TMP_DIR.clone();

        cli.DUPLICATE_SCORING_STRATEGY = ScoringStrategy::SumOfBaseQualities;

        Self {
            pg_ids_seen: HashSet::new(),
            library_id_generator: LibraryIdGenerator::new(),
            optical_duplicate_finder: OpticalDuplicateFinder::default(),
            cli,
            calc_helper: MarkDuplicatesHelper::MarkDuplicatesHelper,
            frag_sort: SortingCollection::new(max_in_memory, false, tmp_dir.clone()),
            pair_sort: SortingCollection::new(max_in_memory, false, tmp_dir),

            duplicate_indexes: Default::default(),
            optical_duplicate_indexes: Default::default(),
            representative_read_indices_for_duplicates: Default::default(),
            num_duplicate_indices: 0,
        }
    }

    #[allow(non_upper_case_globals)]
    const log: &'static str = stringify!(MarkDuplicates);
    const NO_SUCH_INDEX: i64 = i64::MAX;

    pub(super) fn build_read_ends(
        &mut self,
        header: &HeaderView,
        index: i64,
        rec: &Record,
        use_barcode: bool,
    ) -> Result<ReadEndsForMarkDuplicates, Error> {
        match self.calc_helper {
            MarkDuplicatesHelper::MarkDuplicatesHelper => {
                MarkDuplicatesHelper::build_read_ends_normal(self, header, index, rec, use_barcode)
            }
            MarkDuplicatesHelper::MarkDuplicatesForFlowHelper => {
                MarkDuplicatesHelper::build_read_ends_for_flow(
                    self,
                    header,
                    index,
                    rec,
                    use_barcode,
                )
            }
        }
    }

    fn build_sorted_read_end_vecs(&mut self, use_barcodes: bool) -> Result<(), Error> {
        let mut indexed_reader = self.open_inputs()?;

        #[allow(non_upper_case_globals)]
        const size_in_bytes: usize = size_of::<ReadEndsForMarkDuplicates>();

        // TODO: MAX_RECORDS_IN_RAM~

        let mut duplicate_query_name = String::new();
        let mut duplicate_index = Self::NO_SUCH_INDEX;

        let assumed_sort_order = SortOrder::from_header(indexed_reader.header())?;

        indexed_reader.fetch(".")?;

        let mut rec = Record::new();

        let mut tmp = match assumed_sort_order {
            SortOrder::QueryName => ReadEndsForMarkDuplicatesMap::MemoryBased(
                MemoryBasedReadEndsForMarkDuplicatesMap::new(),
            ),
            _ => {
                ReadEndsForMarkDuplicatesMap::DiskBased(DiskBasedReadEndsForMarkDuplicatesMap::new(
                    self.cli.MAX_FILE_HANDLES_FOR_READ_ENDS_MAP,
                )?)
            }
        };

        let mut index = 0;

        let mut progress = ProgressLogger::new(Self::log, 10_usize.pow(6), "Read", "records");

        while let Some(read_res) = indexed_reader.read(&mut rec) {
            // Just unwrap `Result` to check the reading process successful.
            // If Err, it means that bam file may be corrupted. Don't need to recover this error.
            read_res.unwrap();

            let qname = std::str::from_utf8(rec.qname())?;

            // This doesn't have anything to do with building sorted ReadEnd lists, but it can be done in the same pass
            // over the input
            if self.cli.PROGRAM_RECORD_ID.is_some() {
                // Gather all PG IDs seen in merged input files in first pass.  These are gathered for two reasons:
                // - to know how many different PG records to create to represent this program invocation.
                // - to know what PG IDs are already used to avoid collisions when creating new ones.
                // Note that if there are one or more records that do not have a PG tag, then a null value
                // will be stored in this set.
                match rec.aux(SAMTag::PG.name().as_bytes()) {
                    Ok(aux) => match aux {
                        rust_htslib::bam::record::Aux::String(pg) => {
                            if !self.pg_ids_seen.contains(pg) {
                                self.pg_ids_seen.insert(pg.to_string());
                            }
                        }
                        _ => Err(rust_htslib::errors::Error::BamAuxParsingError)?,
                    },
                    Err(err) => {
                        if !matches!(err, rust_htslib::errors::Error::BamAuxTagNotFound) {
                            Err(err)?
                        }
                    }
                }
            }

            // If working in query-sorted, need to keep index of first record with any given query-name.
            if matches!(assumed_sort_order, SortOrder::QueryName)
                && !rec.qname().eq(duplicate_query_name.as_bytes())
            {
                duplicate_query_name.clear();
                duplicate_query_name.push_str(qname);

                duplicate_index = index;
            }

            if rec.is_unmapped() {
                if rec.tid() == -1 && matches!(assumed_sort_order, SortOrder::Coordinate) {
                    // When we hit the unmapped reads with no coordinate, no reason to continue (only in coordinate sort).
                    break;
                }

                // If this read is unmapped but sorted with the mapped reads, just skip it.
            } else if !(rec.is_secondary() || rec.is_supplementary()) {
                let index_for_read = match assumed_sort_order {
                    SortOrder::QueryName => duplicate_index,
                    _ => index,
                };
                let fragment_end =
                    self.build_read_ends(&rec.header().unwrap(), index, &rec, use_barcodes)?;

                if rec.is_paired() && !rec.is_mate_unmapped() {
                    let mut key = String::new();

                    write!(
                        &mut key,
                        "{}{}",
                        rec.get_read_group()?.get_read_group_id(),
                        qname
                    )
                    .unwrap(); // Writing to string is infallable.

                    let paired_ends_opt = tmp.remove(rec.tid(), &key)?;

                    let mut paired_ends;
                    // See if we've already seen the first end or not
                    if paired_ends_opt.is_none() {
                        // at this point pairedEnds and fragmentEnd are the same, but we need to make
                        // a copy since pairedEnds will be modified when the mate comes along.
                        paired_ends = fragment_end.clone();

                        self.frag_sort.add(fragment_end)?;

                        tmp.put(
                            paired_ends.read_ends.read2_reference_index,
                            key,
                            paired_ends,
                        )?;
                    } else {
                        paired_ends = paired_ends_opt.unwrap();

                        let mates_ref_index = fragment_end.read_ends.read1_reference_index;
                        let mates_coordinate = fragment_end.read_ends.read1_coordinate;

                        self.frag_sort.add(fragment_end)?;

                        // Set orientationForOpticalDuplicates, which always goes by the first then the second end for the strands.  NB: must do this
                        // before updating the orientation later.
                        if rec.is_first_in_template() {
                            paired_ends.read_ends.orientation_for_optical_duplicates =
                                ReadEnds::get_orientation_bytes(
                                    rec.is_reverse(),
                                    paired_ends.read_ends.orientation == ReadEnds::R,
                                );

                            if use_barcodes {
                                paired_ends.barcode_data.as_mut().unwrap().read_one_barcode =
                                    self.get_read_one_barcode_value(&rec);
                            }
                        } else {
                            paired_ends.read_ends.orientation_for_optical_duplicates =
                                ReadEnds::get_orientation_bytes(
                                    paired_ends.read_ends.orientation == ReadEnds::R,
                                    rec.is_reverse(),
                                );

                            if use_barcodes {
                                paired_ends.barcode_data.as_mut().unwrap().read_two_barcode =
                                    self.get_read_two_barcode_value(&rec);
                            }
                        }

                        // If the other read is actually later, simply add the other read's data as read2, else flip the reads
                        if mates_ref_index > paired_ends.read_ends.read1_reference_index
                            || (mates_ref_index == paired_ends.read_ends.read1_reference_index
                                && mates_coordinate >= paired_ends.read_ends.read1_coordinate)
                        {
                            paired_ends.read_ends.read2_reference_index = mates_ref_index;
                            paired_ends.read_ends.read2_coordinate = mates_coordinate;
                            paired_ends.read2_index_in_file = index_for_read;
                            paired_ends.read_ends.orientation = ReadEnds::get_orientation_bytes(
                                paired_ends.read_ends.orientation == ReadEnds::R,
                                rec.is_reverse(),
                            );

                            // if the two read ends are in the same position, pointing in opposite directions,
                            // the orientation is undefined and the procedure above
                            // will depend on the order of the reads in the file.
                            // To avoid this, we set it explicitly (to FR):
                            if paired_ends.read_ends.read2_reference_index
                                == paired_ends.read_ends.read1_reference_index
                                && paired_ends.read_ends.read2_coordinate
                                    == paired_ends.read_ends.read1_coordinate
                                && paired_ends.read_ends.orientation == ReadEnds::RF
                            {
                                paired_ends.read_ends.orientation = ReadEnds::FR;
                            }
                        } else {
                            paired_ends.read_ends.read2_reference_index =
                                paired_ends.read_ends.read1_reference_index;
                            paired_ends.read_ends.read2_coordinate =
                                paired_ends.read_ends.read1_coordinate;
                            paired_ends.read2_index_in_file = paired_ends.read1_index_in_file;

                            paired_ends.read_ends.read1_reference_index = mates_ref_index;
                            paired_ends.read_ends.read1_coordinate = mates_coordinate;
                            paired_ends.read1_index_in_file = index_for_read;
                            paired_ends.read_ends.orientation = ReadEnds::get_orientation_bytes(
                                rec.is_reverse(),
                                paired_ends.read_ends.orientation == ReadEnds::R,
                            );
                        }

                        paired_ends.score +=
                            self.calc_helper
                                .get_read_duplicate_score(self, &rec, &paired_ends);
                        self.pair_sort.add(paired_ends)?;
                    }
                } else {
                    self.frag_sort.add(fragment_end)?;
                }
            }

            // Print out some stats every 1m reads
            index += 1;
            if progress.record(&rec) {
                log::info!(target: Self::log,
                    "Tracking {} as yet unmatched pairs. {} records in RAM.", tmp.size(), tmp.size_in_ram(),
                );
            }
        }

        log::info!(target: Self::log, "Read {} records. {} pairs never matched.", index, tmp.size());

        // Tell these collections to free up memory if possible.
        self.pair_sort.done_adding()?;
        self.frag_sort.done_adding()?;

        Ok(())
    }

    pub(super) fn get_read_one_barcode_value(&self, rec: &Record) -> u64 {
        EstimateLibraryComplexity::get_read_barcode_value(rec, &self.cli.READ_ONE_BARCODE_TAG)
    }

    pub(super) fn get_read_two_barcode_value(&self, rec: &Record) -> u64 {
        EstimateLibraryComplexity::get_read_barcode_value(rec, &self.cli.READ_TWO_BARCODE_TAG)
    }

    fn open_inputs(&self) -> Result<IndexedReader, Error> {
        assert!(
            self.cli.INPUT.len() == 1,
            "Running with multiple input bams is not supported yet."
        );

        let mut r = IndexedReader::from_path(self.cli.INPUT.first().unwrap())?;
        r.set_threads(self.cli.THREADS)?;

        Ok(r)
    }

    pub(super) fn get_read_duplicate_score(
        &self,
        rec: &Record,
        paired_ends: &ReadEndsForMarkDuplicates,
    ) -> i16 {
        DuplicateScoringStrategy::compute_duplicate_score(
            rec,
            self.cli.DUPLICATE_SCORING_STRATEGY,
            false,
        )
        .unwrap() // we can use unwrap() with assume_mate_cigar = false;
    }

    /**
     * Main work method.  Reads the SAM file once and collects sorted information about
     * the 5' ends of both ends of each read (or just one end in the case of pairs).
     * Then makes a pass through those determining duplicates before re-reading the
     * input file and writing it out with duplication flags set correctly.
     */
    pub(crate) fn do_work(&mut self) {
        let use_barcodes = !self.cli.BARCODE_TAG.is_empty()
            || !self.cli.READ_ONE_BARCODE_TAG.is_empty()
            || !self.cli.READ_TWO_BARCODE_TAG.is_empty();

        // use flow based calculation helper?
        if self.cli.FLOW_MODE {
            self.calc_helper = MarkDuplicatesHelper::MarkDuplicatesForFlowHelper;
            MarkDuplicatesHelper::validate_flow_parameters(&self);
        }

        log::info!(target: Self::log, "Reading input file and constructing read end information.");
        self.build_sorted_read_end_vecs(use_barcodes).unwrap();

        self.generate_duplicate_indexes(
            use_barcodes,
            self.cli.REMOVE_SEQUENCING_DUPLICATES
                || self.cli.TAGGING_POLICY != DuplicateTaggingPolicy::DontTag,
        );

        mlog::info!(
            "Marking {}  records as duplicates.",
            self.num_duplicate_indices
        );

        if self.cli.READ_NAME_REGEX.is_empty() {
            mlog::warn!("Skipped optical duplicate cluster discovery; library size estimation may be inaccurate!");
        } else {
            mlog::info!(
                "Found {} optical duplicate clusters.",
                self.library_id_generator
                    .get_number_of_optical_duplicate_clusters()
            )
        }

        let mut indexed_reader = self.open_inputs().unwrap();
        let header = indexed_reader.header();
        let sort_order =
            SortOrder::from_str(header.header_map().get_sort_order().unwrap()).unwrap();

        let mut output_header = Header::from_template(header);

        mlog::info!("Reads are assumed to be ordered by: {}", sort_order);

        // Setting the ASSUME_SORT_ORDER to equal queryname is understood to mean that
        // the input is
        // queryname **grouped**. So that's what we set the output order to be, so that
        // the validation will pass
        if let Some(true) = self
            .cli
            .ASSUME_SORT_ORDER
            .as_ref()
            .and_then(|v| Some(matches!(v, SortOrder::QueryName)))
        {
            let mut header_record = HeaderRecord::new(b"A");
            header_record.push_tag(SortOrder::GO.as_bytes(), "query");
            header_record.push_tag(SortOrder::SO.as_bytes(), "unknown");
            output_header.push_record(&header_record);
        }

        if self.cli.ASSUME_SORT_ORDER.is_none()
            && !matches!(sort_order, SortOrder::Coordinate)
            && !matches!(sort_order, SortOrder::QueryName)
            || self.cli.ASSUME_SORT_ORDER.is_some()
                && matches!(
                    self.cli.ASSUME_SORT_ORDER.as_ref().unwrap(),
                    SortOrder::Coordinate
                )
                && matches!(
                    self.cli.ASSUME_SORT_ORDER.as_ref().unwrap(),
                    SortOrder::QueryName
                )
        {
            panic!(
                "This program requires input that are either coordinate or query sorted (according to the header, or at least ASSUME_SORT_ORDER and the content.) \
                Found ASSUME_SORT_ORDER={:?} and header sortorder={}",
                self.cli.ASSUME_SORT_ORDER, sort_order
            );
        }

        self.cli.COMMENT.iter().for_each(|c| {
            output_header.push_comment(c.as_bytes());
        });

        // Key: previous PG ID on a SAM Record (or null). Value: New PG ID to replace
        // it.
        let chained_pg_ids =
            Self::get_chained_pg_ids(&mut self.cli, &self.pg_ids_seen, &mut output_header);

        let mut write_output_closure = || -> Result<(), Error> {
            let mut out_writer = Writer::from_path(
                &self.cli.OUTPUT,
                &output_header,
                rust_htslib::bam::Format::Bam,
            )?;

            out_writer.set_threads(self.cli.THREADS)?;

            let mut optical_duplicate_indexes_drain = self
                .optical_duplicate_indexes
                .as_mut()
                .and_then(|v| Some(v.drain()));

            let mut duplicate_indexed_drain = self.duplicate_indexes.drain();

            // Now copy over the file while marking all the necessary indexes as duplicates
            let mut record_in_file_index = 0_i64;
            let mut next_optical_duplicate_index = if let Some(odi) =
                optical_duplicate_indexes_drain
                    .as_mut()
                    .and_then(|it| it.next())
            {
                odi
            } else {
                Self::NO_SUCH_INDEX
            };

            let mut next_duplicate_index = duplicate_indexed_drain
                .next()
                .unwrap_or(Self::NO_SUCH_INDEX);

            // initialize variables for optional representative read tagging
            let mut representative_read_iterator = None;
            let rri;
            let mut representative_read_index_in_file = -1;
            let mut duplicate_set_size = -1;
            let mut next_read_in_duplicate_set_index = -1;

            if self.cli.TAG_DUPLICATE_SET_MEMBERS {
                representative_read_iterator =
                    Some(self.representative_read_indices_for_duplicates.iter()?);
                if let Some(rri_n) = representative_read_iterator.as_mut().unwrap().next() {
                    rri = rri_n;
                    next_read_in_duplicate_set_index = rri.read_index_in_file;
                    representative_read_index_in_file = rri.representative_read_index_in_file;
                    duplicate_set_size = rri.set_size;
                }
            }

            let mut progress =
                ProgressLogger::new(Self::log, 10_usize.pow(7), "Written", "records");
            let mut duplicate_query_name = Vec::<u8>::new();
            let mut representative_query_name = Vec::<u8>::new();

            indexed_reader.fetch(".")?;

            let mut rec = Record::new();
            while let Some(rl) = indexed_reader.read(&mut rec) {
                rl?;

                let metrics = Self::add_read_to_library_metrics(
                    &rec,
                    indexed_reader.header(),
                    &mut self.library_id_generator,
                    self.cli.FLOW_MODE,
                )?;

                // Now try and figure out the next duplicate index (if going by coordinate. if
                // going by query name, only do this
                // if the query name has changed.
                next_duplicate_index = Self::next_index_if_needed(
                    &sort_order,
                    record_in_file_index,
                    next_duplicate_index,
                    &duplicate_query_name,
                    &rec,
                    &mut duplicate_indexed_drain,
                );

                let is_duplicate = record_in_file_index == next_duplicate_index
                    || (matches!(sort_order, SortOrder::QueryName)
                        && record_in_file_index > next_duplicate_index
                        && rec.qname().eq(&duplicate_query_name));

                if is_duplicate {
                    rec.set_duplicate_read_flag(true);

                    metrics.add_duplicate_read_to_metrics(&rec)
                } else {
                    rec.set_duplicate_read_flag(false);
                }

                if let Some(odid) = optical_duplicate_indexes_drain.as_mut() {
                    next_optical_duplicate_index = Self::next_index_if_needed(
                        &sort_order,
                        record_in_file_index,
                        next_optical_duplicate_index,
                        &duplicate_query_name,
                        &rec,
                        odid,
                    )
                }

                let is_optical_duplicate = matches!(sort_order, SortOrder::QueryName)
                    && record_in_file_index > next_optical_duplicate_index
                    && rec.qname().eq(&duplicate_query_name)
                    || record_in_file_index == next_optical_duplicate_index;

                if self.cli.CLEAR_DT {
                    rec.remove_aux(Self::DUPLICATE_TYPE_TAG.as_bytes())?;
                }

                if !matches!(self.cli.TAGGING_POLICY, DuplicateTaggingPolicy::DontTag)
                    && rec.is_duplicate()
                {
                    if is_optical_duplicate {
                        rec.push_aux(
                            Self::DUPLICATE_TYPE_TAG.as_bytes(),
                            Aux::String(DuplicateType::SEQUENCING.code()),
                        )?;
                    } else if matches!(self.cli.TAGGING_POLICY, DuplicateTaggingPolicy::All) {
                        rec.push_aux(
                            Self::DUPLICATE_TYPE_TAG.as_bytes(),
                            Aux::String(DuplicateType::LIBRARY.code()),
                        )?;
                    }
                }

                // Tag any read pair that was in a duplicate set with the duplicate set size and
                // a representative read name
                if self.cli.TAG_DUPLICATE_SET_MEMBERS {
                    let need_next_representative_index = record_in_file_index
                        > next_read_in_duplicate_set_index as i64
                        && (matches!(sort_order, SortOrder::Coordinate)
                            || rec.qname().eq(&representative_query_name));

                    if need_next_representative_index {
                        if let Some(rri) = representative_read_iterator.as_mut().unwrap().next() {
                            next_read_in_duplicate_set_index = rri.read_index_in_file;
                            representative_read_index_in_file =
                                rri.representative_read_index_in_file;
                            duplicate_set_size = rri.set_size;
                        }
                    }

                    /*
                     * If this record's index is readInDuplicateSetIndex, then it is in a
                     * duplicateset.
                     * For queryname sorted data, we only have one representativeReadIndex entry per
                     * read name, so we need
                     * to also look for additional reads with the same name.
                     */
                    let is_duplicate_set = record_in_file_index
                        == next_read_in_duplicate_set_index as i64
                        || (matches!(sort_order, SortOrder::QueryName)
                            && record_in_file_index > next_read_in_duplicate_set_index as i64
                            && rec.qname().eq(&representative_query_name));
                    if is_duplicate_set {
                        if !rec.is_secondary_or_supplementary() && !rec.is_unmapped() {
                            rec.push_aux(
                                Self::DUPLICATE_SET_INDEX_TAG.as_bytes(),
                                Aux::I32(representative_read_index_in_file),
                            )?;
                            rec.push_aux(
                                Self::DUPLICATE_SET_SIZE_TAG.as_bytes(),
                                Aux::I32(duplicate_set_size),
                            )?;
                            representative_query_name = rec.qname().to_vec();
                        }
                    }
                }

                // Set MOLECULAR_IDENTIFIER_TAG for SAMRecord rec
                if !self.cli.BARCODE_TAG.is_empty() {
                    UmiUtil::set_molecular_identifier(
                        &mut rec,
                        "",
                        &self.cli.MOLECULAR_IDENTIFIER_TAG,
                        self.cli.DUPLEX_UMI,
                    );
                }

                // Note, duplicateQueryName must be incremented after we have marked both
                // optical and sequencing duplicates for queryname sorted files.
                if is_duplicate {
                    duplicate_query_name = rec.qname().to_vec();
                }

                // Output the record if desired and bump the record index
                record_in_file_index += 1;
                if self.cli.REMOVE_DUPLICATES && rec.is_duplicate() {
                    continue;
                }

                if self.cli.REMOVE_SEQUENCING_DUPLICATES && is_optical_duplicate {
                    continue;
                }

                if self.cli.PROGRAM_RECORD_ID.is_some() && self.cli.ADD_PG_TAG_TO_READS {
                    if let Ok(pg) = rec.get_str_aux(SAMTag::PG.name().as_bytes()) {
                        rec.push_aux(
                            SAMTag::PG.name().as_bytes(),
                            Aux::String(chained_pg_ids.as_ref().unwrap().get(pg).unwrap()),
                        )?;
                    }
                }

                out_writer.write(&rec)?;
                progress.record(&rec);
            }

            // remember to close the inputs
            mlog::info!("Writing complete. Closing input iterator.");

            mlog::info!("Duplicate Index cleanup.");

            drop(duplicate_indexed_drain);
            self.duplicate_indexes.clean_up();

            if self.cli.TAG_DUPLICATE_SET_MEMBERS {
                mlog::info!("Representative read Index cleanup.");
                self.representative_read_indices_for_duplicates.clean_up();
            }

            mlog::info!("Getting Memory Stats.");

            Ok(())
        };

        write_output_closure().unwrap();

        mlog::info!("Closed outputs. Getting more Memory Stats.");

        // Write out the metrics
        Self::finalize_and_write_metrics(
            &mut self.library_id_generator,
            self.get_metrics_file(),
            &self.cli.METRICS_FILE,
        );
    }

    /**
     * Goes through the accumulated ReadEndsForMarkDuplicates objects and determines
     * which of them are
     * to be marked as duplicates.
     */
    pub(crate) fn sort_indices_for_duplicates(&mut self, index_optical_duplicates: bool) {
        let entry_overhead = if self.cli.TAG_DUPLICATE_SET_MEMBERS {
            // Memory requirements for RepresentativeReadIndexer:
            // three int entries + overhead: (3 * 4) + 4 = 16 bytes
            16
        } else {
            SortingLongCollection::SIZEOF
        };

        // Keep this number from getting too large even if there is a huge heap.
        let mut max_in_memory =
            (self.cli.MAX_MEMORY as f64 * 0.25 / entry_overhead as f64).min((i32::MAX - 5) as f64);

        // If we're also tracking optical duplicates, reduce maxInMemory, since we'll
        // need two sorting collections
        if index_optical_duplicates {
            max_in_memory /=
                ((entry_overhead + SortingLongCollection::SIZEOF) / entry_overhead) as f64;
            self.optical_duplicate_indexes = Some(SortingLongCollection::new(
                max_in_memory as usize,
                Vec::with_capacity(self.cli.TMP_DIR.len()),
            ));
        }

        mlog::info!(
            "Will retain up to {} duplicate indices before spilling to disk.",
            max_in_memory
        );
        self.duplicate_indexes = SortingLongCollection::new(
            max_in_memory as usize,
            Vec::with_capacity(self.cli.TMP_DIR.len()),
        );

        if self.cli.TAG_DUPLICATE_SET_MEMBERS {
            self.representative_read_indices_for_duplicates =
                SortingCollection::new(max_in_memory as usize, false, self.cli.TMP_DIR.clone());
        }
    }

    pub(crate) fn handle_chunk(&mut self, next_chunk: &mut Vec<ReadEndsForMarkDuplicates>) {
        // if next_chunk.len() >= 1 {
        //     mlog::debug!(
        //         "handle_chunk: next_chunk.len()={} first_index={}",
        //         next_chunk.len(),
        //         next_chunk.first().unwrap().read1_index_in_file
        //     );
        // }

        // mlog::debug!(
        //     "handle_chunk: iif_indices={:#?}",
        //     next_chunk
        //         .iter()
        //         .map(|e| (e.read1_index_in_file, e.read2_index_in_file))
        //         .collect::<Vec<_>>()
        // );

        if next_chunk.len() > 1 {
            self.mark_duplicate_pairs(next_chunk.as_mut_slice());
            if self.cli.TAG_DUPLICATE_SET_MEMBERS {
                self.add_representative_read_index(&next_chunk);
            }
        } else if next_chunk.len() == 1 {
            Self::add_singleton_to_count(&mut self.library_id_generator);
        }
    }

    fn add_singleton_to_count(library_id_generator: &mut LibraryIdGenerator) {
        library_id_generator
            .get_duplicate_count_hist()
            .increment1(1.0);
        library_id_generator
            .get_non_optical_duplicate_count_hist()
            .increment1(1.0);
    }

    /**
     * Takes a list of ReadEndsForMarkDuplicates objects and identify the
     * representative read based on
     * quality score. For all members of the duplicate set, add the read1
     * index-in-file of the representative
     * read to the records of the first and second in a pair. This value becomes is
     * used for
     * the 'DI' tag.
     */
    fn add_representative_read_index(&mut self, read_ends_slice: &[ReadEndsForMarkDuplicates]) {
        let mut max_score = 0_i16;
        let mut best = None;

        // All read ends should have orientation FF, FR, RF, or RR
        for end in read_ends_slice.iter() {
            if end.score > max_score || best.is_none() {
                max_score = end.score;
                best = Some(end);
            }
        }

        let best = best.unwrap();

        // for read name (for representative read name), add the last of the pair that
        // was examined
        for end in read_ends_slice {
            self.add_representative_read_of_duplicate_set(
                best.read1_index_in_file,
                read_ends_slice.len(),
                end.read1_index_in_file,
            );

            if end.read1_index_in_file != end.read2_index_in_file {
                self.add_representative_read_of_duplicate_set(
                    best.read1_index_in_file,
                    read_ends_slice.len(),
                    end.read2_index_in_file,
                );
            }
        }
    }

    fn add_representative_read_of_duplicate_set(
        &mut self,
        representative_read_index_in_file: i64,
        set_size: usize,
        read1_index_in_file: i64,
    ) {
        let mut rri = RepresentativeReadIndexer::default();
        rri.representative_read_index_in_file = representative_read_index_in_file as i32;
        rri.set_size = set_size as i32;
        rri.read_index_in_file = read1_index_in_file as i32;

        self.representative_read_indices_for_duplicates
            .add(rri)
            .unwrap();
    }

    /**
     * Takes a list of ReadEndsForMarkDuplicates objects and removes from it all
     * objects that should
     * not be marked as duplicates. This assumes that the list contains objects
     * representing pairs.
     */
    fn mark_duplicate_pairs(&mut self, read_ends_slice: &mut [ReadEndsForMarkDuplicates]) {
        let mut max_score = 0_i16;
        let mut best = None;

        // All read ends should have orientation FF, FR, RF, or RR
        for end in read_ends_slice.iter() {
            if end.score > max_score || best.is_none() {
                max_score = end.score;
                best = Some(end);
            }
        }

        let best = best.cloned();

        if !self.cli.READ_NAME_REGEX.is_empty() {
            Self::track_optical_duplicates(
                read_ends_slice,
                best.as_ref(),
                &self.optical_duplicate_finder,
                &mut self.library_id_generator,
            )
        }

        if best.is_none() {
            return;
        }

        let best = best.unwrap();
        for end in read_ends_slice.iter() {
            // if end.read1_index_in_file == 4 || end.read2_index_in_file == 4 {
            //     println!("end={:#?}", end);
            //     println!("best={:#?}", end);
            //     panic!("early end.");
            // }

            if end != &best {
                // mlog::debug!("Marked as duplicate: {} at p1.", end.read1_index_in_file);
                self.add_index_as_duplicate(end.read1_index_in_file);

                // in query-sorted case, these will be the same.
                // TODO: also in coordinate sorted, when one read is unmapped
                if end.read2_index_in_file != end.read1_index_in_file {
                    // mlog::debug!("Marked as duplicate: {} at p2.", end.read2_index_in_file);
                    self.add_index_as_duplicate(end.read2_index_in_file);
                }

                if let Some((true, odi)) = self
                    .optical_duplicate_indexes
                    .as_mut()
                    .and_then(|v| Some((end.read_ends.is_optical_duplicate, v)))
                {
                    odi.add(end.read1_index_in_file);
                    // We expect end.read2IndexInFile==read1IndexInFile when we are in queryname
                    // sorted files, as the read-pairs
                    // will be sorted together and nextIndexIfNeeded() will only pull one index from
                    // opticalDuplicateIndexes.
                    // This means that in queryname sorted order we will only pull from the sorting
                    // collection once,
                    // where as we would pull twice for coordinate sorted files.
                    if end.read2_index_in_file != end.read1_index_in_file {
                        odi.add(end.read2_index_in_file);
                    }
                }
            }
        }
    }

    fn add_index_as_duplicate(&mut self, bam_index: i64) {
        self.duplicate_indexes.add(bam_index);
        self.num_duplicate_indices += 1;
    }
    /**
     * Looks through the set of reads and identifies how many of the duplicates are
     * in fact optical duplicates, and stores the data in the instance level histogram.
     * Additionally sets the transient isOpticalDuplicate flag on each read end that is
     * identified as an optical duplicate.
     */
    fn track_optical_duplicates(
        ends: &mut [ReadEndsForMarkDuplicates],
        keeper: Option<&ReadEndsForMarkDuplicates>,
        optical_duplicate_finder: &OpticalDuplicateFinder,
        library_id_generator: &mut LibraryIdGenerator,
    ) {
        #[allow(non_snake_case)]
        let mut has_FR = false;
        #[allow(non_snake_case)]
        let mut has_RF = false;

        // Check to see if we have a mixture of FR/RF
        for end in ends.iter() {
            if ReadEnds::FR == end.read_ends.orientation_for_optical_duplicates {
                has_FR = true;

                if has_RF {
                    break;
                }
            } else if ReadEnds::RF == end.read_ends.orientation_for_optical_duplicates {
                has_RF = true;

                if has_FR {
                    break;
                }
            }
        }

        // Check if we need to partition since the orientations could have changed
        let n_optical_dup;
        if has_FR && has_RF {
            // need to track them independently
            // Variables used for optical duplicate detection and tracking
            #[allow(non_snake_case)]
            let (mut track_optical_duplicates_F, mut track_optical_duplicates_R) =
                (Vec::new(), Vec::new());

            // Split into two lists: first of pairs and second of pairs, since they must have orientation and same starting end
            for end in ends.iter_mut() {
                if ReadEnds::FR == end.read_ends.orientation_for_optical_duplicates {
                    track_optical_duplicates_F.push(end);
                } else if ReadEnds::RF == end.read_ends.orientation_for_optical_duplicates {
                    track_optical_duplicates_R.push(end);
                } else {
                    panic!(
                        "Found an unexpected orientation: {}",
                        end.read_ends.orientation
                    )
                }
            }

            // track the duplicates
            #[allow(non_snake_case)]
            let n_optical_dup_F = Self::track_optical_duplicates_with_histo(
                track_optical_duplicates_F.as_mut_slice(),
                keeper,
                optical_duplicate_finder,
                library_id_generator.get_mut_optical_duplicates_by_library_id_map(),
            );

            #[allow(non_snake_case)]
            let n_optical_dup_R = Self::track_optical_duplicates_with_histo(
                track_optical_duplicates_R.as_mut_slice(),
                keeper,
                optical_duplicate_finder,
                library_id_generator.get_mut_optical_duplicates_by_library_id_map(),
            );

            n_optical_dup = n_optical_dup_F + n_optical_dup_R;
        } else {
            // No need to partition
            n_optical_dup = Self::track_optical_duplicates_with_histo(
                ends,
                keeper,
                optical_duplicate_finder,
                library_id_generator.get_mut_optical_duplicates_by_library_id_map(),
            );
        }

        Self::track_duplicate_counts(ends.len(), n_optical_dup, library_id_generator);
    }

    /**
     * Looks through the set of reads and identifies how many of the duplicates are
     * in fact optical duplicates, and stores the data in the instance level histogram.
     *
     * We expect only reads with FR or RF orientations, not a mixture of both.
     *
     * In PCR duplicate detection, a duplicates can be a have FR and RF when fixing the orientation order to the first end of the mate.  In
     * optical duplicate detection, we do not consider them duplicates if one read as FR and the other RF when we order orientation by the
     * first mate sequenced (read #1 of the pair).
     */
    fn track_optical_duplicates_with_histo(
        ends: &mut [impl BorrowMut<ReadEndsForMarkDuplicates> + Borrow<ReadEndsForMarkDuplicates>],
        keeper: Option<&ReadEndsForMarkDuplicates>,
        optical_duplicate_finder: &OpticalDuplicateFinder,
        optical_duplicates_by_library_id: &mut Histogram<i16>,
    ) -> i32 {
        let ends_immut = ends.iter().map(|e| e.borrow()).collect::<Vec<_>>();

        let optical_duplicate_flags =
            optical_duplicate_finder.find_optical_duplicates(ends_immut.as_slice(), keeper);

        let mut optical_duplicates = 0;

        for (odf, end) in optical_duplicate_flags.into_iter().zip(ends.iter_mut()) {
            if odf {
                optical_duplicates += 1;
                end.borrow_mut().read_ends.is_optical_duplicate = true;
            }
        }

        if optical_duplicates > 0 {
            optical_duplicates_by_library_id.increment(
                ends.first().unwrap().borrow().get_library_id(),
                optical_duplicates as f64,
            );
        }

        optical_duplicates
    }

    fn track_duplicate_counts(
        vec_len: usize,
        opt_dup_cnt: i32,
        library_id_generator: &mut LibraryIdGenerator,
    ) {
        let duplicates_count_hist = &mut library_id_generator.duplicate_count_hist;
        let non_optical_duplicates_count_hist =
            &mut library_id_generator.non_optical_duplicate_count_hist;
        let optical_duplicates_count_hist = &mut library_id_generator.optical_duplicate_count_hist;

        duplicates_count_hist.increment1(vec_len as f64);

        if vec_len as i32 - opt_dup_cnt > 0 {
            non_optical_duplicates_count_hist.increment1((vec_len as i32 - opt_dup_cnt) as f64);
        }

        if opt_dup_cnt > 0 {
            optical_duplicates_count_hist.increment1((opt_dup_cnt + 1) as f64);
        }
    }

    pub(super) fn mark_duplicate_fragments(
        &mut self,
        list: &[ReadEndsForMarkDuplicates],
        contains_pairs: bool,
    ) {
        // mlog::debug!(
        //     "mark_duplicate_fragments: list.len()={} first_index={}",
        //     list.len(),
        //     list.first().unwrap().read1_index_in_file
        // );

        if contains_pairs {
            for end in list {
                if !end.is_paired() {
                    self.add_index_as_duplicate(end.read1_index_in_file);
                }
            }
        } else {
            let mut max_score = 0_i16;
            let mut best = None;
            for end in list {
                if end.score > max_score || best.is_none() {
                    max_score = end.score;
                    best = Some(end);
                }
            }

            let best = best.unwrap();
            for end in list {
                if end != best {
                    self.add_index_as_duplicate(end.read1_index_in_file);
                }
            }
        }
    }

    /**
     * We have to re-chain the program groups based on this algorithm.  This returns the map from existing program group ID
     * to new program group ID.
     */
    fn get_chained_pg_ids<'p>(
        cli: &mut Cli,
        pg_ids_seen: &'p HashSet<String>,
        output_header: &mut Header,
    ) -> Option<HashMap<&'p str, String>> {
        let mut chained_pg_ids = HashMap::new();
        cli.PROGRAM_GROUP_NAME = "MarkDuplicates".to_string();

        // Generate new PG record(s)
        if let Some(pr_id) = cli.PROGRAM_RECORD_ID.as_ref() {
            let mut pg_id_generator = PgIdGenerator::new(output_header);
            if cli.PROGRAM_GROUP_VERSION.is_empty() {
                cli.PROGRAM_GROUP_VERSION = env!("CARGO_PKG_VERSION").to_string();
            }

            if cli.PROGRAM_GROUP_COMMAND_LINE.is_empty() {
                cli.PROGRAM_GROUP_COMMAND_LINE = "MarkDuplicatesRust".to_string();
            }

            for existing_id in pg_ids_seen.iter().map(String::as_str) {
                let new_pg_id = pg_id_generator.get_non_colliding_id(pr_id.to_string());
                chained_pg_ids.insert(existing_id, new_pg_id);

                let mut program_record = HeaderRecord::new(b"Z");
                program_record.push_tag(b"VN", &cli.PROGRAM_GROUP_VERSION);
                program_record.push_tag(b"CL", &cli.PROGRAM_GROUP_COMMAND_LINE);
                program_record.push_tag(b"PN", &cli.PROGRAM_GROUP_NAME);
                program_record.push_tag(b"ID", existing_id);
                output_header.push_record(&program_record);
            }
        } else {
            return None;
        }

        Some(chained_pg_ids)
    }

    fn generate_duplicate_indexes(&mut self, use_barcodes: bool, index_optical_duplicates: bool) {
        match self.calc_helper {
            MarkDuplicatesHelper::MarkDuplicatesHelper => {
                MarkDuplicatesHelper::generate_duplicate_indexes_normal(
                    self,
                    use_barcodes,
                    index_optical_duplicates,
                )
            }
            MarkDuplicatesHelper::MarkDuplicatesForFlowHelper => {
                MarkDuplicatesHelper::generate_duplicate_indexes_for_flow(
                    self,
                    use_barcodes,
                    index_optical_duplicates,
                )
            }
        }
    }

    fn add_read_to_library_metrics<'l>(
        rec: &Record,
        header: &HeaderView,
        library_id_generator: &'l mut LibraryIdGenerator,
        flow_metrics: bool,
    ) -> Result<&'l mut DuplicationMetrics, Error> {
        let library = LibraryIdGenerator::get_library_name(rec)?;

        if library_id_generator
            .get_mut_metrics_by_library(library)
            .is_none()
        {
            let mut metrics = DuplicationMetrics::create_metrics(flow_metrics);
            metrics.library = library.to_string();

            library_id_generator.add_metrics_by_library(library.to_string(), metrics);
        }

        let metrics = library_id_generator
            .get_mut_metrics_by_library(library)
            .unwrap();

        metrics.add_read_to_library_metrics(rec);

        Ok(metrics)
    }

    fn next_index_if_needed(
        sort_order: &SortOrder,
        record_in_file_index: i64,
        mut next_duplicate_index: i64,
        last_query_name: &[u8],
        rec: &Record,
        duplicate_indexes_drain: &mut SortingLongCollectionDrain,
    ) -> i64 {
        // Manage the flagging of optical/sequencing duplicates
        // Possibly figure out the next opticalDuplicate index (if going by coordinate,
        // if going by query name, only do this
        // if the query name has changed)
        let need_next_duplicate_index = record_in_file_index > next_duplicate_index
            && (matches!(sort_order, SortOrder::Coordinate) || !rec.qname().eq(last_query_name));

        if need_next_duplicate_index {
            next_duplicate_index = if let Some(di) = duplicate_indexes_drain.next() {
                di
            } else {
                Self::NO_SUCH_INDEX
            }
        }

        next_duplicate_index
    }

    /**
     * Writes the metrics given by the libraryIdGenerator to the outputFile.
     *
     * @param libraryIdGenerator A {@link LibraryIdGenerator} object that contains the map from library to {@link DuplicationMetrics} for
     *                           that library
     * @param metricsFile        An empty {@link MetricsFile} object that will be filled, with "finalized" metrics and written out.
     *                           It needs to be generated from a non-static context so that various commandline information is
     *                           added to the header when {@link CommandLineProgram#getMetricsFile()} is called.
     * @param outputFile         The file to write the metrics to
     */
    fn finalize_and_write_metrics(
        library_id_generator: &mut LibraryIdGenerator,
        mut metrics_file: MetricsFile<DuplicationMetrics, f64H>,
        output_file: impl AsRef<Path>,
    ) {
        let metrics_by_library = library_id_generator.get_mut_metrics_by_library_map();
        let optical_duplicates_by_library_id = library_id_generator.get_optical_duplicates_by_library_id_map();
        let library_ids = library_id_generator.get_library_ids_map();

        // Write out the metrics
        for (library_name, metrics) in metrics_by_library.into_iter() {
            metrics.READ_PAIRS_EXAMINED = metrics.READ_PAIRS_EXAMINED / 2;
            metrics.READ_PAIR_DUPLICATES = metrics.READ_PAIR_DUPLICATES / 2;

            // Add the optical dupes to the metrics
            if let Some(library_id) = library_ids.get(library_name) {
                if let Some(bin) = optical_duplicates_by_library_id.get(library_id) {
                    metrics.READ_PAIR_OPTICAL_DUPLICATES = bin.get_value() as usize;
                }
            }

            metrics.calculate_derived_fields();
            metrics_file.add_metric(metrics.clone());
        }

        if metrics_by_library.len() == 1 {
            metrics_file.set_histogram(metrics_by_library.pop_first().unwrap().1.calculate_roi_histogram());
        }

        // Add set size histograms - the set size counts are printed on adjacent columns to the ROI metric.
        metrics_file.add_histogram(library_id_generator.get_duplicate_count_hist().clone());
        metrics_file.add_histogram(library_id_generator.get_optical_duplicate_count_hist().clone());
        metrics_file.add_histogram(library_id_generator.get_non_optical_duplicate_count_hist().clone());

        metrics_file.write(output_file.as_ref());
    }

    /** Gets a MetricsFile with default headers already written into it. */
    fn get_metrics_file<B, H>(&self) -> MetricsFile<B, H> {
        let file = MetricsFile::new();

        for h in self.get_default_metrics_headers() {
            file.add_header(h.clone());
        }

        file
    }

    fn get_default_metrics_headers(&mut self) -> &Vec<MetricsHeader>  {
        if self.default_headers.is_none() {
            // set default header. code based on CommandLineProgram.java.

            // Build the default headers
            self.default_headers = Some(vec![
                MetricsHeader::StringHeader(self.cli.commandline.clone()),
                MetricsHeader::StringHeader(format!("Started on: {}", chrono::offset::Local::now().format("%Y%m%d %H:%M"))),
            ]);
        }

        self.default_headers.as_ref().unwrap()
    }
}

pub(super) trait MarkDuplicatesExt {
    fn are_comparable_for_duplicates(
        lhs: &ReadEndsForMarkDuplicates,
        rhs: &ReadEndsForMarkDuplicates,
        compare_read2: bool,
        use_barcodes: bool,
    ) -> bool {
        let mut are_comparable = lhs.read_ends.library_id == rhs.read_ends.library_id;

        // areComparable is useful here to avoid the casts below
        if use_barcodes && are_comparable {
            let lhs_with_barcodes = lhs.barcode_data.as_ref().unwrap();
            let rhs_with_barcodes = rhs.barcode_data.as_ref().unwrap();

            are_comparable = lhs_with_barcodes.barcode == rhs_with_barcodes.barcode
                && lhs_with_barcodes.read_one_barcode == rhs_with_barcodes.read_one_barcode
                && lhs_with_barcodes.read_two_barcode == rhs_with_barcodes.read_two_barcode;
        }

        if are_comparable {
            are_comparable = lhs.read_ends.read1_reference_index
                == rhs.read_ends.read1_reference_index
                && lhs.read_ends.read1_coordinate == rhs.read_ends.read1_coordinate
                && lhs.read_ends.orientation == rhs.read_ends.orientation;
        }

        if are_comparable && compare_read2 {
            are_comparable = lhs.read_ends.read2_reference_index
                == rhs.read_ends.read2_reference_index
                && lhs.read_ends.read2_coordinate == rhs.read_ends.read2_coordinate
        }

        are_comparable
    }
}

impl MarkDuplicatesExt for MarkDuplicates {}

#[derive(ValueEnum, Clone, PartialEq, Debug)]
pub(crate) enum DuplicateTaggingPolicy {
    DontTag,
    OpticalOnly,
    All,
}

pub(crate) enum DuplicateType {
    LIBRARY,
    SEQUENCING,
}

impl DuplicateType {
    fn code(&self) -> &str {
        match self {
            DuplicateType::LIBRARY => MarkDuplicates::DUPLICATE_TYPE_LIBRARY,
            DuplicateType::SEQUENCING => MarkDuplicates::DUPLICATE_TYPE_SEQUENCING,
        }
    }
}

#[cfg(test)]
mod test {
    use std::{borrow::Cow, path::Path};

    use clap::Parser;

    use crate::tests::{
        from_java_read_ends, JavaReadEndIterator, JsonFileIterator, ToJsonFileForTest,
    };

    use super::*;

    #[test]
    fn test_do_work() {
        init_global_logger(log::LevelFilter::Debug);

        let bam_path =
            Path::new("tests/data/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam");

        let out_path =
            Path::new("NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.markdup.bam");

        let args = std::env::var("MDR_ARGS").unwrap();
        let args_splited = args.trim_end().split(" ").collect::<Vec<_>>();

        println!("{:?}", args_splited);
        let mut cli = Cli::parse_from(args_splited.into_iter());
        cli.after_action();

        // println!("{:?}", std::env::args());

        // println!("Start to test.");

        // let cli = Cli::parse();

        let mut md = MarkDuplicates::new(cli);

        md.do_work();

        // let duplicate_indexes_drain = md.duplicate_indexes.drain();

        // duplicate_indexes_drain
        //     .save_object_to_json(format!("{}.DI.rust.json", md.cli.INPUT.get(0).unwrap()));

        // JsonFileIterator::<i64>::new(format!("{}.DI.java.json", md.cli.INPUT.get(0).unwrap()))
        //     .zip(duplicate_indexes_drain)
        //     .enumerate()
        //     .for_each(|(i, (j, r))| {
        //         assert_eq!(j, r, ">> i={}\njava={:#?}\n\nrust={:#?}", i, j, r);
        //     });

        // if let Some(odi) = md.optical_duplicate_indexes.as_mut() {
        //     JsonFileIterator::<i64>::new(format!("{}.ODI.java.json", md.cli.INPUT.get(0).unwrap()))
        //         .zip(odi.drain())
        //         .enumerate()
        //         .for_each(|(i, (j, r))| {
        //             assert_eq!(j, r, ">> i={}\njava={:#?}\n\nrust={:#?}", i, j, r);
        //         });
        // }

        // JsonFileIterator::<RepresentativeReadIndexer>::new(format!("{}.RRI.java.json", md.cli.INPUT.get(0).unwrap()))
        //         .zip(md.representative_read_indices_for_duplicates.drain().unwrap())
        //         .enumerate()
        //         .for_each(|(i, (j, r))| {
        //             assert_eq!(j, r, ">> i={}\njava={:#?}\n\nrust={:#?}", i, j, r);
        //         });
    }

    #[test]
    fn test_build_read_ends() {
        init_global_logger(log::LevelFilter::Debug);
        // mlog::warn!("works?");

        // NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
        // NA12878.chrom11.100.bam
        let bam_path =
            Path::new("tests/data/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam");

        let cli = Cli::parse_from([
            "",
            "--INPUT",
            bam_path.to_str().unwrap(), //
        ]);

        let mut md = MarkDuplicates::new(cli);

        md.build_sorted_read_end_vecs(false).unwrap();

        let mut ps = md.pair_sort.iter().unwrap();
        let fs = md.frag_sort.iter().unwrap();

        // let fs_queue_order_path = format!(
        //     "{}.frag_sort.queue.rust.json",
        //     bam_path.file_stem().unwrap().to_str().unwrap()
        // );
        // match fs {
        //     crate::hts::utils::sorting_collection::SortingCollectionIter::InMemoryIter(_) => {}
        //     crate::hts::utils::sorting_collection::SortingCollectionIter::MergingIterator(ref v) => {
        //         v.queue
        //             .iter()
        //             .map(|e| e.peeked().unwrap())
        //             .save_object_to_json(&fs_queue_order_path);
        //     }
        // }
        // println!("ps_item={:#?}", ps.next().unwrap());
        // println!("pair_sort(N={})={:#?}", ps.len(), ps);
        // println!("frag_sort(N={})={:#?}", fs.len(), fs);

        // let ps_json_path = format!("{}.pair_sort.rust.json", bam_path.file_stem().unwrap().to_str().unwrap());
        // let fs_json_path = format!("{}.frag_sort.rust.json", bam_path.file_stem().unwrap().to_str().unwrap());
        // ps.save_object_to_json(&ps_json_path);
        // fs.map(|e| e.into_owned()).take(2000).save_object_to_json(&fs_json_path);

        // println!("ps_json_path={}\nfs_json_pat={}", ps_json_path, fs_json_path)

        // JavaReadEndIterator::new(format!("{}.pairSort.java.json", bam_path.to_str().unwrap()))
        //     .map(from_java_read_ends)
        //     .zip(ps)
        //     .enumerate()
        //     .for_each(|(i, (j, r))| {
        //         assert_eq!(j, r, ">> i={}\njava={:#?}\n\nrust={:#?}", i, j, r);
        //     });

        // let fs = md.frag_sort.iter().unwrap();
        // JavaReadEndIterator::new(format!("{}.fragSort.java.json", bam_path.to_str().unwrap()))
        //     .map(from_java_read_ends)
        //     .zip(fs)
        //     .enumerate()
        //     .for_each(|(i, (j, r))| {
        //         assert_eq!(j, *r, ">> i={}\njava={:#?}\n\nrust={:#?}", i, j, r);
        //     });
        // assert_eq!(ps, from_java_read_ends_to_rust_read_ends("tests/data/NA12878.chrom11.20.bam.20.pairSort.read.ends"));
        // assert_eq!(fs, from_java_read_ends_to_rust_read_ends("tests/data/NA12878.chrom11.20.bam.20.fragSort.read.ends"));
    }

    #[test]
    fn test_are_comparable() {
        let a = r#"
        {
            "read1_qname": "SRR622461.127",
            "score": 4482,
            "read1_index_in_file": 2,
            "read2_index_in_file": 728,
            "duplicate_set_size": -1,
            "read_ends": {
                "library_id": 1,
                "orientation": 3,
                "read1_reference_index": 0,
                "read1_coordinate": 10001,
                "read2_reference_index": 0,
                "read2_coordinate": 10388,
                "read_group": -1,
                "orientation_for_optical_duplicates": 3,
                "pls": {
                    "tile": -1,
                    "x": -1,
                    "y": -1
                }
            },
            "barcode_data": null
        }"#;
        // let a = r#"{
        //     "read1_qname": "SRR622461.133",
        //     "score": 4593,
        //     "read1_index_in_file": 3,
        //     "read2_index_in_file": 734,
        //     "duplicate_set_size": -1,
        //     "read_ends": {
        //         "library_id": 1,
        //         "orientation": 3,
        //         "read1_reference_index": 0,
        //         "read1_coordinate": 10001,
        //         "read2_reference_index": 0,
        //         "read2_coordinate": 10388,
        //         "read_group": -1,
        //         "orientation_for_optical_duplicates": 5,
        //         "pls": {
        //             "tile": -1,
        //             "x": -1,
        //             "y": -1
        //         }
        //     },
        //     "barcode_data": null
        // }"#;

        let a = serde_json::from_str::<ReadEndsForMarkDuplicates>(a).unwrap();

        let b = r#"
        {
            "read1_qname": "SRR622461.145",
            "score": 3327,
            "read1_index_in_file": 4,
            "read2_index_in_file": 744,
            "duplicate_set_size": -1,
            "read_ends": {
                "library_id": 1,
                "orientation": 3,
                "read1_reference_index": 0,
                "read1_coordinate": 10001,
                "read2_reference_index": 0,
                "read2_coordinate": 10388,
                "read_group": -1,
                "orientation_for_optical_duplicates": 3,
                "pls": {
                    "tile": -1,
                    "x": -1,
                    "y": -1
                }
            },
            "barcode_data": null
        }"#;

        let b = serde_json::from_str::<ReadEndsForMarkDuplicates>(b).unwrap();

        println!(
            "{}",
            MarkDuplicates::are_comparable_for_duplicates(&a, &b, true, false)
        );
    }
}
