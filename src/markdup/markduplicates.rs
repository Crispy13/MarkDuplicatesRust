use std::{collections::HashSet, fmt::Write, mem::size_of};

use rust_htslib::bam::{
    ext::BamRecordExtensions,
    record::{Aux, RecordExt, SAMReadGroupRecord},
    HeaderView, IndexedReader, Read, Record,
};

use crate::{
    cmdline::cli::Cli,
    hts::{
        duplicate_scoring_strategy::DuplicateScoringStrategy,
        utils::sorting_collection::SortingCollection, SAMTag, SortOrder,
    },
    markdup::{
        estimate_library_complexity::EstimateLibraryComplexity,
        umi_utils::UmiUtil,
        utils::{
            library_id_generator::LibraryIdGenerator,
            read_ends::{ReadEnds, ReadEndsExt},
            read_ends_md_map::{
                DiskBasedReadEndsForMarkDuplicatesMap, MemoryBasedReadEndsForMarkDuplicatesMap,
                ReadEndsForMarkDuplicatesMap,
            },
            read_name_parser::ReadNameParserExt,
        },
    },
    utils::{hash_code, logging::ProgressLogger, CommonHasher},
};

use super::{
    markduplicates_helper::MarkDuplicatesHelper,
    utils::{
        optical_duplicate_finder::OpticalDuplicateFinder, physical_location::PhysicalLocation,
        read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates,
        read_ends_for_mark_duplicates_with_barcodes::ReadEndsBarcodeData,
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

    frag_sort: SortingCollection<ReadEndsForMarkDuplicates>,
    pair_sort: SortingCollection<ReadEndsForMarkDuplicates>,
}

impl MarkDuplicates {
    pub(crate) fn new(
        cli: Cli,
    ) -> Self {
        const size_in_bytes: usize = size_of::<ReadEndsForMarkDuplicates>();
        let max_in_memory = ((24 * 10_usize.pow(9)) as f64 * cli.SORTING_COLLECTION_SIZE_RATIO
            / size_in_bytes as f64) as usize;
        let tmp_dir = cli.TMP_DIR.clone();

        Self {
            pg_ids_seen: HashSet::new(),
            library_id_generator: LibraryIdGenerator::new(),
            optical_duplicate_finder: OpticalDuplicateFinder::default(),
            cli,
            calc_helper: MarkDuplicatesHelper::MarkDuplicatesHelper,
            frag_sort: SortingCollection::new(max_in_memory, false, tmp_dir.clone()),
            pair_sort: SortingCollection::new(max_in_memory, false, tmp_dir),
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
                MarkDuplicatesHelper::build_read_ends_for_flow(self, header, index, rec, use_barcode)
            }
        }
    }

    fn build_sorted_read_end_vecs(&mut self, use_barcodes: bool) -> Result<(), Error> {
        let mut indexed_reader = self.open_inputs()?;

        let size_in_bytes = size_of::<ReadEndsForMarkDuplicates>();

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
                    .unwrap();

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

                        paired_ends.score += self.get_read_duplicate_score(&rec, &paired_ends);
                        self.pair_sort.add(paired_ends)?;
                    }
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

        Ok(IndexedReader::from_path(self.cli.INPUT.first().unwrap())?)
    }

    fn get_read_duplicate_score(
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

        todo!()
    }
}

trait MarkDuplicatesExt {
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
