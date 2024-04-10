use std::{
    borrow::{Borrow, BorrowMut},
    collections::HashSet,
    fmt::Write,
    mem::size_of,
};

use rust_htslib::bam::{
    ext::BamRecordExtensions,
    record::{Aux, RecordExt, SAMReadGroupRecord},
    HeaderView, IndexedReader, Read, Record,
};

use macro_sup::set_mlog;

set_mlog!(MarkDuplicates::log);

use crate::{
    cmdline::cli::Cli,
    hts::{
        duplicate_scoring_strategy::{DuplicateScoringStrategy, ScoringStrategy},
        utils::{
            histogram::Histogram,
            sorting_collection::{CowForSC, SortingCollection},
            sorting_long_collection::SortingLongCollection,
        },
        SAMTag, SortOrder,
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

    pub(super) representative_read_indices_for_duplicates: SortingCollection<RepresentativeReadIndexer>,

    num_duplicate_indices: usize,
}

impl MarkDuplicates {
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

        Ok(IndexedReader::from_path(self.cli.INPUT.first().unwrap())?)
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
        init_global_logger(self.cli.LOGGING_LEVEL);

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
        if next_chunk.len() > 1 {
            self.mark_duplicate_pairs(next_chunk.as_mut_slice());
            if self.cli.TAG_DUPLICATE_SET_MEMBERS {
                self.add_representative_read_index(&next_chunk);
            }
        } else if next_chunk.len() == 1{
            Self::add_singleton_to_count(&mut self.library_id_generator);
        }
    }

    fn add_singleton_to_count(library_id_generator: &mut LibraryIdGenerator) {
        library_id_generator.get_duplicate_count_hist().increment1(1.0);
        library_id_generator.get_non_optical_duplicate_count_hist().increment1(1.0);
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
            if end != &best {
                self.add_index_as_duplicate(end.read1_index_in_file);
            }

            // in query-sorted case, these will be the same.
            // TODO: also in coordinate sorted, when one read is unmapped
            if end.read2_index_in_file != end.read1_index_in_file {
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
                library_id_generator.get_optical_duplicates_by_library_id_map(),
            );

            #[allow(non_snake_case)]
            let n_optical_dup_R = Self::track_optical_duplicates_with_histo(
                track_optical_duplicates_R.as_mut_slice(),
                keeper,
                optical_duplicate_finder,
                library_id_generator.get_optical_duplicates_by_library_id_map(),
            );

            n_optical_dup = n_optical_dup_F + n_optical_dup_R;
        } else {
            // No need to partition
            n_optical_dup = Self::track_optical_duplicates_with_histo(
                ends,
                keeper,
                optical_duplicate_finder,
                library_id_generator.get_optical_duplicates_by_library_id_map(),
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

#[cfg(test)]
mod test {
    use std::{borrow::Cow, path::Path};

    use clap::Parser;

    use crate::tests::{from_java_read_ends, JavaReadEndIterator, ToJsonFileForTest};

    use super::*;

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
}
