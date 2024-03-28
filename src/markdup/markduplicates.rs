use std::{collections::HashSet, fmt::Write, mem::size_of};

use rust_htslib::bam::{
    ext::BamRecordExtensions,
    record::{Aux, RecordExt, SAMReadGroupRecord},
    HeaderView, IndexedReader, Read, Record,
};

use crate::{
    cmdline::cli::Cli,
    hts::{duplicate_scoring_strategy::DuplicateScoringStrategy, SAMTag, SortOrder},
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

use super::utils::{
    optical_duplicate_finder::OpticalDuplicateFinder, physical_location::PhysicalLocation,
    read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates,
    read_ends_for_mark_duplicates_with_barcodes::ReadEndsBarcodeData,
};

// type Error = Box<dyn std::error::Error + Send + Sync>;
use anyhow::{anyhow, Error};

pub(crate) struct MarkDuplicates {
    pg_ids_seen: HashSet<String>,
    library_id_generator: LibraryIdGenerator,
    optical_duplicate_finder: OpticalDuplicateFinder,
    cli: Cli,

    frag_sort: Vec<ReadEndsForMarkDuplicates>,
}

impl MarkDuplicates {
    const NO_SUCH_INDEX: i64 = i64::MAX;

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

                        self.frag_sort.push(fragment_end);

                        tmp.put(
                            paired_ends.read_ends.read2_reference_index,
                            key,
                            paired_ends,
                        )?;
                    } else {
                        self.frag_sort.push(fragment_end);
                        
                        paired_ends = paired_ends_opt.unwrap();

                        let mates_ref_index = fragment_end.read_ends.read1_reference_index;
                        let mates_coordinate = fragment_end.read_ends.read1_coordinate;

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
                        self.frag_sort.push(paired_ends);
                    }
                }
            }

            // Print out some stats every 1m reads
            index+=1;
            

        }

        todo!()
    }

    /**
     * Builds a read ends object that represents a single read.
     */
    fn build_read_ends(
        &mut self,
        header: &HeaderView,
        index: i64,
        rec: &Record,
        use_barcode: bool,
    ) -> Result<ReadEndsForMarkDuplicates, Error> {
        let mut read_ends = ReadEnds::default();

        read_ends.read1_reference_index = rec.tid();
        read_ends.read1_coordinate = if rec.is_reverse() {
            rec.reference_end() as i32
        } else {
            rec.reference_start() as i32
        };

        read_ends.orientation = if rec.is_reverse() {
            ReadEnds::R
        } else {
            ReadEnds::F
        };

        let read1_index_in_file = index;
        let score = DuplicateScoringStrategy::compute_duplicate_score(
            rec,
            self.cli.DUPLICATE_SCORING_STRATEGY,
            false,
        )
        .unwrap(); // we can assure that no error occurs when `assume_mate_cigar=false`.

        // Doing this lets the ends object know that it's part of a pair
        if rec.is_paired() && !rec.is_mate_unmapped() {
            read_ends.read2_reference_index = rec.mtid()
        }

        // Fill in the library ID
        read_ends.library_id = self.library_id_generator.get_library_id(rec)?;

        // Fill in the location information for optical duplicates
        if self.optical_duplicate_finder.add_location_information(
            std::str::from_utf8(rec.qname())?.to_string(),
            &mut read_ends,
        ) {
            // calculate the RG number (nth in list)
            read_ends.read_group = 0;

            let rg = rec
                .aux(b"RG")
                .or_else(|err| Err(Error::from(err)))
                .and_then(|aux| match aux {
                    Aux::String(v) => Ok(v),
                    _ => Err(anyhow!("Invalid type for RG tag")),
                });
            let read_groups = header.get_read_groups();

            if let (Ok(rg), Ok(read_groups)) = (rg, read_groups) {
                for read_group in read_groups {
                    if read_group.get_read_group_id().eq(rg) {
                        break;
                    } else {
                        read_ends.read_group += 1;
                    }
                }
            }
        }

        let mut read_ends_for_md = ReadEndsForMarkDuplicates {
            score,
            read1_index_in_file,
            read_ends,
            ..Default::default()
        };

        if use_barcode {
            let top_strand_normalized_umi = UmiUtil::get_top_strand_normalized_umi(
                rec,
                &self.cli.BARCODE_TAG,
                self.cli.DUPLEX_UMI,
            )?;

            let barcode = hash_code(top_strand_normalized_umi);

            let mut barcode_data = ReadEndsBarcodeData {
                barcode,
                ..Default::default()
            };

            if (!rec.is_paired() || rec.is_first_in_template()) {
                barcode_data.read_one_barcode = self.get_read_one_barcode_value(rec);
            } else {
                barcode_data.read_two_barcode = self.get_read_two_barcode_value(rec);
            }

            read_ends_for_md.barcode_data = Some(barcode_data);
        }

        Ok(read_ends_for_md)
    }

    fn get_read_one_barcode_value(&self, rec: &Record) -> u64 {
        EstimateLibraryComplexity::get_read_barcode_value(rec, &self.cli.READ_ONE_BARCODE_TAG)
    }

    fn get_read_two_barcode_value(&self, rec: &Record) -> u64 {
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
        ).unwrap() // we can use unwrap() with assume_mate_cigar = false;
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
