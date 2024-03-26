use std::{collections::HashSet, mem::size_of};

use rust_htslib::bam::{
    ext::BamRecordExtensions,
    record::{Aux, SAMReadGroupRecord},
    HeaderView, IndexedReader, Read, Record,
};

use crate::{
    cmdline::cli::Cli,
    hts::{duplicate_scoring_strategy::DuplicateScoringStrategy, SAMTag, SortOrder},
    markdup::{
        umi_utils::UmiUtil,
        utils::{
            library_id_generator::LibraryIdGenerator,
            read_ends::ReadEnds,
            read_ends_md_map::{
                DiskBasedReadEndsForMarkDuplicatesMap, MemoryBasedReadEndsForMarkDuplicatesMap,
                ReadEndsForMarkDuplicatesMap,
            },
            read_name_parser::ReadNameParserExt,
        },
    }, utils::{hash_code, CommonHasher},
};

use super::utils::{
    optical_duplicate_finder::OpticalDuplicateFinder, physical_location::PhysicalLocation,
    read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates,
    read_ends_for_mark_duplicates_with_barcodes::ReadEndsBarcodeData,
};

// type Error = Box<dyn std::error::Error + Send + Sync>;
use anyhow::{anyhow, Error};

pub(crate) struct MarkDuplicates {
    indexed_reader: IndexedReader,
    pg_ids_seen: HashSet<String>,
    library_id_generator: LibraryIdGenerator,
    optical_duplicate_finder: OpticalDuplicateFinder,
    cli: Cli,
}

impl MarkDuplicates {
    const NO_SUCH_INDEX: i64 = i64::MAX;

    fn build_sorted_read_end_vecs(&mut self, use_barcodes: bool) -> Result<(), Error> {
        let indexed_reader = &mut self.indexed_reader;

        let size_in_bytes = size_of::<ReadEndsForMarkDuplicates>();

        // TODO: MAX_RECORDS_IN_RAM~

        let mut duplicate_query_name = String::new();
        let mut duplicate_index = Self::NO_SUCH_INDEX;

        let assumed_sort_order = SortOrder::from_header(indexed_reader.header())?;

        indexed_reader.fetch(".")?;

        let mut record = Record::new();

        let tmp = match assumed_sort_order {
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

        while let Some(read_res) = indexed_reader.read(&mut record) {
            // Just unwrap `Result` to check the reading process successful.
            // If Err, it means that bam file may be corrupted. Don't need to recover this error.
            read_res.unwrap();

            // This doesn't have anything to do with building sorted ReadEnd lists, but it can be done in the same pass
            // over the input
            if self.cli.PROGRAM_RECORD_ID.is_some() {
                // Gather all PG IDs seen in merged input files in first pass.  These are gathered for two reasons:
                // - to know how many different PG records to create to represent this program invocation.
                // - to know what PG IDs are already used to avoid collisions when creating new ones.
                // Note that if there are one or more records that do not have a PG tag, then a null value
                // will be stored in this set.
                match record.aux(SAMTag::PG.name().as_bytes()) {
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
                && !record.qname().eq(duplicate_query_name.as_bytes())
            {
                duplicate_query_name.clear();
                duplicate_query_name.push_str(std::str::from_utf8(record.qname())?);

                duplicate_index = index;
            }

            if record.is_unmapped() {
                if record.tid() == -1 && matches!(assumed_sort_order, SortOrder::Coordinate) {
                    // When we hit the unmapped reads with no coordinate, no reason to continue (only in coordinate sort).
                    break;
                }

                // If this read is unmapped but sorted with the mapped reads, just skip it.
            } else if !(record.is_secondary() || record.is_supplementary()) {
                let index_for_read = match assumed_sort_order {
                    SortOrder::QueryName => duplicate_index,
                    _ => index,
                };
                // let fragment_end =
            }
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

        todo!()
    }

    fn get_read_one_barcode_value(&self, rec: &Record) -> i32 {
        todo!()
    }

    fn get_read_two_barcode_value(&self, rec: &Record) -> i32 {
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
