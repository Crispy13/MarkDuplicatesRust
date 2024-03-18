use std::collections::HashSet;

use rust_htslib::bam::{HeaderView, IndexedReader, Read, Record};

use crate::{cmdline::cli::Cli, hts::SAMTag};

use super::utils::{physical_location::PhysicalLocation, read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates};

type Error = Box<dyn std::error::Error + Send + Sync>;
pub(crate) struct MarkDuplicates {
    indexed_reader: IndexedReader,
    pg_ids_seen: HashSet<String>,
    cli: Cli,
}

impl MarkDuplicates {
    const NO_SUCH_INDEX: i64 = i64::MAX;

    fn build_sorted_read_end_vecs(&mut self, use_barcodes: bool) -> Result<(), Error> {
        let indexed_reader = &mut self.indexed_reader;

        let duplicate_query_name = String::new();
        let duplicate_index = Self::NO_SUCH_INDEX;

        indexed_reader.fetch(".")?;

        let mut record = Record::new();
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
        }

        todo!()
    }
}

trait MarkDuplicatesExt {
    fn are_comparable_for_duplicates(
        lhs: &ReadEndsForMarkDuplicates,
        rhs: &ReadEndsForMarkDuplicates,
        compare_read2: bool,
        use_barcodes: bool,
    ) {
        let are_comparable = lhs.read_ends.library_id == rhs.read_ends.library_id;

        if use_barcodes && are_comparable { // areComparable is useful here to avoid the casts below

        }
    }
}

