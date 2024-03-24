use std::collections::HashMap;

use rust_htslib::bam::{record::RecordExt, HeaderView, Record};

// type Error = Box<dyn std::error::Error + Send + Sync>;
use anyhow::Error;
pub(crate) struct LibraryIdGenerator {
    library_ids: HashMap<String, i16>, // from library string to library id
    next_library_id: i16,
}

impl LibraryIdGenerator {
    const UNKNOWN_LIBRARY: &'static str = "Unknown Library";

    /** Get the library ID for the given SAM record. */
    pub(crate) fn get_library_id(&mut self, rec: &Record) -> Result<i16, Error> {
        let library = Self::get_library_name(rec)?;

        let library_id = match self.library_ids.get(library) {
            Some(v) => v,
            None => {
                self.library_ids
                    .insert(library.to_string(), self.next_library_id);
                self.next_library_id += 1;
                self.library_ids.get(library).unwrap()
            }
        }
        .clone();

        Ok(library_id)
    }

    /**
     * Gets the library name from the header for the record. If the RG tag is not present on
     * the record, or the library isn't denoted on the read group, a constant string is
     * returned.
     */
    pub(crate) fn get_library_name<'r>(rec: &'r Record) -> Result<&'r str, Error> {
        match rec.get_read_group()?.get_library() {
            Some(lb) => Ok(lb),
            None => Ok(Self::UNKNOWN_LIBRARY),
        }
    }
}

#[cfg(example)]
mod examples {
    use super::*;

    struct A {
        string: String,
    }

    struct B<'s> {
        str_slice: &'s str,
    }

    impl<'s> B<'s> {
        fn from_a(a: &'s A) -> Self {
            Self {
                str_slice: &a.string,
            }
        }

        fn get_str(&self) -> &'s str {
            self.str_slice
        }
    }

    fn get_str_slice<'a>(a: &'a A, s: &str) -> &'a str {
        let b = B::from_a(a);
        // let b = B {
        //     str_slice: &a.string,
        // };

        b.get_str()
    }
}
