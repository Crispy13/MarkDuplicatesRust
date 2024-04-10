use serde::{Deserialize, Serialize};

/**
 * Little struct-like class to hold a record index, the index of the corresponding representative read, and duplicate set size information.
 */
#[derive(PartialEq, PartialOrd, Eq, Ord, Serialize, Deserialize)]
pub(crate) struct RepresentativeReadIndexer {
    pub(crate) read_index_in_file: i32,
    pub(crate) set_size: i32,
    pub(crate) representative_read_index_in_file: i32,
}

impl Default for RepresentativeReadIndexer {
    fn default() -> Self {
        Self {
            read_index_in_file: -1,
            set_size: -1,
            representative_read_index_in_file: -1,
        }
    }
}
