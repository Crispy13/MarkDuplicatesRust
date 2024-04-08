use serde::{Deserialize, Serialize};

/**
 * Little struct-like class to hold a record index, the index of the corresponding representative read, and duplicate set size information.
 */
#[derive(PartialEq, PartialOrd, Eq, Ord, Serialize, Deserialize)]
pub(crate) struct RepresentativeReadIndexer {
    read_index_in_file: i32,
    set_size: i32,
    representative_read_index_in_file: i32,
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
