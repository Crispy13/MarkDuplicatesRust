use serde::{Deserialize, Serialize};

/**
 * Little struct-like class to hold a record index, the index of the corresponding representative read, and duplicate set size information.
 */
#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
pub(crate) struct RepresentativeReadIndexer {
    pub(crate) read_index_in_file: i32,
    pub(crate) set_size: i32,
    pub(crate) representative_read_index_in_file: i32,
}


impl PartialEq for RepresentativeReadIndexer {
    fn eq(&self, other: &Self) -> bool {
        self.read_index_in_file == other.read_index_in_file
    }
}

impl PartialOrd for RepresentativeReadIndexer {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.read_index_in_file
            .partial_cmp(&other.read_index_in_file)
    }
}

impl Eq for RepresentativeReadIndexer {}
impl Ord for RepresentativeReadIndexer {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
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
