use serde::{Deserialize, Serialize};

use super::{physical_location::PhysicalLocation, read_name_parser::{impl_read_name_parser_ext, ReadNameParser, ReadNameParserExt}};
use regex::Regex;

type Error = Box<dyn std::error::Error + Send + Sync>;

#[derive(Serialize, Deserialize)]
pub(crate) struct OpticalDuplicateFinder {
    pub(crate) optical_duplicate_pixel_distance: i32,
    big_duplicate_set_size: i32,
    max_duplicate_set_size: i64,

    rnp: ReadNameParser,
}

impl Default for OpticalDuplicateFinder {
    fn default() -> Self {
        Self {
            optical_duplicate_pixel_distance: Self::DEFAULT_OPTICAL_DUPLICATE_DISTANCE,
            big_duplicate_set_size: Self::DEFAULT_BIG_DUPLICATE_SET_SIZE,
            max_duplicate_set_size: Self::DEFAULT_MAX_DUPLICATE_SET_SIZE,
            rnp: ReadNameParser::new(),
        }
    }
}

impl OpticalDuplicateFinder {
    pub(crate) const DEFAULT_OPTICAL_DUPLICATE_DISTANCE: i32 = 100;
    pub(crate) const DEFAULT_MAX_DUPLICATE_SET_SIZE: i64 = 300000;
    pub(crate) const DEFAULT_BIG_DUPLICATE_SET_SIZE: i32 = 1000;
}

impl ReadNameParserExt for OpticalDuplicateFinder {
    impl_read_name_parser_ext!(self, self.rnp);
}

