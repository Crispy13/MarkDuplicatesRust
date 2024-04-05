use core::fmt;
use std::cmp::Ordering;

use rust_htslib::bam::{self, ext::BamRecordExtensions, HeaderView, Record};
use serde::{Deserialize, Serialize};

use super::{
    physical_location::{impl_physical_location_core, PhysicalLocation, PhysicalLocationShort},
    read_ends::{impl_physical_location_read_ends, impl_read_ends_ext, ReadEnds, ReadEndsExt},
    read_ends_for_mark_duplicates_with_barcodes::ReadEndsBarcodeData,
};

/**
 * Little struct-like class to hold read pair (and fragment) end data for MarkDuplicatesWithMateCigar
 *
 * @author Nils Homer
 */
#[derive(Serialize, Deserialize, Clone, Debug)]
pub(crate) struct ReadEndsForMarkDuplicates {
    pub(crate) score: i16,
    pub(crate) read1_index_in_file: i64,
    pub(crate) read2_index_in_file: i64,
    pub(crate) duplicate_set_size: i32,

    pub(crate) read_ends: ReadEnds,
    pub(crate) barcode_data: Option<ReadEndsBarcodeData>,
}

impl ReadEndsExt for ReadEndsForMarkDuplicates {
    impl_read_ends_ext!(self, self.read_ends);
}

impl PhysicalLocation for ReadEndsForMarkDuplicates {
    impl_physical_location_core!(self, self.read_ends.pls, i16);
    impl_physical_location_read_ends!(self, self.read_ends);
}

impl ReadEndsForMarkDuplicates {
    /*
    What do we need to store you ask?  Well, we need to store:
       - byte: orientation
       - short: libraryId, readGroup, tile, x, y, score
       - int: read1ReferenceIndex, read1Coordinate, read2ReferenceIndex, read2Coordinate, duplicateSetSize
       - long: read1IndexInFile, read2IndexInFile
     */
    const SIZE_OF: i32 = (1 * 1) + (5 * 2) + (5 * 4) + (8 * 2) + 1
        + 8 // last 8 == reference overhead
        + 13; // This is determined experimentally with JProfiler

    // pub(crate) fn new() -> Self {
    //     Self::default()
    // }

    pub(crate) fn with_read(read: Self) -> Self {
        let physical_loc_short = PhysicalLocationShort {
            tile: read.get_tile(),
            x: read.get_x(),
            y: read.get_y(),
        };

        let read_ends = ReadEnds {
            library_id: read.get_library_id(),
            orientation: read.read_ends.orientation,
            read1_reference_index: read.read_ends.read1_reference_index,
            read1_coordinate: read.read_ends.read1_coordinate,
            read2_reference_index: read.read_ends.read2_reference_index,
            read2_coordinate: read.read_ends.read2_coordinate,
            read_group: read.get_read_group(),
            orientation_for_optical_duplicates: read.read_ends.orientation_for_optical_duplicates,
            is_optical_duplicate: Default::default(),
            pls: physical_loc_short,
        };

        Self {
            score: read.score,
            read1_index_in_file: read.read1_index_in_file,
            read2_index_in_file: read.read2_index_in_file,
            duplicate_set_size: Default::default(),
            read_ends,
            barcode_data: None,
        }
    }
}

impl fmt::Display for ReadEndsForMarkDuplicates {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} {} {}",
            self.read1_index_in_file, self.read_ends.read1_coordinate, self.score
        )
    }
}

impl Default for ReadEndsForMarkDuplicates {
    fn default() -> Self {
        Self {
            score: 0,
            read1_index_in_file: -1,
            read2_index_in_file: -1,
            duplicate_set_size: -1,
            read_ends: Default::default(),
            barcode_data: Default::default(),
        }
    }
}

impl PartialOrd for ReadEndsForMarkDuplicates {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self
            .read_ends
            .library_id
            .partial_cmp(&other.read_ends.library_id)
        {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        match (&self.barcode_data, &self.barcode_data) {
            (Some(lb), Some(rb)) => {
                match lb.barcode.partial_cmp(&rb.barcode) {
                    Some(Ordering::Equal) => {}
                    ord => return ord,
                }

                match lb.read_one_barcode.partial_cmp(&rb.read_one_barcode) {
                    Some(Ordering::Equal) => {}
                    ord => return ord,
                }

                match lb.read_two_barcode.partial_cmp(&rb.read_two_barcode) {
                    Some(Ordering::Equal) => {}
                    ord => return ord,
                }
            }
            (None, None) => {}
            _ => {
                panic!(
                    "lhs.is_some()={}, rhs.is_some()={}",
                    self.barcode_data.is_some(),
                    other.barcode_data.is_some()
                );
            }
        }

        match self
            .read_ends
            .read1_reference_index
            .partial_cmp(&other.read_ends.read1_reference_index)
        {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        match self
            .read_ends
            .read1_coordinate
            .partial_cmp(&other.read_ends.read1_coordinate)
        {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        match self
            .read_ends
            .orientation
            .partial_cmp(&other.read_ends.orientation)
        {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        match self
            .read_ends
            .read2_reference_index
            .partial_cmp(&other.read_ends.read2_reference_index)
        {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        match self
            .read_ends
            .read2_coordinate
            .partial_cmp(&other.read_ends.read2_coordinate)
        {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        match self.get_tile().partial_cmp(&other.get_tile()) {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        match self.get_x().partial_cmp(&other.get_x()) {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        match self.get_y().partial_cmp(&other.get_y()) {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        // The following is arbitrary and is only included for completeness.
        // Other implementations may chose to forgo this tiebreak if they do not have
        // access to the index-in-file of the records (e.g. SPARK implmentations)

        match self
            .read1_index_in_file
            .partial_cmp(&other.read1_index_in_file)
        {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        self.read2_index_in_file
            .partial_cmp(&other.read2_index_in_file)
    }
}

impl PartialEq for ReadEndsForMarkDuplicates {
    fn eq(&self, other: &Self) -> bool {
        match self.partial_cmp(other) {
            Some(Ordering::Equal) => true,
            _ => false,
        }
    }
}

impl Eq for ReadEndsForMarkDuplicates {}
impl Ord for ReadEndsForMarkDuplicates {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or_else(|| panic!("impl Ord failed for `ReadEndsForMarkDuplicates`."))
    }
}