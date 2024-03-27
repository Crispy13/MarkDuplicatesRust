use serde::{Deserialize, Serialize};

use super::physical_location::{
    impl_physical_location_core, PhysicalLocation, PhysicalLocationShort,
};

#[derive(Serialize, Deserialize, Clone)]
pub(crate) struct ReadEnds {
    pub(crate) library_id: i16,
    pub(crate) orientation: i8,
    pub(crate) read1_reference_index: i32,
    pub(crate) read1_coordinate: i32,
    pub(crate) read2_reference_index: i32,
    pub(crate) read2_coordinate: i32,

    // Additional information used to detect optical dupes
    pub(crate) read_group: i16,

    /** For optical duplicate detection the orientation matters regard to 1st or 2nd end of a mate */
    pub(crate) orientation_for_optical_duplicates: i8,

    /** A *transient* flag marking this read end as being an optical duplicate. */
    #[serde(skip)]
    pub(crate) is_optical_duplicate: bool,

    pub(crate) pls: PhysicalLocationShort,
}

impl Default for ReadEnds {
    /// > **Warning:**
    /// > Do not use this function directly.
    ///
    /// This method was implemented to designate only **some** fields' default values.
    fn default() -> Self {
        Self {
            read1_reference_index: -1,
            read1_coordinate: -1,
            read2_reference_index: -1,
            read2_coordinate: -1,
            read_group: -1,
            orientation_for_optical_duplicates: -1,
            is_optical_duplicate: false,

            library_id: Default::default(),
            orientation: Default::default(),
            pls: Default::default(),
        }
    }
}

impl ReadEnds {
    pub(crate) const F: i8 = 0;
    pub(crate) const R: i8 = 1;
    pub(crate) const FF: i8 = 2;
    pub(crate) const FR: i8 = 3;
    pub(crate) const RR: i8 = 4;
    pub(crate) const RF: i8 = 5;
}

impl PhysicalLocation for ReadEnds {
    impl_physical_location_core!(self, self.pls, i16);
    impl_physical_location_read_ends!(self, self);
    
}

macro_rules! impl_physical_location_read_ends {
    ($self:ident, $read_ends:expr) => {
        fn get_read_group(&$self) -> i16 {
            $read_ends.read_group
        }
    
        fn set_read_group(&mut $self, read_group: i16) {
            $read_ends.read_group = read_group;
        }
    
        fn get_library_id(&$self) -> i16 {
            $read_ends.library_id
        }
    
        fn set_library_id(&mut $self, library_id: i16) {
            $read_ends.library_id = library_id
        }
    };
}

pub(crate) use impl_physical_location_read_ends;

pub(crate) trait ReadEndsExt {
    fn is_paired(&self) -> bool;

    /**
     * Returns a single byte that encodes the orientation of the two reads in a pair.
     */
    fn get_orientation_bytes(read1_negative_strand: bool, read2_negative_strand: bool)
        -> i8;
}

impl ReadEndsExt for ReadEnds {
    impl_read_ends_ext!(self, self);
}

macro_rules! impl_read_ends_ext {
    ($self:ident, $read_end:expr) => {
        fn is_paired(&$self) -> bool {
            $read_end.read2_reference_index != -1
        }

        fn get_orientation_bytes(read1_negative_strand: bool, read2_negative_strand: bool) -> i8 {
            // maybe match would be better.
            if read1_negative_strand {
                if read2_negative_strand {
                    ReadEnds::RR
                } else {
                    ReadEnds::RF
                }
            } else {
                if read2_negative_strand {
                    ReadEnds::FR
                } else {
                    ReadEnds::FF
                }
            }
        }

    };
}

pub(crate) use impl_read_ends_ext;


