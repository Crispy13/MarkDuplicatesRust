use serde::{Deserialize, Serialize};

use super::{
    physical_location::{impl_physical_location_core, PhysicalLocation},
    read_ends::{impl_physical_location_read_ends, impl_read_ends_ext, ReadEnds, ReadEndsExt},
    read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates,
};

#[derive(Serialize, Deserialize, Clone, Debug)]
pub(crate) struct ReadEndsForMarkDuplicatesWithBarcodes {
    pub(crate) barcode: i32,          // primary barcode for this read (and pair)
    pub(crate) read_one_barcode: i32, // read one barcode, 0 if not present
    pub(crate) read_two_barcode: i32, // read two barcode, 0 if not present or not paired

    // pub(crate) read_ends_for_mark_duplicates: ReadEndsForMarkDuplicates,
}

impl ReadEndsForMarkDuplicatesWithBarcodes {
    const BARCODE_DEFAULT:i32 = 0;
    const READ_ONE_BARCODE_DEFAULT:i32 = 0;
    const READ_TWO_BARCODE_DEFAULT:i32 = 0;
}



// impl ReadEndsExt for ReadEndsForMarkDuplicatesWithBarcodes {
//     impl_read_ends_ext!(self, self.read_ends_for_mark_duplicates.read_ends);
// }

// impl PhysicalLocation for ReadEndsForMarkDuplicatesWithBarcodes {
//     impl_physical_location_core!(self, self.read_ends_for_mark_duplicates.read_ends.pls, i16);
//     impl_physical_location_read_ends!(self, self.read_ends_for_mark_duplicates.read_ends);
// }

// impl From<ReadEndsForMarkDuplicates> for ReadEndsForMarkDuplicatesWithBarcodes {
//     fn from(value: ReadEndsForMarkDuplicates) -> Self {
//         Self {
//             barcode: Self::BARCODE_DEFAULT,
//             read_one_barcode: Self::READ_ONE_BARCODE_DEFAULT,
//             read_two_barcode: Self::READ_TWO_BARCODE_DEFAULT,
//             read_ends_for_mark_duplicates: ReadEndsForMarkDuplicates::with_read(value),
//         }

//     }
// }


