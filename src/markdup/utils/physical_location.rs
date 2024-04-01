use serde::{Deserialize, Serialize};

/**
 * Small interface that provides access to the physical location information about a cluster.
 * All values should be defaulted to -1 if unavailable.  ReadGroup and Tile should only allow
 * non-zero positive integers, x and y coordinates may be negative.
 */
pub(crate) trait PhysicalLocation: Serialize {
    const NO_VALUE: i32 = -1;

    fn get_read_group(&self) -> i16;

    fn set_read_group(&mut self, read_group: i16);

    fn get_tile(&self) -> i16;

    fn set_tile(&mut self, tile: i16);

    fn get_x(&self) -> i32;

    fn set_x(&mut self, x: i32);

    fn get_y(&self) -> i32;

    fn set_y(&mut self, y: i32);

    fn get_library_id(&self) -> i16;

    fn set_library_id(&mut self, library_id: i16);

    /** Default implementation of a method to check whether real location data has been set. */
    fn has_location(&self) -> bool {
        self.get_tile() as i32 != Self::NO_VALUE
    }

    // fn eq(&self, other: &Self) -> bool
    // where
    //     Self: Sized + PartialEq,
    // {
    //     self == other
    // }
}

/**
 * Small class that provides access to the physical location information about a cluster.
 * All values should be defaulted to -1 if unavailable.  Tile should only allow
 * non-zero positive integers, x and y coordinates must be non-negative.
 * This is different from PhysicalLocationShort in that the x and y positions are ints, not shorts
 * thus, they do not overflow within a HiSeqX tile.
 */
#[derive(Serialize, Deserialize)]
pub(crate) struct PhysicalLocationInt {
    pub(crate) tile: i16,
    pub(crate) x: i32,
    pub(crate) y: i32,
}

macro_rules! impl_default_for_physical_loc {
    ($target:ident) => {
        impl Default for $target {
            fn default() -> Self {
                Self {
                    tile: -1,
                    x: -1,
                    y: -1,
                }
            }
        }
    };
}

impl_default_for_physical_loc!(PhysicalLocationInt);

macro_rules! impl_physical_location_unimpl {
    ($self:ident, $pl:expr, $xyt:ty) => {
        fn get_read_group(&$self) -> i16 {
            unimplemented!()
        }

        fn set_read_group(&mut $self, read_group:i16) {
            unimplemented!()
        }

        fn get_library_id(&$self) -> i16 {
            unimplemented!()
        }

        fn set_library_id(&mut $self, library_id: i16) {
            unimplemented!()
        }
    }
}

macro_rules! impl_physical_location_core {
    ($self:ident, $pl:expr, $xyt:ty) => {
        fn get_tile(&$self) -> i16 {
            $pl.tile
        }

        fn set_tile(&mut $self, tile: i16) {
            $pl.tile = tile;
        }

        fn get_x(&$self) -> i32 {
            $pl.x
        }

        fn set_x(&mut $self, x: i32) {
            $pl.x = x;
        }

        fn get_y(&$self) -> i32 {
            $pl.y
        }

        fn set_y(&mut $self, y: i32) {
            $pl.y = y;
        }


    };
}

pub(crate) use impl_physical_location_core;

impl PhysicalLocation for PhysicalLocationInt {
    impl_physical_location_core!(self, self, i32);
    impl_physical_location_unimpl!(self, self, i32);
}

#[derive(Serialize, Deserialize, Clone, PartialEq)]
pub(crate) struct PhysicalLocationShort {
    pub(crate) tile: i16,
    pub(crate) x: i32,
    pub(crate) y: i32,
}

impl_default_for_physical_loc!(PhysicalLocationShort);

impl PhysicalLocationShort {
    fn set_x(&mut self, x: i32) {
        self.x = x as i16 as i32;
    }

    fn set_y(&mut self, y: i32) {
        self.y = y as i16 as i32;
    }
}

impl PhysicalLocation for PhysicalLocationShort {
    impl_physical_location_core!(self, self, i16);
    impl_physical_location_unimpl!(self, self, i16);
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn ovr() {
        let a = PhysicalLocationInt {
            tile: todo!(),
            x: todo!(),
            y: todo!(),
        };

        // a.get_read_group();
    }
}
