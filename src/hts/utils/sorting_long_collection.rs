use std::{path::PathBuf, sync::OnceLock};

pub(crate) struct SortingLongCollection {
    max_values_in_ram: usize,
    tmp_dir: Vec<PathBuf>,

    ram_values: Vec<i64>,
}



impl SortingLongCollection {
    pub(crate) const SIZEOF: usize = 8;
    
    #[allow(non_snake_case)]
    fn MAX_ITEMS_IN_RAM() -> usize {
        static N:OnceLock<usize> = OnceLock::new();
        *N.get_or_init(|| ((i32::MAX / 8) as f64 * 0.999).floor() as usize)
    }

    pub(crate) fn new(mut max_values_in_ram: usize, tmp_dir: Vec<PathBuf>) -> Self {
        if max_values_in_ram <= 0 {
            panic!("maxValuesInRam must be > 0");
        }

        max_values_in_ram = max_values_in_ram.min(Self::MAX_ITEMS_IN_RAM());
        Self {
            max_values_in_ram,
            tmp_dir,
            ram_values:Vec::with_capacity(max_values_in_ram)
        }
    }
}

impl Default for SortingLongCollection {
    fn default() -> Self {
        Self {
            max_values_in_ram: Default::default(),
            tmp_dir: Default::default(),
            ram_values: Default::default(),
        }
    }
}
