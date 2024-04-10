use std::{
    collections::BTreeSet,
    fs::File,
    io::{BufReader, Write},
    path::PathBuf,
    sync::OnceLock,
};

use tempfile::{NamedTempFile, TempPath};

use crate::hts::utils::save_as_byte_to_file;

pub(crate) struct SortingLongCollection {
    max_values_in_ram: usize,
    tmp_dir: Vec<PathBuf>,

    ram_values: Vec<i64>,
    num_values_in_ram: usize,

    files: Vec<TempPath>,

    done_adding: bool,
    cleaned_up: bool,
    // priority_queue: Option<BTreeSet<PeekFileRecordIterator<BufReader<File>, i64>>>,
}

impl SortingLongCollection {
    pub(crate) const SIZEOF: usize = 8;

    #[allow(non_snake_case)]
    fn MAX_ITEMS_IN_RAM() -> usize {
        static N: OnceLock<usize> = OnceLock::new();
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
            ram_values: Vec::with_capacity(max_values_in_ram),
            ..Default::default()
        }
    }

    pub(crate) fn add(&mut self, value: i64) {
        if self.done_adding {
            panic!("Cannot add after calling doneAddingStartIteration()");
        }

        if self.num_values_in_ram == self.max_values_in_ram {
            self.spill_to_disk();
        }

        self.ram_values.push(value);
        self.num_values_in_ram += 1;
    }

    /**
     * Sort the values in memory, write them to a file, and clear the buffer of values in memory.
     */
    fn spill_to_disk(&mut self) {
        self.ram_values
            .get_mut(..self.num_values_in_ram)
            .unwrap()
            .sort();

        let mut f = tempfile::Builder::new()
            .prefix("sortingcollection.")
            .suffix(".tmp")
            .tempfile_in(self.tmp_dir.first().unwrap())
            .unwrap(); // exit program when it fails to create temp file.

        let os = f.as_file_mut();

        save_as_byte_to_file(&self.num_values_in_ram, os).unwrap();

        for (_i, value) in (0..self.num_values_in_ram).zip(self.ram_values.drain(..)) {
            save_as_byte_to_file(&value, os).unwrap();
        }

        os.flush().unwrap();

        self.num_values_in_ram = 0;
        self.files.push(f.into_temp_path());
    }

    pub(crate) fn done_adding_start_iteration(&mut self) {
        if self.cleaned_up || self.done_adding {
            panic!("Cannot call doneAddingStartIteration() after cleanup() was called.");
        }

        self.done_adding = true;

        if self.files.is_empty() {
            self.ram_values
                .get_mut(..self.num_values_in_ram)
                .unwrap()
                .sort();
            return;
        }

        if self.num_values_in_ram > 0 {
            self.spill_to_disk();
        }

        // let mut priority_queue = BTreeSet::new();

        todo!()
    }
}

struct FileValueIterator {
    file: PathBuf,
    is: BufReader<File>,
    current_record: i64,
    is_current_record: bool,
}

impl FileValueIterator {
    fn new(
        file: PathBuf,
        is: BufReader<File>,
        current_record: i64,
        is_current_record: bool,
    ) -> Self {
        Self {
            file,
            is,
            current_record,
            is_current_record,
        }
    }
}

impl Default for SortingLongCollection {
    fn default() -> Self {
        Self {
            max_values_in_ram: Default::default(),
            tmp_dir: Default::default(),
            ram_values: Default::default(),
            num_values_in_ram: 0,
            files: Vec::new(),
            done_adding: false,
            cleaned_up: false,
            priority_queue: None,
        }
    }
}
