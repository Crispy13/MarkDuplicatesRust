use std::{
    collections::BTreeSet,
    fs::File,
    io::{BufReader, Write},
    path::PathBuf,
    sync::OnceLock,
};

use tempfile::{NamedTempFile, TempPath};

use crate::hts::utils::{
    load_byte_file_as_obj, save_as_byte_to_file, sorting_collection::MergingIterator,
};

pub(crate) struct SortingLongCollection {
    max_values_in_ram: usize,
    tmp_dir: Vec<PathBuf>,

    ram_values: Vec<i64>,
    num_values_in_ram: usize,

    files: Vec<TempPath>,

    done_adding: bool,
    cleaned_up: bool,
    priority_queue: Option<MergingIterator<BufReader<File>, i64>>,

    iteration_index: usize,
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

    pub(crate) fn done_adding(&mut self) {
        if self.cleaned_up || self.done_adding {
            panic!("Cannot call doneAddingStartIteration() after cleanup() was called.");
        }

        self.done_adding = true;
    }

    pub(crate) fn drain(&mut self) -> SortingLongCollectionDrain<'_> {
        if self.files.is_empty() {
            self.ram_values
                .get_mut(..self.num_values_in_ram)
                .unwrap()
                .sort();
            let r = SortingLongCollectionDrain::InMemoryDrain(
                self.ram_values.drain(..self.num_values_in_ram),
            );

            self.num_values_in_ram = 0;

            r
        } else {
            if self.num_values_in_ram > 0 {
                self.spill_to_disk();
            }

            self.ram_values.clear();

            SortingLongCollectionDrain::MergingIterator(MergingIteratorLong::new(&self.files))
        }
    }
    
    pub(crate) fn clean_up(&mut self) {
        self.done_adding = true;
        self.cleaned_up = true;
        self.ram_values.clear();


    }
}

pub(crate) enum SortingLongCollectionDrain<'a> {
    InMemoryDrain(std::vec::Drain<'a, i64>),
    MergingIterator(MergingIteratorLong),
}

impl<'a> Iterator for SortingLongCollectionDrain<'a> {
    type Item=i64;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            SortingLongCollectionDrain::InMemoryDrain(it) => it.next(),
            SortingLongCollectionDrain::MergingIterator(it) => it.next(),
        }
    }
}

pub(crate) struct MergingIteratorLong {
    priority_queue: BTreeSet<PeekFileValueIterator>,
}

impl MergingIteratorLong {
    fn new(files: &[TempPath]) -> Self {
        let mut priority_queue = BTreeSet::new();
        for f in files.iter() {
            let it = PeekFileValueIterator::new(FileValueIterator::new(f.to_path_buf()));
            if it.peeked().is_some() {
                priority_queue.insert(it);
            }
        }
        Self { priority_queue }
    }
}

impl Iterator for MergingIteratorLong {
    type Item = i64;

    fn next(&mut self) -> Option<Self::Item> {
        match self.priority_queue.pop_first() {
            Some(mut peek_iter) => {
                let r = peek_iter.next();
                if peek_iter.peeked().is_some() {
                    self.priority_queue.insert(peek_iter);
                }
                r
            }
            None => None,
        }
    }
}

struct FileValueIterator {
    file: PathBuf,
    is: BufReader<File>,
    current_record: i64,
    is_current_record: bool,

    n_remained: usize,
}

impl FileValueIterator {
    fn new(file: PathBuf) -> Self {
        let mut is = BufReader::new(File::open(&file).unwrap());

        let n_remained = load_byte_file_as_obj::<usize>(&mut is).unwrap();

        let mut s = Self {
            file,
            is,
            current_record: 0,
            is_current_record: true,
            n_remained,
        };

        s.next().unwrap();

        s
    }
}

impl Iterator for FileValueIterator {
    type Item = i64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.n_remained < 1 {
            self.current_record = 0;
            self.is_current_record = false;

            return None;
        }

        let ret = self.current_record;

        self.current_record = load_byte_file_as_obj(&mut self.is).unwrap();

        Some(ret)
    }
}

struct PeekFileValueIterator {
    iter: FileValueIterator,
    peeked: Option<Option<i64>>,
}

impl PeekFileValueIterator {
    fn new(iter: FileValueIterator) -> Self {
        // assert!(iter.is_current_record);

        let mut s = Self { iter, peeked: None };

        s.peek().unwrap(); // assume `iter` has at least one item.

        s
    }

    fn peek(&mut self) -> Option<i64> {
        let iter = &mut self.iter;
        self.peeked.get_or_insert_with(|| iter.next()).clone()
    }

    /// get peeked item.
    ///
    /// It's the first item of the file in this iterator.
    pub(crate) fn peeked(&self) -> Option<&i64> {
        self.peeked.as_ref().unwrap().as_ref()
    }
}

impl Iterator for PeekFileValueIterator {
    type Item = i64;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let r = match self.peeked.take() {
            Some(v) => v,
            None => self.iter.next(),
        };

        self.peek(); // assure that self always has Some() in peeked field.

        r
    }
}

impl PartialOrd for PeekFileValueIterator {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.peeked()
            .unwrap()
            .partial_cmp(other.peeked().as_ref().unwrap())
    }
}

impl PartialEq for PeekFileValueIterator {
    fn eq(&self, other: &Self) -> bool {
        self.peeked().unwrap().eq(other.peeked().unwrap())
    }
}
impl Eq for PeekFileValueIterator {}

impl Ord for PeekFileValueIterator {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
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
            iteration_index: 0,
        }
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn save_and_load() {}
}
