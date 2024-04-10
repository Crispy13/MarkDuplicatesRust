use std::{
    borrow::{Borrow, BorrowMut, Cow},
    collections::{BTreeSet, HashSet},
    fmt,
    fs::File,
    io::{BufReader, BufWriter, Read},
    iter::Peekable,
    marker::PhantomData,
    mem::size_of,
    ops::{Deref, DerefMut},
    path::PathBuf,
};

use anyhow::Error;
use serde::{Deserialize, Serialize};
use tempfile::{tempfile, NamedTempFile, TempPath};

use crate::{
    hts::utils::save_as_byte_to_file,
    markdup::utils::read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates,
    utils::{get_allocated_and_resident_mem_of_app, human_readable_byte_count, mem_stats},
};

use super::load_byte_file_as_obj;

/**
 * Collection to which many records can be added.  After all records are added, the collection can be
 * iterated, and the records will be returned in order defined by the comparator.  Records may be spilled
 * to a temporary directory if there are more records added than will fit in memory.  As a result of this,
 * the objects returned may not be identical to the objects added to the collection, but they should be
 * equal as determined by the codec used to write them to disk and read them back.
 * <p>
 * When iterating over the collection, the number of file handles required is numRecordsInCollection/maxRecordsInRam.
 * If this becomes a limiting factor, a file handle cache could be added.
 * <p>
 * If Snappy DLL is available and snappy.disable system property is not set to true, then Snappy is used
 * to compress temporary files.
 */
#[derive(Default)]
pub(crate) struct SortingCollection<T> {
    /**
     * Directories where files of sorted records go.
     */
    tmp_dirs: Vec<PathBuf>,

    // /**
    //  * Used to write records to file, and used as a prototype to create codecs for reading.
    //  */
    // codec: Codec,
    ram_records: Vec<T>,
    done_adding: bool,
    iteration_started: bool,

    num_records_in_ram: usize,
    max_records_in_ram: usize,

    print_record_size_sampling: bool,

    files: Vec<TempPath>,

    cleaned_up: bool,
}

impl<T> SortingCollection<T>
where
    T: Ord + Serialize + for<'de> Deserialize<'de>,
{
    pub(crate) fn new(
        max_records_in_ram: usize,
        print_record_size_sampling: bool,
        tmp_dirs: Vec<PathBuf>,
    ) -> Self {
        if max_records_in_ram <= 0 {
            panic!("maxRecordsInRam must be > 0");
        }

        if tmp_dirs.is_empty() {
            panic!("At least one temp directory must be provided.");
        }

        let ram_records = Vec::with_capacity(max_records_in_ram);

        Self {
            tmp_dirs,
            ram_records,
            done_adding: false,
            iteration_started: false,
            num_records_in_ram: 0,
            max_records_in_ram,
            print_record_size_sampling,
            files: Vec::new(),
            cleaned_up: false,
        }
    }

    #[allow(non_upper_case_globals)]
    const log: &'static str = stringify!(SortingCollection);

    pub(crate) fn add(&mut self, rec: T) -> Result<(), Error> {
        if self.done_adding {
            panic!("Cannot add after calling done_adding()");
        }

        if self.iteration_started {
            panic!("Cannot add after calling iterator()")
        }

        if self.num_records_in_ram == self.max_records_in_ram {
            let mut start_mem = 0;
            if self.print_record_size_sampling {
                start_mem = mem_stats().record().resident;
            }

            self.spill_to_disk()?;

            if self.print_record_size_sampling {
                let end_mem = mem_stats().record().resident;

                let used_bytes = end_mem - start_mem;

                log::debug!(target: Self::log, "{} records in ram required approximately {} memory or {} per record. ",
                    self.max_records_in_ram,
                    human_readable_byte_count(used_bytes),
                    human_readable_byte_count(used_bytes / self.max_records_in_ram));
            }
        }

        self.ram_records.push(rec);
        self.num_records_in_ram += 1;

        Ok(())
    }

    ///
    /// Sort the records in memory, write them to a file, and clear the buffer of records in memory.
    ///
    /// # Errors
    /// - failed to create temp file.
    /// - failed to save as byte to file
    ///
    ///
    fn spill_to_disk(&mut self) -> Result<(), Error> {
        self.ram_records.get_mut(0..self.num_records_in_ram).unwrap_or_else(|| {
            panic!("Code error. Failed to get slice of `ram_records` with index 0..self.num_records_in_ram.");
        }).sort();

        let mut f = tempfile::Builder::new()
            .prefix("sortingcollection.")
            .suffix(".tmp")
            .tempfile_in(self.tmp_dirs.first().unwrap())?;

        let mut buf_writer = BufWriter::new(f.as_file_mut());

        // Not in original java code.
        // write number of records to the front of the file.
        save_as_byte_to_file(&self.num_records_in_ram, &mut buf_writer)?;

        self.ram_records
            .drain(..self.num_records_in_ram)
            .map(|r| save_as_byte_to_file(&r, &mut buf_writer))
            .collect::<Result<(), Error>>()?;

        drop(buf_writer);

        self.num_records_in_ram = 0;

        self.files.push(f.into_temp_path());

        Ok(())
    }

    /**
     * This method can be called after caller is done adding to collection, in order to possibly free
     * up memory.  If iterator() is called immediately after caller is done adding, this is not necessary,
     * because iterator() triggers the same freeing.
     */
    pub(crate) fn done_adding(&mut self) -> Result<(), Error> {
        if self.cleaned_up {
            panic!("Cannot call done_adding() after cleanup() was called.");
        }

        if self.done_adding {
            return Ok(());
        }

        self.done_adding = true;

        if self.files.is_empty() {
            return Ok(());
        }

        if self.num_records_in_ram > 0 {
            self.spill_to_disk()?
        }

        Ok(())
    }

    /**
     * Prepare to iterate through the records in order.  This method may be called more than once,
     * but add() may not be called after this method has been called.
     */
    pub(crate) fn iter(&mut self) -> Result<SortingCollectionIter<T>, Error> {
        self.get_ready_for_iter()?;

        if self.files.is_empty() {
            Ok(SortingCollectionIter::InMemoryIter(self.in_memory_iter()))
        } else {
            Ok(SortingCollectionIter::MergingIterator(
                MergingIterator::new(&self.files),
            ))
        }
    }

    fn get_ready_for_iter(&mut self) -> Result<(), Error> {
        if self.cleaned_up {
            panic!("Cannot call iterator() after cleanup() was called.");
        }

        self.done_adding()?;

        self.iteration_started = true;

        Ok(())
    }

    pub(crate) fn drain(&mut self) -> Result<SortingCollectionDrain<T>, Error> {
        self.get_ready_for_iter()?;

        if self.files.is_empty() {
            Ok(SortingCollectionDrain::InMemoryDrain(self.in_memory_drain()))
        } else {
            Ok(SortingCollectionDrain::MergingIterator(
                MergingIterator::new(&self.files),
            ))
        }

    }

    fn sort_current(&mut self) {
        self.ram_records.get_mut(0..self.num_records_in_ram).unwrap_or_else(|| {
            panic!("Code error. Failed to get slice of `ram_records` with index 0..self.num_records_in_ram.");
        }).sort();
    }

    fn in_memory_drain(&mut self) -> std::vec::Drain<'_, T> {
        self.sort_current();

        self.ram_records.drain(..self.num_records_in_ram)
    }

    fn in_memory_iter(&mut self) -> std::slice::Iter<'_, T> {
        self.sort_current();

        // self.ram_records.drain(..self.num_records_in_ram)
        self.ram_records.get(..self.num_records_in_ram).unwrap().iter()
    }
    
    pub(crate) fn clean_up(&mut self) {
        self.iteration_started = true;
        self.cleaned_up = true;

    }
}

pub(crate) enum SortingCollectionIter<'a, T>
where
    T: for<'de> Deserialize<'de>,
{
    InMemoryIter(std::slice::Iter<'a, T>),
    MergingIterator(MergingIterator<BufReader<File>, T>),
}


impl<'a, T> Iterator for SortingCollectionIter<'a, T>
where
    T: for<'de> Deserialize<'de> + ToOwned<Owned = T>,
    PeekFileRecordIterator<BufReader<File>, T>: Ord,
{
    type Item = CowForSC<'a, T>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            SortingCollectionIter::InMemoryIter(it) => it.next().map(CowForSC::Borrowed),
            SortingCollectionIter::MergingIterator(it) => it.next().map(CowForSC::Owned),
        }
    }
}

pub(crate) enum SortingCollectionDrain<'a, T>
where
    T: for<'de> Deserialize<'de>,
{
    InMemoryDrain(std::vec::Drain<'a, T>),
    MergingIterator(MergingIterator<BufReader<File>, T>),
}


impl<'a, T> Iterator for SortingCollectionDrain<'a, T>
where
    T: for<'de> Deserialize<'de> + ToOwned<Owned = T>,
    PeekFileRecordIterator<BufReader<File>, T>: Ord,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            SortingCollectionDrain::InMemoryDrain(it) => it.next(),
            SortingCollectionDrain::MergingIterator(it) => it.next(),
        }
    }
}

struct FileRecordIterator<R, T> {
    inner: R,
    buf: Vec<u8>,

    n_remaining_record: usize,
    phantom_data: PhantomData<T>,
}

impl<R: Read, T> FileRecordIterator<R, T> {
    fn new(mut inner: R) -> Self {
        let n_remaining_record = Self::get_saved_count(&mut inner);

        Self {
            inner,
            buf: Vec::with_capacity(size_of::<T>()),
            phantom_data: PhantomData,
            n_remaining_record,
        }
    }

    fn get_saved_count(inner: &mut R) -> usize {
        // The file's first bytes indicates number of records.
        // if the following fails, it means this software fails.
        bincode::deserialize_from::<_, usize>(inner).unwrap()
    }
}

impl<T: for<'de> Deserialize<'de>, R: Read> Iterator for FileRecordIterator<R, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        // let buf = &mut self.buf;

        // buf.clear();

        // match self.inner.read_exact(buf) {
        //     Ok(_) => {}
        //     Err(err) => match err.kind() {
        //         std::io::ErrorKind::UnexpectedEof => return None,
        //         err => panic!("{:?}", err),
        //     },
        // }

        // match bincode::deserialize(&buf) {
        //     Ok(v) => Some(v),
        //     Err(err) => panic!("{:?}", err),
        // }
        if self.n_remaining_record > 0 {
            self.n_remaining_record -= 1;
            Some(bincode::deserialize_from::<_, T>(&mut self.inner).unwrap())
        } else {
            None
        }
    }
}

pub(crate) struct PeekFileRecordIterator<R, T>
where
    T: for<'de> Deserialize<'de>,
    R: Read,
{
    iter: FileRecordIterator<R, T>,
    n: i32,
    peeked: Option<Option<T>>,
    // n_remaining_record: usize,
}

impl<R, T> Iterator for PeekFileRecordIterator<R, T>
where
    T: for<'de> Deserialize<'de>,
    R: Read,
{
    type Item = T;

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

impl<R, T> PeekFileRecordIterator<R, T>
where
    T: for<'de> Deserialize<'de>,
    R: Read,
{
    fn new(file_iter: FileRecordIterator<R, T>, n: i32) -> Self {
        // let n_remaining_record = Self::get_saved_count(&mut file_iter);

        let mut s = Self {
            iter: file_iter,
            n,
            peeked: None,
            // n_remaining_record
        };

        s.peek();

        s
    }

    pub(crate) fn peek(&mut self) -> Option<&T> {
        let iter = &mut self.iter;
        self.peeked.get_or_insert_with(|| iter.next()).as_ref()
    }

    /// get peeked item.
    ///
    /// It's the first item of the file in this iterator.
    pub(crate) fn peeked(&self) -> Option<&T> {
        self.peeked.as_ref().unwrap().as_ref()
    }
}

pub(crate) struct MergingIterator<R, T>
where
    T: for<'de> Deserialize<'de>,
    R: Read,
{
    pub(crate) queue: BTreeSet<PeekFileRecordIterator<R, T>>,
}

impl<R, T> Ord for PeekFileRecordIterator<R, T>
where
    T: for<'de> Deserialize<'de> + PartialEq + PartialOrd,
    R: Read,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl<R, T> Eq for PeekFileRecordIterator<R, T>
where
    T: for<'de> Deserialize<'de> + PartialEq + PartialOrd,
    R: Read,
{
}

impl<R, T> PartialOrd for PeekFileRecordIterator<R, T>
where
    T: for<'de> Deserialize<'de> + PartialEq + PartialOrd,
    R: Read,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.peeked().partial_cmp(&other.peeked()) {
            Some(std::cmp::Ordering::Equal) => {}
            ord => return ord,
        }

        self.n.partial_cmp(&other.n)
    }
}

impl<R, T> PartialEq for PeekFileRecordIterator<R, T>
where
    T: for<'de> Deserialize<'de> + PartialEq,
    R: Read,
{
    fn eq(&self, other: &Self) -> bool {
        self.peeked() == other.peeked()
    }
}

impl<T> MergingIterator<BufReader<File>, T>
where
    T: for<'de> Deserialize<'de> + Ord,
{
    pub(crate) fn new(file_paths: &[TempPath]) -> Self {
        let mut queue = BTreeSet::new();

        let mut n = 0;
        for f in file_paths {
            let it = PeekFileRecordIterator::new(
                FileRecordIterator::new(BufReader::new(File::open(f).unwrap())),
                n,
            );

            if it.peeked().is_some() {
                queue.insert(it);
            }

            n += 1;
        }

        Self { queue }
    }
}

impl<R, T> Iterator for MergingIterator<R, T>
where
    PeekFileRecordIterator<R, T>: Ord,
    R: std::io::Read,
    T: for<'de> serde::Deserialize<'de>,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        let mut file_iterator = match self.queue.pop_first() {
            Some(v) => v,
            None => return None,
        };

        let ret = file_iterator.next().unwrap(); // only file_iters having a next item are contained in `self.queue`

        if let Some(_) = file_iterator.peek() {
            self.queue.insert(file_iterator);
        }

        Some(ret)
    }
}

pub(crate) enum CowForSC<'a, B: ?Sized + 'a>
where
    B: ToOwned,
{
    /// Borrowed data.
    Borrowed(&'a B),

    /// Owned Data
    Owned(<B as ToOwned>::Owned),
}

impl<B: ?Sized + ToOwned> Clone for CowForSC<'_, B> {
    fn clone(&self) -> Self {
        match *self {
            CowForSC::Borrowed(b) => CowForSC::Borrowed(b),
            CowForSC::Owned(ref o) => {
                let b: &B = o.borrow();
                CowForSC::Owned(b.to_owned())
            }
        }
    }

    fn clone_from(&mut self, source: &Self) {
        match (self, source) {
            (&mut CowForSC::Owned(ref mut dest), &CowForSC::Owned(ref o)) => {
                o.borrow().clone_into(dest)
            }
            (t, s) => *t = s.clone(),
        }
    }
}

impl<B: ?Sized + ToOwned> CowForSC<'_, B> {
    pub fn into_owned(self) -> <B as ToOwned>::Owned {
        match self {
            CowForSC::Borrowed(v) => v.to_owned(),
            CowForSC::Owned(v) => v,
        }
    }
}

impl<'a, B> fmt::Debug for CowForSC<'a, B>
where
    B: fmt::Debug + ToOwned,
    <B as ToOwned>::Owned: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CowForSC::Borrowed(v) => v.fmt(f),
            CowForSC::Owned(ref v) => v.fmt(f),
        }
    }
}

impl<B> Deref for CowForSC<'_, B>
where
    B: ?Sized + ToOwned,
    B::Owned: Borrow<B>,
{
    type Target = B;

    fn deref(&self) -> &Self::Target {
        match *self {
            CowForSC::Borrowed(v) => v,
            CowForSC::Owned(ref v) => v.borrow(),
        }
    }
}

impl<'a, B: ?Sized> Borrow<B> for CowForSC<'a, B>
where
    B: ToOwned,
{
    fn borrow(&self) -> &B {
        &**self
    }
}

/// MutBorrow or Owned enum.
/// 
/// This enum was made for being return values of `SortingCollection::iter_mut()`
pub(crate) enum MooForSC<'a, B: ?Sized + 'a>
where
    B: ToOwned,
{
    /// MutBorrowed data.
    MutBorrowed(&'a mut B),

    /// Owned Data
    Owned(<B as ToOwned>::Owned),
}


#[cfg(test)]
mod test {
    use std::{fs::File, str::FromStr};

    use crate::markdup::utils::read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates;

    use super::*;

    struct A<T: ?Sized> {
        m: T,
    }

    struct B {
        s: String,
    }

    impl B {
        fn consume(self) -> String {
            "ASDD".into()
        }

        fn ref_mut(&mut self) -> &mut String {
            &mut self.s
        }
    }

    #[test]
    fn mut_borrow() {
        let mut b = B { s: "ASD".into() };
        let mut v = vec![];
        let a = A { m: b.ref_mut() };

        v.push(b.consume());
    }

    #[test]
    fn save_and_load() {
        let mut sc = SortingCollection::<ReadEndsForMarkDuplicates>::new(
            1,
            false,
            vec![PathBuf::from_str(".tmp").unwrap()],
        );

        let mut a = ReadEndsForMarkDuplicates {
            score: 12,
            read1_index_in_file: 1,
            read2_index_in_file: 2,
            ..Default::default()
        };

        a.read_ends.read2_coordinate = 55534;

        let mut b = ReadEndsForMarkDuplicates {
            score: 25,
            read1_index_in_file: 33,
            read2_index_in_file: 22,
            ..Default::default()
        };

        b.read_ends.read2_coordinate = 12231;

        let mut c: ReadEndsForMarkDuplicates = ReadEndsForMarkDuplicates {
            score: 25,
            read1_index_in_file: 112,
            read2_index_in_file: 22,
            ..Default::default()
        };

        c.read_ends.read2_coordinate = 9978;

        sc.add(a).unwrap();

        sc.add(b).unwrap();

        sc.add(c).unwrap();

        println!("{:#?}", sc.iter().unwrap().collect::<Vec<_>>());
    }
}
