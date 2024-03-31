use std::{collections::{BTreeSet, HashSet}, fs::File, io::{BufReader, BufWriter, Read}, iter::Peekable, marker::PhantomData, mem::size_of, ops::{Deref, DerefMut}, path::PathBuf};

use serde::{Deserialize, Serialize};
use tempfile::{tempfile, NamedTempFile, TempPath};
use anyhow::Error;

use crate::{hts::utils::save_as_byte_to_file, utils::{get_allocated_and_resident_mem_of_app, human_readable_byte_count, mem_stats}};

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

    files:Vec<TempPath>,

    cleaned_up:bool,
}

impl<T> SortingCollection<T>
where
    T: Ord + Serialize + for<'de> Deserialize<'de>
{
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

                log::debug!(target: Self::log, 
                    "{} records in ram required approximately {} memory or {} per record. ", 
                    self.max_records_in_ram, 
                    human_readable_byte_count(used_bytes), 
                    human_readable_byte_count(used_bytes / self.max_records_in_ram));
            }
        }

        self.ram_records.push(rec);
        self.num_records_in_ram += 1;

        Ok(())
    }

    /**
     * Sort the records in memory, write them to a file, and clear the buffer of records in memory.
     */
    fn spill_to_disk(&mut self) -> Result<(), Error>{
        self.ram_records.get_mut(0..self.num_records_in_ram).unwrap_or_else(|| {
            panic!("Code error. Failed to get slice of `ram_records` with index 0..self.num_records_in_ram.");
        }).sort();

        let mut f = NamedTempFile::new()?;

        let mut buf_writer = BufWriter::new(f.as_file_mut());

        self.ram_records.drain(..self.num_records_in_ram).map(|r| {
            save_as_byte_to_file(&r, &mut buf_writer)
        }).collect::<Result<(), Error>>()?;

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
    fn iter(&mut self) -> Result<SortingCollectionIter<T>, Error> {
        if self.cleaned_up {
            panic!("Cannot call iterator() after cleanup() was called.");
        }

        self.done_adding()?;

        self.iteration_started = true;

        if self.files.is_empty() {
            Ok(self.in_memory_iter())
        } else {
            Ok(MergingIterator::new(&self.files))
        }
    }
    
    

    fn in_memory_iter(&mut self) -> std::slice::Iter<'_, T> {
        self.ram_records.get_mut(0..self.num_records_in_ram).unwrap_or_else(|| {
            panic!("Code error. Failed to get slice of `ram_records` with index 0..self.num_records_in_ram.");
        }).sort();

        self.ram_records.get(..self.num_records_in_ram).unwrap().iter()

    }
}


pub(crate) enum SortingCollectionIter<'a, T>
where
    T: for<'de> Deserialize<'de>,
{
    InMemoryIter(std::slice::Iter<'a, T>),
    MergingIterator(MergingIterator<BufReader<File>, T>)
}

impl<'a, T> Iterator for SortingCollectionIter<'a, T>
where
    T: for<'de> Deserialize<'de>,
    PeekFileRecordIterator<BufReader<File>, T>: Ord
{
    type Item=&'a T;

    fn next(&mut self) -> Option<Self::Item> {
        todo!();

        match self {
            SortingCollectionIter::InMemoryIter(it) => it.next(),
            SortingCollectionIter::MergingIterator(it) => it.next().as_ref(),
        }

    }
}

struct FileRecordIterator<R, T> {
    inner: R,
    buf: Vec<u8>,

    phantom_data: PhantomData<T>,
}

impl<R, T> FileRecordIterator<R, T> {
    fn new(inner: R, ) -> Self {
        Self { inner, buf:Vec::with_capacity(size_of::<T>()), phantom_data: PhantomData }
    }
}

impl<T:for<'de> Deserialize<'de>, R:Read> Iterator for FileRecordIterator<R, T> {
    type Item=T;

    fn next(&mut self) -> Option<Self::Item> {
        let buf = &mut self.buf;
        
        buf.clear();

        match self.inner.read_exact(buf) {
            Ok(_) => {},
            Err(err) => {
                match err.kind() {
                    std::io::ErrorKind::UnexpectedEof => return None,
                    err => panic!("{:?}", err),
                }
            },
        }

        match bincode::deserialize(&buf) {
            Ok(v) => Some(v),
            Err(err) => panic!("{:?}", err),
        }
    }
}


struct PeekFileRecordIterator<R, T> where
T: for<'de> Deserialize<'de>,
R: Read {
    inner: Peekable<FileRecordIterator<R, T>>,
    n: i32,
}

impl<R, T> PeekFileRecordIterator<R, T>
where
T: for<'de> Deserialize<'de>,
R: Read
{
    fn new(inner: Peekable<FileRecordIterator<R, T>>, n: i32) -> Self {
        Self { inner, n }
    }
}

impl<R,T> Deref for PeekFileRecordIterator<R, T>
where
    T: for<'de> Deserialize<'de>,
    R: Read,
{
    type Target=Peekable<FileRecordIterator<R, T>>;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<R,T> DerefMut for PeekFileRecordIterator<R, T>
where
    T: for<'de> Deserialize<'de>,
    R: Read,
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

pub(crate) struct MergingIterator<R, T>
where
    T: for<'de> Deserialize<'de>,
    R: Read,
{
    queue: BTreeSet<PeekFileRecordIterator<R, T>>,

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
    R: Read, {
        
    }

impl<R, T> PartialOrd for PeekFileRecordIterator<R, T>
where
    T: for<'de> Deserialize<'de> + PartialEq + PartialOrd,
    R: Read,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.peek().partial_cmp(&other.peek()) {
            Some(ord) => return Some(ord),
            None => {},
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
        self.peek() == other.peek()
    }
}

impl<T> MergingIterator<BufReader<File>, T>
where
    T: for<'de> Deserialize<'de> + Ord,
{
    pub(crate) fn new(file_paths:&[TempPath]) -> Self {
        let mut queue = BTreeSet::new();

        let mut n =0;
        for f in file_paths {
            let mut it = FileRecordIterator::new(BufReader::new(File::open(f).unwrap())).peekable();

            if it.peek().is_some() {
                queue.insert(
                    
                        PeekFileRecordIterator::new(
                            it,
                            n,
                        )
                    
                );

                n+=1;
            }
        }

        Self {
            queue,
        }
    }
}


impl<R, T> Iterator for MergingIterator<R, T>
where
    PeekFileRecordIterator<R, T>: Ord,
    R: std::io::Read,
    T: for<'de> serde::Deserialize<'de>
{
    type Item=T;

    fn next(&mut self) -> Option<Self::Item> {
        let mut file_iterator = match self.queue.pop_first() {
            Some(v) => v,
            None => return None,
        };

        let ret = file_iterator.next().unwrap();

        if let Some(_) = file_iterator.peek() {
            self.queue.insert(file_iterator);
        } 
        
        Some(ret)
    }
}





#[cfg(test)]
mod test {
    use std::fs::File;

    use super::*;

    struct A<T:?Sized> {
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
        let mut b = B{s:"ASD".into()};
        let mut v = vec![];
        let a = A {m: b.ref_mut()};

        v.push(b.consume());
        
    }
}
