use anyhow::anyhow;
use anyhow::Error;
use lru::LruCache;
use serde::Deserialize;
use serde::Serialize;
use std::fmt::Display;
use std::num::NonZeroUsize;
use std::{
    collections::HashMap,
    fs::{File, OpenOptions},
    io::BufWriter,
    path::{Path, PathBuf},
};

use std::hash::Hash;

use crate::hts::coordinate_sorted_pair_info_map::CoordinateSortedPairInfoMap;
use crate::hts::coordinate_sorted_pair_info_map::CoordinateSortedPairInfoMapExt;

pub(crate) struct FileAppendStreamLRUCache {
    cache_size: usize,
    inner: LruCache<PathBuf, BufWriter<File>>,
}

impl FileAppendStreamLRUCache {
    pub(crate) fn new(cache_size: usize) -> Self {
        assert!(cache_size != 0, "cache size must be non-zero");

        Self {
            cache_size,
            inner: LruCache::new(NonZeroUsize::new(cache_size).unwrap()),
        }
    }

    pub(crate) fn contains_key(&self, k: &PathBuf) -> bool {
        self.inner.contains(k)
    }

    // pub(crate) fn get_or_insert(&mut self, k: PathBuf) -> Result<&mut BufWriter<File>, Error> {
    //     match self.inner.entry(k) {
    //         std::collections::hash_map::Entry::Occupied(e) => Ok(e.into_mut()),
    //         std::collections::hash_map::Entry::Vacant(e) => {
    //             let path = e.key();

    //             let f = OpenOptions::new()
    //                 .write(true)
    //                 .append(true)
    //                 .create(true)
    //                 .open(path)
    //                 .or_else(|err| Err(anyhow!("err={:?}, path={:?}", err, path)))?;

    //             Ok(e.insert(BufWriter::new(f)))
    //         }
    //     }
    // }

    pub(crate) fn get_or_insert(&mut self, k: PathBuf) -> Result<&mut BufWriter<File>, Error> {
        if !self.inner.contains(&k) {
            let f = OpenOptions::new()
                .write(true)
                .append(true)
                .create(true)
                .open(&k)
                .or_else(|err| Err(anyhow!("err={:?}, path={:?}", err, &k)))?;

            self.inner.push(k.clone(), BufWriter::new(f));
        }

        Ok(self.inner.get_mut(&k).unwrap())
    }

    pub(crate) fn remove(&mut self, k: &PathBuf) -> Option<BufWriter<File>> {
        self.inner.pop(k)
    }

    // pub(crate) fn get_output_stream_for_sequence<K, R>(
    //     &mut self,
    //     pim: &mut CoordinateSortedPairInfoMap<K, R>,
    //     mate_sequence_index: i32,
    // ) -> Result<&mut BufWriter<File>, Error>
    // where
    //     K: std::cmp::Eq + PartialEq + Hash + Serialize + for<'de> Deserialize<'de> + Display,
    //     R: Serialize + for<'de> Deserialize<'de>,
    // {
    //     let path = pim.make_file_for_sequence(mate_sequence_index)?;

    //     self.get_or_insert(path)
    // }
}
