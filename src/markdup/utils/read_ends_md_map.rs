use std::collections::HashMap;

use crate::hts::coordinate_sorted_pair_info_map::{
    CoordinateSortedPairInfoMap, CoordinateSortedPairInfoMapExt,
};

use super::read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates;

// common error type
// type Error = Box<dyn std::error::Error + Send + Sync>;
use anyhow::{anyhow, Error};

pub(crate) struct DiskBasedReadEndsForMarkDuplicatesMap {
    pair_info_map: CoordinateSortedPairInfoMap<String, ReadEndsForMarkDuplicates>,
}

impl DiskBasedReadEndsForMarkDuplicatesMap {
    pub(crate) fn new(max_open_files: i32) -> Result<Self, Error> {
        Ok(Self {
            pair_info_map: CoordinateSortedPairInfoMap::new(max_open_files)?,
        })
    }

    fn remove(
        &mut self,
        mate_sequence_index: i32,
        key: &str,
    ) -> Result<Option<ReadEndsForMarkDuplicates>, Error> {
        self.pair_info_map.remove(mate_sequence_index, key)
    }

    fn put(
        &mut self,
        mate_sequence_index: i32,
        key: String,
        read_ends: ReadEndsForMarkDuplicates,
    ) -> Result<(), Error> {
        self.pair_info_map.put(mate_sequence_index, key, read_ends)
    }
}

pub(crate) struct MemoryBasedReadEndsForMarkDuplicatesMap {
    map_per_sequence: Vec<HashMap<String, ReadEndsForMarkDuplicates>>,
}

impl MemoryBasedReadEndsForMarkDuplicatesMap {
    pub(crate) fn new() -> MemoryBasedReadEndsForMarkDuplicatesMap {
        Self {
            map_per_sequence: Vec::new(),
        }
    }

    fn remove(
        &mut self,
        mate_sequence_index: usize,
        key: &str,
    ) -> Option<ReadEndsForMarkDuplicates> {
        self.map_per_sequence
            .get_mut(mate_sequence_index)
            .and_then(|m| m.remove(key))
    }

    fn put(
        &mut self,
        mate_sequence_index: usize,
        key: String,
        read_ends: ReadEndsForMarkDuplicates,
    ) {
        while mate_sequence_index > self.map_per_sequence.len() {
            self.map_per_sequence.push(HashMap::new())
        }

        self.map_per_sequence
            .push(HashMap::from([(key, read_ends)]));
    }

    fn size(&self) -> usize {
        self.map_per_sequence.iter().fold(0, |a, b| a + b.len())
    }

    /**
     * @return number of elements stored in RAM.  Always <= size()
     */
    fn size_in_ram(&self) -> usize {
        self.size()
    }
}

pub(crate) enum ReadEndsForMarkDuplicatesMap {
    MemoryBased(MemoryBasedReadEndsForMarkDuplicatesMap),
    DiskBased(DiskBasedReadEndsForMarkDuplicatesMap),
}

impl ReadEndsForMarkDuplicatesMap {
    pub(crate) fn remove(
        &mut self,
        mate_sequence_index: i32,
        key: &str,
    ) -> Result<Option<ReadEndsForMarkDuplicates>, Error> {
        match self {
            ReadEndsForMarkDuplicatesMap::MemoryBased(m) => {
                Ok(m.remove(mate_sequence_index as usize, key))
            }

            ReadEndsForMarkDuplicatesMap::DiskBased(m) => Ok(m.remove(mate_sequence_index, key)?),
        }
    }

    pub(crate) fn put(
        &mut self,
        mate_sequence_index: i32,
        key: String,
        read_ends: ReadEndsForMarkDuplicates,
    ) -> Result<(), Error> {
        match self {
            ReadEndsForMarkDuplicatesMap::MemoryBased(m) => {
                Ok(m.put(mate_sequence_index as usize, key, read_ends))
            }
            ReadEndsForMarkDuplicatesMap::DiskBased(m) => {
                m.put(mate_sequence_index, key, read_ends)
            }
        }
    }
}

pub(crate) trait ReadEndsForMarkDuplicatesMapExt {}
