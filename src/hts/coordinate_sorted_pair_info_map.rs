use std::{
    borrow::Borrow,
    collections::{HashMap, HashSet},
    fmt::Display,
    fs::{create_dir, remove_file, File, OpenOptions},
    hash::Hash,
    io::{Read, Write},
    mem::size_of,
    path::{Path, PathBuf},
    str::FromStr,
};

use serde::{Deserialize, Serialize};

use crate::utils::PathExt;

// type Error = Box<dyn std::error::Error + Send + Sync>;
use anyhow::{anyhow, Error};

struct FileAppendStreamLRUCache {
    cache_size: i32,
}

pub(crate) struct CoordinateSortedPairInfoMap<K, R> {
    work_dir: PathBuf,
    sequence_index_of_map_in_ram: i32,
    map_in_ram: Option<HashMap<K, R>>,
    output_streams: FileAppendStreamLRUCache,
    // element_codec: Codec<K, R>,
    size_of_map_on_disk: HashMap<i32, usize>,
    iteration_in_progress: bool,

    files_to_delete_on_drop: HashSet<String>,
}

impl<K, R> CoordinateSortedPairInfoMap<K, R>
where
    K: std::cmp::Eq + PartialEq + Hash + Serialize + for<'de> Deserialize<'de> + Display,
    R: Serialize + for<'de> Deserialize<'de>,
{
    const INVALID_SEQUENCE_INDEX: i32 = -2;

    pub(crate) fn new(
        max_open_file: i32,
        // element_codec: Codec<K, R>,
    ) -> Result<Self, Error> {
        let work_dir = PathBuf::from_str("CSPI.tmp")?;

        create_dir(&work_dir)?;

        Ok(Self {
            work_dir,
            sequence_index_of_map_in_ram: -2,
            map_in_ram: None,
            size_of_map_on_disk: HashMap::new(),
            iteration_in_progress: false,
            output_streams: FileAppendStreamLRUCache {
                cache_size: max_open_file,
            },
            files_to_delete_on_drop: HashSet::new(),
        })
    }

    fn size(&self) -> usize {
        self.size_of_map_on_disk
            .values()
            .fold(self.size_in_ram(), |a, b| a + b)
    }

    fn size_in_ram(&self) -> usize {
        match self.map_in_ram {
            Some(ref m) => m.len(),
            None => 0,
        }
    }
}

pub(crate) trait CoordinateSortedPairInfoMapExt<K, R>
where
    K: std::cmp::Eq + PartialEq + Hash + Serialize + for<'de> Deserialize<'de> + Display,
    R: Serialize + for<'de> Deserialize<'de>,
{
    type DataCodec;

    fn remove<Q: ?Sized>(&mut self, sequence_index: i32, key: &Q) -> Result<R, Error>
    where
        K: Borrow<Q>,
        Q: Hash + Eq;

    fn get_output_stream_for_sequence(&mut self, mate_sequence_index: i32) -> Result<File, Error>;

    fn ensure_sequence_loaded(&mut self, sequence_index: i32) -> Result<(), Error>;

    fn make_file_for_sequence(&mut self, index: i32) -> Result<PathBuf, Error>;

    fn put(&mut self, sequence_index: i32, key: K, record: R) -> Result<(), Error>;
}

impl<K, R> CoordinateSortedPairInfoMapExt<K, R> for CoordinateSortedPairInfoMap<K, R>
where
    K: std::cmp::Eq + PartialEq + Hash + Serialize + for<'de> Deserialize<'de> + Display,
    R: Serialize + for<'de> Deserialize<'de>,
{
    type DataCodec = (K, R);

    fn remove<Q: ?Sized>(&mut self, sequence_index: i32, key: &Q) -> Result<R, Error>
    where
        K: Borrow<Q>,
        Q: Hash + Eq,
    {
        if self.iteration_in_progress {
            Err(anyhow!("Cannot be called when iteration is in progress"))?
        } else {
            self.ensure_sequence_loaded(sequence_index)?;
            Ok(self.map_in_ram.as_mut().unwrap().remove(key).unwrap())
        }
    }

    fn get_output_stream_for_sequence(&mut self, mate_sequence_index: i32) -> Result<File, Error> {
        Ok(OpenOptions::new()
            .write(true)
            .append(true)
            .open(self.make_file_for_sequence(mate_sequence_index)?)?)
    }

    fn ensure_sequence_loaded(&mut self, sequence_index: i32) -> Result<(), Error> {
        if self.sequence_index_of_map_in_ram != sequence_index {
            let mut map_on_disk: PathBuf;
            if self.map_in_ram.is_some() {
                map_on_disk = self.make_file_for_sequence(self.sequence_index_of_map_in_ram)?;

                // file existence is checked in `make_file_for_sequence`

                if !self.map_in_ram.as_ref().unwrap().is_empty() {
                    let mut os =
                        self.get_output_stream_for_sequence(self.sequence_index_of_map_in_ram)?;

                    let map_in_ram = self.map_in_ram.as_mut().unwrap();

                    self.size_of_map_on_disk
                        .insert(self.sequence_index_of_map_in_ram, map_in_ram.len());

                    // write read ends to file
                    map_in_ram
                        .drain()
                        .map(|(k, v)| save_to_byte_file(&(k, v), &mut os))
                        .collect::<Result<Vec<_>, Error>>()?;
                }
            } else {
                self.map_in_ram = Some(HashMap::new());
            }

            self.sequence_index_of_map_in_ram = sequence_index;
            map_on_disk = self.make_file_for_sequence(sequence_index)?;

            let num_records = self.size_of_map_on_disk.remove(&sequence_index);
            if map_on_disk.exists() {
                if num_records.is_none() {
                    Err(anyhow!("No num records for {}", map_on_disk.try_to_str()?))?
                }

                let mut input_stream;

                input_stream = File::open(&map_on_disk)?;

                for i in (0..num_records.unwrap()) {
                    let (key, record) =
                        load_byte_file_as_obj::<Self::DataCodec>(&mut input_stream)?;
                    let map_in_ram = self.map_in_ram.as_mut().unwrap();
                    if map_in_ram.contains_key(&key) {
                        Err(anyhow!(
                            "Value was put into PairInfoMap more than once.  {}: {}",
                            sequence_index, key
                        ))?
                    }

                    map_in_ram.insert(key, record);
                }

                drop(input_stream);

                remove_file(&map_on_disk)?;
            } else if let Some(true) = num_records.and_then(|nr| Some(nr > 0)) {
                Err(anyhow!(
                    "Non-zero numRecords but {} does not exist",
                    map_on_disk.try_to_str()?
                ))?
            }
        }

        Ok(())
    }

    fn make_file_for_sequence(&mut self, index: i32) -> Result<PathBuf, Error> {
        let fp = format!(
            "{}/{}.tmp",
            self.work_dir
                .to_str()
                .ok_or_else(|| anyhow!("Cannot convert `work_dir` into str."))?,
            index
        );

        self.files_to_delete_on_drop.insert(fp.clone());

        Ok(PathBuf::from_str(&fp)?)
    }

    fn put(&mut self, sequence_index: i32, key: K, record: R) -> Result<(), Error> {
        if self.iteration_in_progress {
            Err(anyhow!("Cannot be called when iteration is in progress"))?
        } else {
            if sequence_index == self.sequence_index_of_map_in_ram {
                // just unwrap(). we can assure it is not None when this method is called.
                let map_in_ram = self.map_in_ram.as_mut().unwrap();
                if map_in_ram.contains_key(&key) {
                    Err(anyhow!(
                        "Putting value into PairInfoMap that already existed. {}: {}",
                        sequence_index, key
                    ))?
                }

                map_in_ram.insert(key, record);
            } else {
                let mut os = self.get_output_stream_for_sequence(sequence_index)?;
                save_to_byte_file(&(key, record), &mut os)?;

                let prev_count = match self.size_of_map_on_disk.get(&sequence_index) {
                    Some(v) => *v,
                    None => 0,
                };

                self.size_of_map_on_disk
                    .insert(sequence_index, prev_count + 1);
            }
        }

        todo!()
    }
}

impl<K, R> Drop for CoordinateSortedPairInfoMap<K, R> {
    fn drop(&mut self) {
        self.files_to_delete_on_drop
            .iter()
            .for_each(|fp| remove_file(fp).unwrap_or_else(|err| ()));
        //TODO: log a message when removing fails as debug level
    }
}

fn save_to_byte_file<T: Serialize>(value: &T, file: &mut File) -> Result<(), Error> {
    let encoded = bincode::serialize(value)?;

    file.write_all(&encoded)?;

    Ok(())
}

fn load_byte_file_as_obj<T: for<'de> Deserialize<'de>>(file: &mut File) -> Result<T, Error> {
    let mut byte_buf = Vec::with_capacity(size_of::<T>());

    file.read_exact(&mut byte_buf)?;

    Ok(bincode::deserialize(&byte_buf)?)
}

// pub(crate) trait SaveToByteFile {
//     fn save_to_byte_file() {

//     }
// }

#[cfg(test)]
mod test {
    use std::mem::size_of;

    use crate::markdup::utils::read_ends_for_mark_duplicates_with_barcodes::ReadEndsBarcodeData;

    use super::*;
    #[test]
    fn save_and_load() {
        let a = ReadEndsBarcodeData {
            barcode: 1,
            read_one_barcode: 2,
            read_two_barcode: 44,
        };

        let a_ser = bincode::serialize(&a).unwrap();

        let fp = "test_a.bytes";
        save_to_byte_file(&a, &mut File::create(fp).unwrap()).unwrap();
        println!("{}", size_of::<ReadEndsBarcodeData>());
        println!("{}", a_ser.len());

        let a_deser: ReadEndsBarcodeData = {
            let mut byte_buf = Vec::new();

            File::open(fp).unwrap().read_to_end(&mut byte_buf).unwrap();

            println!("{}", byte_buf.len());

            bincode::deserialize(&byte_buf).unwrap()
        };

        println!("{:#?}", a_deser);
    }
}
