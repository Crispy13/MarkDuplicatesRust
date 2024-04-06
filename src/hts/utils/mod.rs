use anyhow::Error;
use serde::{Deserialize, Serialize};
use std::{
    fs::File,
    io::{Read, Write},
    mem::size_of,
};

pub(crate) mod murmur3;
pub(crate) mod sorting_collection;
pub(crate) mod file_append_stream_lru_cache;

pub(crate) fn save_as_byte_to_file<T: Serialize>(
    value: &T,
    file: &mut impl Write,
) -> Result<(), Error> {
    // let ser_bytes = bincode::serialize(value)?;

    // assert_eq!(ser_bytes.len(), size_of::<T>());

    // file.write_all(&ser_bytes)?;
    bincode::serialize_into(file, value)?;

    Ok(())
}

pub(crate) fn load_byte_file_as_obj_v1<T: for<'de> Deserialize<'de>>(
    file: &mut impl Read,
) -> Result<T, Error> {
    let size_t: usize = size_of::<T>();

    let mut byte_buf = Vec::with_capacity(size_t);

    file.read_exact(&mut byte_buf)?;

    if byte_buf.len() != size_t {
        Err(crate::utils::errors::Error::FailedDeserializeFromByteFile {
            eb:size_t,
            rb:byte_buf.len(),
            input_file: String::new(),
        })?
    }

    Ok(bincode::deserialize(&byte_buf)?)
}

pub(crate) fn load_byte_file_as_obj<T: for<'de> Deserialize<'de>>(
    file: &mut impl Read,
) -> Result<T, Error> {
    Ok(bincode::deserialize_from::<_, T>(file)?)
}
