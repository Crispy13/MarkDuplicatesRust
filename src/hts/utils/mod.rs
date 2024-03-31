use std::{fs::File, io::{Read, Write}, mem::size_of};
use anyhow::Error;
use serde::{Deserialize, Serialize};

pub(crate) mod murmur3;
pub(crate) mod sorting_collection;

pub(crate) fn save_as_byte_to_file<T: Serialize>(value: &T, file: &mut impl Write) -> Result<(), Error> {
    bincode::serialize_into(file, value)?;
    
    Ok(())
}


pub(crate) fn load_byte_file_as_obj<T: for<'de> Deserialize<'de>>(file: &mut impl Read) -> Result<T, Error> {
    let mut byte_buf = Vec::with_capacity(size_of::<T>());

    file.read_exact(&mut byte_buf)?;

    Ok(bincode::deserialize(&byte_buf)?)
}