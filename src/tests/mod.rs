use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read},
    path::Path,
};

use serde::ser::{SerializeSeq, Serializer};
use serde::{Deserialize, Serialize};

use crate::markdup::utils::{
    read_ends::ReadEnds, read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates,
};

// pub(crate) fn save_object_to_json<T: Serialize>(value: &T, path: impl AsRef<Path>) {
//     let pf = serde_json::ser::PrettyFormatter::with_indent(b"    ");

//     let buf_writer = BufWriter::new(File::create(path.as_ref()).unwrap());

//     let mut ser = serde_json::Serializer::with_formatter(buf_writer, pf);

//     serde_json::to_writer(ser, value);
// }

pub(crate) trait ToJsonFileForTest {
    fn save_object_to_json(self, path: impl AsRef<Path>);
}

impl<I> ToJsonFileForTest for I
where
    I: Iterator,
    <I as Iterator>::Item: Serialize,
{
    fn save_object_to_json(self, path: impl AsRef<Path>) {
        let pf = serde_json::ser::PrettyFormatter::with_indent(b"    ");

        let buf_writer = BufWriter::new(File::create(path.as_ref()).unwrap());

        let mut ser = serde_json::Serializer::with_formatter(buf_writer, pf);

        let mut ser_seq = ser.serialize_seq(None).unwrap();

        self.for_each(|e| ser_seq.serialize_element(&e).unwrap());

        ser_seq.end().unwrap();

        // self.iter().for_each(|re| {re.serialize(&mut ser).unwrap();});
    }
}




#[allow(non_snake_case)]
#[derive(Deserialize, Debug)]
pub(crate) struct JavaMDReadEnds {
    qname: String,
    score: i16,
    read1IndexInFile: i64,
    read2IndexInFile: i64,
    duplicateSetSize: i32,
    libraryId: i16,
    orientation: i8,
    read1ReferenceIndex: i32,
    read1Coordinate: i32,
    read2ReferenceIndex: i32,
    read2Coordinate: i32,
    readGroup: i16,
    orientationForOpticalDuplicates: i8,
    tile: i16,
    x: i32,
    y: i32,
}

pub(crate) fn from_java_read_ends(jmd: JavaMDReadEnds) -> ReadEndsForMarkDuplicates {
    let mut read_ends = ReadEnds::default();

    read_ends.library_id = jmd.libraryId;
    read_ends.orientation = jmd.orientation;
    read_ends.read1_reference_index = jmd.read1ReferenceIndex;
    read_ends.read2_reference_index = jmd.read2ReferenceIndex;

    read_ends.read1_coordinate = jmd.read1Coordinate;
    read_ends.read2_coordinate = jmd.read2Coordinate;
    read_ends.read_group = jmd.readGroup;

    read_ends.orientation_for_optical_duplicates = jmd.orientationForOpticalDuplicates;

    read_ends.pls.tile = jmd.tile;
    read_ends.pls.x = jmd.x;
    read_ends.pls.y = jmd.y;

    let read_ends_for_md = ReadEndsForMarkDuplicates {
        read1_qname:jmd.qname,
        score: jmd.score,
        read1_index_in_file: jmd.read1IndexInFile,
        read2_index_in_file: jmd.read2IndexInFile,
        duplicate_set_size: jmd.duplicateSetSize,
        read_ends,
        barcode_data: None,
    };

    read_ends_for_md
}

fn from_java_read_ends_to_rust_read_ends(path: impl AsRef<Path>) -> Vec<ReadEndsForMarkDuplicates> {
    let buf_reader = BufReader::new(File::open(path.as_ref()).unwrap());

    let jmds = serde_json::from_reader::<_, Vec<JavaMDReadEnds>>(buf_reader).unwrap();

    jmds.into_iter().map(from_java_read_ends).collect()
}

pub(crate) struct JavaReadEndIterator {
    buf_reader: BufReader<File>,
    skipped: [u8; 1],
    read_buf: Vec<u8>,
    end: bool,
}

impl JavaReadEndIterator {
    pub(crate) fn new(file_path: impl AsRef<Path>) -> Self {
        let mut buf_reader = BufReader::new(File::open(file_path.as_ref()).unwrap());

        let mut skipped = [0_u8; 1];

        if Self::read_one_byte(&mut buf_reader, &mut skipped) != b'[' {
            panic!()
        }

        // let serde_iter = serde_json::from_reader::<_, JavaMDReadEnds>(file).into_iter();

        Self {
            buf_reader,
            skipped,
            end: false,
            read_buf: Vec::with_capacity(128),
        }
    }

    fn read_one_byte(file: &mut BufReader<File>, skipped: &mut [u8]) -> u8 {
        file.read_exact(skipped).unwrap();

        skipped[0]
    }
}

impl Iterator for JavaReadEndIterator {
    type Item = JavaMDReadEnds;

    fn next(&mut self) -> Option<Self::Item> {
        if self.end {
            return None;
        }

        self.buf_reader
            .read_until(b'}', &mut self.read_buf)
            .unwrap();

        let item = serde_json::from_slice::<JavaMDReadEnds>(self.read_buf.as_slice()).unwrap();

        self.read_buf.clear();

        if Self::read_one_byte(&mut self.buf_reader, &mut self.skipped) != b',' {
            self.end = true;

            self.buf_reader.read_to_end(&mut self.read_buf).unwrap();

            if self.read_buf.contains(&b'{') {
                panic!(
                    "Parsing failed. remained_bytes={}",
                    std::str::from_utf8(&self.read_buf).unwrap()
                )
            }
        }

        Some(item)
    }
}
