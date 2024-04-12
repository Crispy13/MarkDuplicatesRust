pub(crate) mod coordinate_sorted_pair_info_map;
pub(crate) mod duplicate_scoring_strategy;
pub(crate) mod utils;
pub(crate) mod header;

use clap::ValueEnum;
use rust_htslib::bam::HeaderView;

// type Error = Box<dyn std::error::Error + Send + Sync>;
use anyhow::{anyhow, Error};

pub(crate) enum SAMTag {
    PG,
}

impl SAMTag {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            SAMTag::PG => "PG",
        }
    }
}

#[derive(Debug, Clone, ValueEnum)]
pub(crate) enum SortOrder {
    Unknown,
    Unsorted,
    QueryName,
    Coordinate,
}

impl SortOrder {
    pub(crate) const SO: &'static str = "SO";
    pub(crate) const GO: &'static str = "GO";

    const UNKNOWN: &'static str = "unknown";
    const UNSORTED: &'static str = "unsorted";
    const QUERYNAME: &'static str = "queryname";
    const COORDINATE: &'static str = "coordinate";


    pub(crate) fn from_str(s: &str) -> Result<Self, Error> {
        let v = match s {
            Self::UNKNOWN => Self::Unknown,
            Self::UNSORTED => Self::Unsorted,
            Self::QUERYNAME => Self::QueryName,
            Self::COORDINATE => Self::Coordinate,
            _ => Err(anyhow!("Invalid value: {}", s))?,
        };

        Ok(v)
    }

    pub(crate) fn from_header(h: &HeaderView) -> Result<Self, Error> {
        let v = match h.header_map().get_sort_order() {
            Some(s) => SortOrder::from_str(s)?,
            None => Err(anyhow!("Failed to get 'SO' value in sam file header."))?,
        };

        Ok(v)
    }
}

impl std::fmt::Display for SortOrder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SortOrder::Unknown => write!(f, "{}", Self::UNKNOWN),
            SortOrder::Unsorted => write!(f, "{}", Self::UNSORTED),
            SortOrder::QueryName => write!(f, "{}", Self::QUERYNAME),
            SortOrder::Coordinate => write!(f, "{}", Self::COORDINATE),
        }
    }
}

#[cfg(test)]
mod test {
    use rust_htslib::bam::{IndexedReader, Read, Record};

    use super::SortOrder;

    #[test]
    fn get_sort_order() {
        let reader = IndexedReader::from_path(
            "/home/eck/workspace/markdup_rust/tests/data/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
        ).unwrap();

        // println!("{:#?}", reader.header().header_map().get("HD").unwrap().first().unwrap().get("SO").unwrap());
        println!(
            "{:?}",
            SortOrder::from_str(reader.header().header_map().get_sort_order().unwrap()).unwrap()
        );
    }

    #[test]
    fn set_value() {
        let mut reader = IndexedReader::from_path(
            "/home/eck/workspace/markdup_rust/tests/data/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
        ).unwrap();

        let mut record = Record::new();

        reader.fetch("11").unwrap();

        while let Some(r) = reader.read(&mut record) {
            r.unwrap();

            println!("{}", std::str::from_utf8(record.qname()).unwrap());
            println!("{}", record.mapq());

            record.set_mapq(59);

            println!("{}", record.mapq());

            break;
        }

        reader.fetch("11").unwrap();

        while let Some(r) = reader.read(&mut record) {
            r.unwrap();

            println!("{}", std::str::from_utf8(record.qname()).unwrap());
            println!("{}", record.mapq());

            break;
        }
    }
}
