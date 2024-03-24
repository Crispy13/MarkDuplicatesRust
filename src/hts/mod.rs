pub(crate) mod coordinate_sorted_pair_info_map;
pub(crate) mod duplicate_scoring_strategy;
pub(crate) mod utils;

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

#[derive(Debug)]
pub(crate) enum SortOrder {
    Unknown,
    Unsorted,
    QueryName,
    Coordinate,
}

impl SortOrder {
    pub(crate) fn from_str(s: &str) -> Result<Self, Error> {
        let v = match s {
            "unknown" => Self::Unknown,
            "unsorted" => Self::Unsorted,
            "queryname" => Self::QueryName,
            "coordinate" => Self::Coordinate,
            _ => Err(anyhow!("Invalid value: {}", s))?,
        };

        Ok(v)
    }
    
    pub(crate) fn from_header(h: &HeaderView) -> Result<Self, Error> {
        let v =match h.header_map().get_sort_order() {
            Some(s) => SortOrder::from_str(s)?,
            None => Err(anyhow!("Failed to get 'SO' value in sam file header."))?,
        };

        Ok(v)
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
