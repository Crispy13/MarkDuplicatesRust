pub(crate) mod coordinate_sorted_pair_info_map;

use rust_htslib::bam::HeaderView;

type Error = Box<dyn std::error::Error + Send + Sync>;

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
            _ => Err(format!("Invalid value: {}", s))?,
        };

        Ok(v)
    }
    
    pub(crate) fn from_header(h: &HeaderView) -> Result<Self, Error> {
        let v =match h.header_map().get_sort_order() {
            Some(s) => SortOrder::from_str(s)?,
            None => Err("Failed to get 'SO' value in sam file header.")?,
        };

        Ok(v)
    }
}

#[cfg(test)]
mod test {
    use rust_htslib::bam::{IndexedReader, Read};

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
}
