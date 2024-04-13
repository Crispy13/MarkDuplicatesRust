use std::collections::{BTreeMap, HashMap};

use rust_htslib::bam::{
    record::{RecordExt, SAMReadGroupRecord},
    HeaderView, Record,
};

// type Error = Box<dyn std::error::Error + Send + Sync>;
use anyhow::Error;

use crate::hts::utils::histogram::{f64H, Histogram};
pub(crate) struct LibraryIdGenerator {
    library_ids: HashMap<String, i16>, // from library string to library id
    next_library_id: i16,
    optical_duplicates_by_library_id: Histogram<i16>,
    pub(crate) duplicate_count_hist: Histogram<f64H>,
    pub(crate) non_optical_duplicate_count_hist: Histogram<f64H>,
    pub(crate) optical_duplicate_count_hist: Histogram<f64H>,

    metrics_by_library: BTreeMap<String, DuplicationMetrics>,
}

impl LibraryIdGenerator {
    const UNKNOWN_LIBRARY: &'static str = "Unknown Library";

    pub(crate) fn new() -> Self {
        Self {
            library_ids: HashMap::new(),
            next_library_id: 1,
            optical_duplicates_by_library_id: Histogram::new(),
            duplicate_count_hist: Histogram::from_labels("set_size", "all_sets"),
            non_optical_duplicate_count_hist: Histogram::from_labels(
                "set_size",
                "non_optical_sets",
            ),
            optical_duplicate_count_hist: Histogram::from_labels("set_size", "optical_sets"),
            metrics_by_library: BTreeMap::new(),
        }
    }

    pub(crate) fn get_read_group_libary_name<'m>(read_group: &SAMReadGroupRecord<'m>) -> &'m str {
        read_group.get_library().unwrap_or(Self::UNKNOWN_LIBRARY)
    }

    /** Get the library ID for the given SAM record. */
    pub(crate) fn get_library_id(&mut self, rec: &Record) -> Result<i16, Error> {
        let library = Self::get_library_name(rec)?;

        let library_id = match self.library_ids.get(library) {
            Some(v) => v,
            None => {
                self.library_ids
                    .insert(library.to_string(), self.next_library_id);
                self.next_library_id += 1;
                self.library_ids.get(library).unwrap()
            }
        }
        .clone();

        Ok(library_id)
    }

    /**
     * Gets the library name from the header for the record. If the RG tag is not present on
     * the record, or the library isn't denoted on the read group, a constant string is
     * returned.
     */
    pub(crate) fn get_library_name<'r>(rec: &'r Record) -> Result<&'r str, Error> {
        match rec.get_read_group()?.get_library() {
            Some(lb) => Ok(lb),
            None => Ok(Self::UNKNOWN_LIBRARY),
        }
    }

    pub(crate) fn get_optical_duplicates_by_library_id_map(&mut self) -> &mut Histogram<i16> {
        &mut self.optical_duplicates_by_library_id
    }

    pub(crate) fn get_non_optical_duplicate_count_hist(&mut self) -> &mut Histogram<f64H> {
        &mut self.non_optical_duplicate_count_hist
    }

    pub(crate) fn get_optical_duplicate_count_hist(&mut self) -> &mut Histogram<f64H> {
        &mut self.optical_duplicate_count_hist
    }

    pub(crate) fn get_duplicate_count_hist(&mut self) -> &mut Histogram<f64H> {
        &mut self.duplicate_count_hist
    }

    pub(crate) fn get_number_of_optical_duplicate_clusters(&self) -> i64 {
        self.optical_duplicates_by_library_id.get_sum_of_values() as i64
    }

    pub(crate) fn get_metrics_by_library(&self, library: &str) -> Option<&DuplicationMetrics> {
        self.metrics_by_library.get(library)
    }

    pub(crate) fn get_mut_metrics_by_library(&mut self, library: &str) -> Option<&mut DuplicationMetrics> {
        self.metrics_by_library.get_mut(library)
    }
    
    pub(crate) fn add_metrics_by_library(&mut self, library: String, metrics: DuplicationMetrics) {
        self.metrics_by_library.insert(library, metrics);
    }
}

#[allow(non_snake_case)]
#[derive(Default)]
pub(crate) struct DuplicationMetrics {
    UNMAPPED_READS: usize,
    SECONDARY_OR_SUPPLEMENTARY_RDS: usize,
    UNPAIRED_READS_EXAMINED: usize,
    READ_PAIRS_EXAMINED: usize,

    flow_based_metrics: Option<FlowBasedMetrics>,

    pub(crate) library:String,

    UNPAIRED_READ_DUPLICATES: usize,
    READ_PAIR_DUPLICATES: usize,
}

#[allow(non_snake_case)]
#[derive(Default)]
struct FlowBasedMetrics {
    UNPAIRED_WITH_TLEN: usize,
}

impl DuplicationMetrics {
    pub(crate) fn create_metrics(flow_metrics: bool) -> Self {
        // create based on the presence of flow order
        if !flow_metrics {
            Self::default()
        } else {
            Self {
                flow_based_metrics: Some(FlowBasedMetrics::default()),
                ..Default::default()
            }
        }
    }

    /**
     * Adds a read to the metrics
     */
    pub(crate) fn add_read_to_library_metrics(&mut self, rec: &Record) {
        // First bring the simple metrics up to date
        if rec.is_unmapped() {
            self.UNMAPPED_READS += 1;
        } else if rec.is_secondary() || rec.is_supplementary() {
            self.SECONDARY_OR_SUPPLEMENTARY_RDS += 1;
        } else if !rec.is_paired() || rec.is_mate_unmapped() {
            self.UNPAIRED_READS_EXAMINED += 1;
        } else {
            self.READ_PAIRS_EXAMINED += 1; // will need to be divided by 2 at the end
        }
    }
    
    /**
     * Adds duplicated read to the metrics
     */
    pub(crate) fn add_duplicate_read_to_metrics(&mut self, rec: &Record) {
        // only update duplicate counts for "decider" reads, not tag-a-long reads
        if !rec.is_secondary_or_supplementary() && !rec.is_unmapped() {
            // Update the duplication metrics
            if !rec.is_paired() || rec.is_mate_unmapped() {
                self.UNPAIRED_READ_DUPLICATES += 1;
            } else {
                self.READ_PAIRS_EXAMINED += 1; // will need to be divided by 2 at the end
            }
        }
    }
}

#[cfg(example)]
mod examples {
    use super::*;

    struct A {
        string: String,
    }

    struct B<'s> {
        str_slice: &'s str,
    }

    impl<'s> B<'s> {
        fn from_a(a: &'s A) -> Self {
            Self {
                str_slice: &a.string,
            }
        }

        fn get_str(&self) -> &'s str {
            self.str_slice
        }
    }

    fn get_str_slice<'a>(a: &'a A, s: &str) -> &'a str {
        let b = B::from_a(a);
        // let b = B {
        //     str_slice: &a.string,
        // };

        b.get_str()
    }
}
