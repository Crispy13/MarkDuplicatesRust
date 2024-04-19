use std::collections::{BTreeMap, HashMap};

use rust_htslib::bam::{extra_ext::SAMReadGroupRecord, record, HeaderView, Record};

use rust_htslib::bam::extra_ext::RecordExt;

// type Error = Box<dyn std::error::Error + Send + Sync>;
use anyhow::{anyhow, Error};

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

    pub(crate) fn get_mut_optical_duplicates_by_library_id_map(&mut self) -> &mut Histogram<i16> {
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

    pub(crate) fn get_mut_metrics_by_library(
        &mut self,
        library: &str,
    ) -> Option<&mut DuplicationMetrics> {
        self.metrics_by_library.get_mut(library)
    }

    pub(crate) fn add_metrics_by_library(&mut self, library: String, metrics: DuplicationMetrics) {
        self.metrics_by_library.insert(library, metrics);
    }

    pub(crate) fn get_metrics_by_library_map(
        &self,
    ) -> &BTreeMap<std::string::String, DuplicationMetrics> {
        &self.metrics_by_library
    }

    pub(crate) fn get_mut_metrics_by_library_map(
        &mut self,
    ) -> &mut BTreeMap<std::string::String, DuplicationMetrics> {
        &mut self.metrics_by_library
    }

    pub(crate) fn get_optical_duplicates_by_library_id_map(&self) -> &Histogram<i16> {
        &self.optical_duplicates_by_library_id
    }

    pub(crate) fn get_library_ids_map(&self) -> &HashMap<String, i16> {
        &self.library_ids
    }
}

#[allow(non_snake_case)]
#[derive(Default, Clone)]
pub(crate) struct DuplicationMetrics {
    UNMAPPED_READS: usize,
    SECONDARY_OR_SUPPLEMENTARY_RDS: usize,
    UNPAIRED_READS_EXAMINED: usize,
    pub(crate) READ_PAIRS_EXAMINED: usize,

    flow_based_metrics: Option<FlowBasedMetrics>,

    pub(crate) library: String,

    UNPAIRED_READ_DUPLICATES: usize,
    pub(crate) READ_PAIR_DUPLICATES: usize,
    pub(crate) READ_PAIR_OPTICAL_DUPLICATES: usize,

    ESTIMATED_LIBRARY_SIZE: Option<i64>,

    PERCENT_DUPLICATION: f64,
}

#[allow(non_snake_case)]
#[derive(Default, Clone)]
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

    /**
     * Fills in the ESTIMATED_LIBRARY_SIZE based on the paired read data examined where
     * possible and the PERCENT_DUPLICATION.
     */
    pub(crate) fn calculate_derived_fields(&self) -> Result<(), Error> {
        self.ESTIMATED_LIBRARY_SIZE = Self::estimate_library_size(
            self.READ_PAIRS_EXAMINED - self.READ_PAIR_OPTICAL_DUPLICATES,
            self.READ_PAIRS_EXAMINED - self.READ_PAIR_DUPLICATES,
        )?;

        if self.UNPAIRED_READS_EXAMINED + self.READ_PAIRS_EXAMINED != 0 {
            self.PERCENT_DUPLICATION = (self.UNPAIRED_READ_DUPLICATES
                + self.READ_PAIR_DUPLICATES * 2) as f64
                / (self.UNPAIRED_READS_EXAMINED + self.READ_PAIRS_EXAMINED * 2) as f64;
        } else {
            self.PERCENT_DUPLICATION = 0.0;
        }

        Ok(())
    }

    /**
     * Estimates the size of a library based on the number of paired end molecules observed
     * and the number of unique pairs observed.
     * <p>
     * Based on the Lander-Waterman equation that states:
     * C/X = 1 - exp( -N/X )
     * where
     * X = number of distinct molecules in library
     * N = number of read pairs
     * C = number of distinct fragments observed in read pairs
     */
    fn estimate_library_size(
        read_pairs: usize,
        unique_read_pairs: usize,
    ) -> Result<Option<i64>, Error> {
        let read_pair_duplicates = read_pairs - unique_read_pairs;

        if read_pairs > 0 && read_pair_duplicates > 0 {
            let unique_read_pairs_f64 = unique_read_pairs as f64;
            let read_pairs_f64 = read_pairs as f64;

            let mut m = 1.0;

            #[allow(non_snake_case)]
            let mut M = 100.0;

            if unique_read_pairs >= read_pairs
                || Self::f(
                    m * unique_read_pairs_f64,
                    unique_read_pairs_f64,
                    read_pairs_f64,
                ) < 0.0
            {
                Err(anyhow!(
                    "Invalid values for pairs and unique pairs: {}, {}",
                    read_pairs,
                    unique_read_pairs
                ))?;
            }

            // find value of M, large enough to act as other side for bisection method
            while Self::f(
                M * unique_read_pairs_f64,
                unique_read_pairs_f64,
                read_pairs_f64,
            ) > 0.0
            {
                M *= 10.0;
            }

            // use bisection method (no more than 40 times) to find solution
            for i in (0..40) {
                let r = (m + M) / 2.0;
                let u = Self::f(
                    r * unique_read_pairs_f64,
                    unique_read_pairs_f64,
                    read_pairs_f64,
                );

                if u == 0.0 {
                    break;
                } else if u > 0.0 {
                    m = r;
                } else if u < 0.0 {
                    M = r;
                }
            }

            Ok(Some((unique_read_pairs_f64 * (m + M) / 2.0) as i64))
        } else {
            Ok(None)
        }
    }

    /**
     * Method that is used in the computation of estimated library size.
     */
    fn f(x: f64, c: f64, n: f64) -> f64 {
        c / x - 1.0 + (-n / x).exp()
    }

    /**
     * Calculates a histogram using the estimateRoi method to estimate the effective yield
     * doing x sequencing for x=1..10.
     */
    pub(crate) fn calculate_roi_histogram(&self) -> Option<Histogram<f64H>> {
        if self.ESTIMATED_LIBRARY_SIZE.is_none() {
            match self.calculate_derived_fields() {
                Ok(_) => {
                    if self.ESTIMATED_LIBRARY_SIZE.is_none() {
                        return None;
                    }
                }
                Err(err) => {
                    return None;
                }
            }
        }

        let unique_pairs = self.READ_PAIRS_EXAMINED - self.READ_PAIR_DUPLICATES;
        let mut histo = Histogram::<f64H>::new();

        for x in (1..=100).map(|e| e as f64) {
            histo.increment(
                x.into(),
                Self::estimate_roi(
                    self.ESTIMATED_LIBRARY_SIZE.clone().unwrap(),
                    x,
                    self.READ_PAIRS_EXAMINED,
                    unique_pairs,
                ),
            );
        }

        histo.set_value_label("CoverageMult".to_string());

        Some(histo)
    }

    /**
     * Estimates the ROI (return on investment) that one would see if a library was sequenced to
     * x higher coverage than the observed coverage.
     *
     * @param estimatedLibrarySize the estimated number of molecules in the library
     * @param x                    the multiple of sequencing to be simulated (i.e. how many X sequencing)
     * @param pairs                the number of pairs observed in the actual sequencing
     * @param uniquePairs          the number of unique pairs observed in the actual sequencing
     * @return a number z <= x that estimates if you had pairs*x as your sequencing then you
     * would observe uniquePairs*z unique pairs.
     */
    fn estimate_roi(estimated_library_size: i64, x: f64, pairs: usize, unique_pairs: usize) -> f64 {
        estimated_library_size as f64
            * (1.0 - (-x * pairs as f64 / estimated_library_size as f64).exp())
            / unique_pairs as f64
    }
}

pub(crate) trait TypeName {
    fn get_name() -> &'static str;
}

impl TypeName for DuplicationMetrics {
    fn get_name() -> &'static str {
        stringify!(DuplicationMetrics)
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
