use std::{fs::File, path::Path};

use super::utils::histogram::Histogram;
use anyhow::Error;

pub(crate) struct MetricsFile<B, H> {
    metrics: Vec<B>,
    histograms: Vec<Histogram<H>>,
}

impl<B, H> MetricsFile<B, H> {
    pub(crate) fn new() -> Self {
        Self {
            metrics: Vec::new(),
            histograms: Vec::new(),
        }
    }

    #[allow(non_upper_case_globals)]
    const serialVersionUID: i64 = 1;
    const MAJOR_HEADER_PREFIX: &'static str = "## ";
    const MINOR_HEADER_PREFIX: &'static str = "# ";
    const SEPARATOR: &'static str = "\t";
    const HISTO_HEADER: &'static str = "## HISTOGRAM\t";
    const METRIC_HEADER: &'static str = "## METRICS CLASS\t";

    /** Adds a bean to the collection of metrics. */
    fn add_metric(&mut self, bean: B) {
        self.metrics.push(bean)
    }

    /** Sets the histogram contained in the metrics file. */
    fn set_histogram(&mut self, histogram: Histogram<H>) {
        if self.histograms.is_empty() {
            self.histograms.push(histogram);
        } else {
            self.histograms
                .get_mut(0)
                .and_then(|v| Some(*v = histogram));
        }
    }

    /** Adds a histogram to the list of histograms in the metrics file. */
    fn add_histogram(&mut self, histogram: Histogram<H>) {
        self.histograms.push(histogram)
    }

    /**
     * Writes out the metrics file to the supplied file. The file is written out
     * headers first, metrics second and histogram third.
     *
     * @param f a File into which to write the metrics
     */
    fn write(&self, f: impl AsRef<Path>) {
        let closure = || -> Result<(), Error> {
            let mut w = File::create(f.as_ref())?;
            self._write(&mut w)?;

            todo!()
        };
    }

    /**
     * Writes out the metrics file to the supplied writer. The file is written out
     * headers first, metrics second and histogram third.
     *
     * @param w a Writer into which to write the metrics
     */
    fn _write(&self, w: &mut File) -> Result<(), Error> {
        // let formatter

        todo!()
    }
}
