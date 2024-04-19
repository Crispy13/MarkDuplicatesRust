use std::{
    fmt::Display,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use crate::markdup::utils::library_id_generator::TypeName;

use super::utils::histogram::Histogram;
use anyhow::Error;
use rust_htslib::bam::Header;

pub(crate) struct MetricsFile<B, H> {
    metrics: Vec<B>,
    histograms: Vec<Histogram<H>>,

    headers: Vec<MetricsHeader>,
}

impl<B, H> MetricsFile<B, H>
where
    B: TypeName,
{
    pub(crate) const MAJOR_HEADER_PREFIX: &'static str = "## ";
    pub(crate) const MINOR_HEADER_PREFIX: &'static str = "# ";

    pub(crate) fn new() -> Self {
        Self {
            metrics: Vec::new(),
            histograms: Vec::new(),
            headers: Vec::new(),
        }
    }

    #[allow(non_upper_case_globals)]
    const serialVersionUID: i64 = 1;
    const SEPARATOR: &'static str = "\t";
    const HISTO_HEADER: &'static str = "## HISTOGRAM\t";
    const METRIC_HEADER: &'static str = "## METRICS CLASS\t";

    /** Adds a bean to the collection of metrics. */
    pub(crate) fn add_metric(&mut self, bean: B) {
        self.metrics.push(bean)
    }

    /** Sets the histogram contained in the metrics file. */
    pub(crate) fn set_histogram(&mut self, histogram: Option<Histogram<H>>) {
        if self.histograms.is_empty() {
            if let Some(h) = histogram {
                self.histograms.push(h);
            }
        } else {
            self.histograms
                .get_mut(0)
                .and_then(|v| Some(*v = histogram.unwrap()));
        }
    }

    /** Adds a histogram to the list of histograms in the metrics file. */
    pub(crate) fn add_histogram(&mut self, histogram: Histogram<H>) {
        self.histograms.push(histogram)
    }

    /**
     * Writes out the metrics file to the supplied file. The file is written out
     * headers first, metrics second and histogram third.
     *
     * @param f a File into which to write the metrics
     */
    pub(crate) fn write(&self, f: impl AsRef<Path>) {
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
        // let formatter = FormatUtil // maybe formatter is not needed in rust. use Display trait instead.
        let mut out = BufWriter::new(w);
        self.print_header(&mut out);

        write!(&mut out, "\n")?;

        todo!()
    }

    /** Prints the headers into the provided PrintWriter. */
    fn print_header(&self, out: impl Write) -> Result<(), std::io::Error> {
        self.headers
            .iter()
            .map(|h| {
                write!(
                    out,
                    "{}{}\n{}{}\n",
                    Self::MAJOR_HEADER_PREFIX,
                    h.variant_name(),
                    Self::MINOR_HEADER_PREFIX,
                    h,
                )
            })
            .collect::<Result<(), std::io::Error>>()
    }

    /** Prints each of the metrics entries into the provided PrintWriter. */
    fn print_bean_metrics(&self, out: impl Write) -> Result<(), Error> {
        if self.metrics.is_empty() {
            return Ok(());
        }

        // Write out a header row with the type of the metric class
        write!(out, "{}{}\n", Self::METRIC_HEADER, self.get_bean_type_name())?;

        todo!()
    }

    /** Gets the type of the metrics bean being used. */
    fn get_bean_type_name(&self) -> &'static str {
        if self.metrics.is_empty() {
            return "null";
        }

        B::get_name()
    }

    /** Adds a header to the collection of metrics. */
    pub(crate) fn add_header(&self, h: MetricsHeader) {
        self.headers.push(h)
    }
}

#[derive(Clone)]
pub(crate) enum MetricsHeader {
    StringHeader(String),
    VersionHeader,
}

impl MetricsHeader {
    pub(crate) fn variant_name(&self) -> &'static str {
        match self {
            MetricsHeader::StringHeader(_) => "StringHeader",
            MetricsHeader::VersionHeader => "VersionHeader",
        }
    }
}

impl Display for MetricsHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            MetricsHeader::StringHeader(v) => v,
            MetricsHeader::VersionHeader => self.variant_name(),
        };
        write!(f, "{}", s)
    }
}
