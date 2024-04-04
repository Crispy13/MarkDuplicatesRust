use std::path::PathBuf;

use clap::{ArgAction, Args, Parser, Subcommand};

use crate::hts::duplicate_scoring_strategy::ScoringStrategy;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[allow(non_snake_case)]
pub(crate) struct Cli {
    /// The program record ID for the @PG record(s) created by this program. Set to null to disable
    /// PG record creation.  This string may have a suffix appended to avoid collision with other
    /// program record IDs.
    #[arg(
        // short=  StandardOptionDefinitions::ASSUME_SORT_ORDER_SHORT_NAME, 
        long = "PROGRAM_RECORD_ID",
        value_name = "String",
        default_value = "MarkDuplicates",
    )]
    pub(crate) PROGRAM_RECORD_ID: Option<String>,

    /// Maximum number of file handles to keep open when spilling read ends to disk.
    /// Set this number a little lower than the per-process maximum number of file that may be open.
    /// This number can be found by executing the 'ulimit -n' command on a Unix system.
    #[arg(
        long = "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP",
        value_name = "i32",
        default_value_t = 8000
    )]
    pub(crate) MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: i32,

    /// The scoring strategy for choosing the non-duplicate among candidates.
    #[arg(
        long = "DUPLICATE_SCORING_STRATEGY",
        value_name = "ScoringStrategy",
        value_enum,
        default_value_t = ScoringStrategy::TOTAL_MAPPED_REFERENCE_LENGTH,
    )]
    pub(crate) DUPLICATE_SCORING_STRATEGY: ScoringStrategy,

    /// Barcode SAM tag (ex. BC for 10X Genomics)
    #[arg(
        long = "BARCODE_TAG",
        value_name = "String",
        default_value_t = String::new(),
    )]
    pub(crate) BARCODE_TAG: String,

    /// Treat UMIs as being duplex stranded.  This option requires that the UMI consist of two equal length
    /// strings that are separated by a hyphen (e.g. 'ATC-GTC'). Reads are considered duplicates if, in addition to standard
    /// definition, have identical normalized UMIs.  A UMI from the 'bottom' strand is normalized by swapping its content
    /// around the hyphen (eg. ATC-GTC becomes GTC-ATC).  A UMI from the 'top' strand is already normalized as it is.
    /// Both reads from a read pair considered top strand if the read 1 unclipped 5' coordinate is less than the read
    /// 2 unclipped 5' coordinate. All chimeric reads and read fragments are treated as having come from the top strand.
    /// With this option is it required that the BARCODE_TAG hold non-normalized UMIs. Default false.
    #[arg(long = "DUPLEX_UMI", value_name = "bool", default_value_t = false)]
    pub(crate) DUPLEX_UMI: bool,

    /// Read one barcode SAM tag (ex. BX for 10X Genomics)
    #[arg(
        long = "READ_ONE_BARCODE_TAG",
        value_name = "String",
        default_value_t = String::new(),
    )]
    pub(crate) READ_ONE_BARCODE_TAG: String,

    /// Read two barcode SAM tag (ex. BX for 10X Genomics)
    #[arg(
        long = "READ_TWO_BARCODE_TAG",
        value_name = "String",
        default_value_t = String::new(),
    )]
    pub(crate) READ_TWO_BARCODE_TAG: String,

    /// One or more input SAM, BAM or CRAM files to analyze. Must be coordinate sorted.
    #[arg(long = "INPUT", value_name = "String", required=true)]
    pub(crate) INPUT: Vec<String>,

    /// enable parameters and behavior specific to flow based reads.
    #[arg(long = "FLOW_MODE", value_name = "bool", default_value_t=false)]
    pub(crate) FLOW_MODE: bool,

    /// Use position of the clipping as the end position, when considering duplicates (or use the unclipped end position) 
    /// (for this argument, \"read end\" means 3' end).
    #[arg(long = "USE_UNPAIRED_CLIPPED_END", value_name = "bool", default_value_t=false)]
    pub(crate) USE_UNPAIRED_CLIPPED_END: bool,

    /// Skip first N flows, starting from the read's start, when considering duplicates. Useful for flow based reads where sometimes there 
    /// is noise in the first flows 
    /// (for this argument, \"read start\" means 5' end).
    #[arg(long = "FLOW_SKIP_FIRST_N_FLOWS", value_name = "i32", default_value_t=0)]
    pub(crate) FLOW_SKIP_FIRST_N_FLOWS: i32,

    /// Treat position of read trimming based on quality as the known end (relevant for flow based reads). Default false - if the read 
    /// is trimmed on quality its end is not defined and the read is duplicate of any read starting at the same place.
    #[arg(long = "FLOW_Q_IS_KNOWN_END", value_name = "bool", default_value_t=false)]
    pub(crate) FLOW_Q_IS_KNOWN_END: bool,

    /// Use specific quality summing strategy for flow based reads. The strategy ensures that the same 
    /// (and correct) quality value is used for all bases of the same homopolymer.
    #[arg(long = "FLOW_QUALITY_SUM_STRATEGY", value_name = "bool", default_value_t=false)]
    pub(crate) FLOW_QUALITY_SUM_STRATEGY: bool,

    /// Threshold for considering a quality value high enough to be included when calculating FLOW_QUALITY_SUM_STRATEGY calculation.
    #[arg(long = "FLOW_EFFECTIVE_QUALITY_THRESHOLD", value_name = "i32", default_value_t=15)]
    pub(crate) FLOW_EFFECTIVE_QUALITY_THRESHOLD: i32,

    /// Maximal difference of the read end position that counted as equal. Useful for flow based 
    /// reads where the end position might vary due to sequencing errors. 
    /// (for this argument, \"read end\" means 3' end)
    #[arg(long = "UNPAIRED_END_UNCERTAINTY", value_name = "i32", default_value_t=0)]
    pub(crate) UNPAIRED_END_UNCERTAINTY: i32,

    /// Make the end location of single end read be significant when considering duplicates, 
    /// in addition to the start location, which is always significant (i.e. require single-ended reads to start and
    /// end on the same position to be considered duplicate) 
    /// (for this argument, \"read end\" means 3' end).
    #[arg(long = "USE_END_IN_UNPAIRED_READS", value_name = "bool", default_value_t=false)]
    pub(crate) USE_END_IN_UNPAIRED_READS: bool,

    /// This number, plus the maximum RAM available to the JVM, determine the memory footprint used by 
    /// some of the sorting collections.  If you are running out of memory, try reducing this number.
    #[arg(long = "SORTING_COLLECTION_SIZE_RATIO", value_name = "f64", default_value_t=0.25)]
    pub(crate) SORTING_COLLECTION_SIZE_RATIO: f64,

    /// One or more directories with space available to be used by this program for temporary storage of working files
    #[arg(long = "TMP_DIR", value_name = "Vec<PathBuf>", default_value="[]")]
    pub(crate) TMP_DIR: Vec<PathBuf>,


}
