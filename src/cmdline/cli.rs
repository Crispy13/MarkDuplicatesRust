use clap::{ArgAction, Args, Parser, Subcommand};

use crate::hts::duplicate_scoring_strategy::ScoringStrategy;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
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
        default_value_t = 8000,
    )]
    pub(crate) MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: i32,

    /// The scoring strategy for choosing the non-duplicate among candidates.
    #[arg(
        long = "DUPLICATE_SCORING_STRATEGY",
        value_name = "ScoringStrategy",
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
    #[arg(
        long = "DUPLEX_UMI",
        value_name = "bool",
        default_value_t = false,
    )]
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
    #[arg(
        long = "INPUT",
        value_name = "String",
    )]
    pub(crate) INPUT: Vec<String>,
}