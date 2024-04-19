use std::{fmt::Write, path::PathBuf};

use clap::{Arg, ArgAction, Args, Parser, Subcommand};
use log::LevelFilter;

use crate::{
    hts::{duplicate_scoring_strategy::ScoringStrategy, SortOrder},
    markdup::{
        markduplicates::DuplicateTaggingPolicy,
        utils::optical_duplicate_finder::OpticalDuplicateFinder,
    },
};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[allow(non_snake_case)]
pub(crate) struct Cli {
    /// One or more input SAM, BAM or CRAM files to analyze. Must be coordinate sorted.
    #[arg(long = "INPUT", value_name = "Vec<String>", value_delimiter=',', num_args = 1.. )]
    pub(crate) INPUT: Vec<String>,

    /// The output file to write marked records to
    #[arg(long = "OUTPUT", value_name = "PathBuf", required = true)]
    pub(crate) OUTPUT: PathBuf,

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
        value_name = "usize",
        default_value_t = 8000
    )]
    pub(crate) MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: usize,

    /// The scoring strategy for choosing the non-duplicate among candidates.
    #[arg(
        long = "DUPLICATE_SCORING_STRATEGY",
        value_name = "ScoringStrategy",
        value_enum,
        default_value_t = ScoringStrategy::TotalMappedReferenceLength,
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
    #[arg(long = "DUPLEX_UMI", value_name = "bool", default_value_t = false, action=clap::ArgAction::Set)]
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

    /// enable parameters and behavior specific to flow based reads.
    #[arg(long = "FLOW_MODE", value_name = "bool", default_value_t = false, action=clap::ArgAction::Set)]
    pub(crate) FLOW_MODE: bool,

    /// Use position of the clipping as the end position, when considering duplicates (or use the unclipped end position)
    /// (for this argument, \"read end\" means 3' end).
    #[arg(
        long = "USE_UNPAIRED_CLIPPED_END",
        value_name = "bool",
        default_value_t = false, action=clap::ArgAction::Set
    )]
    pub(crate) USE_UNPAIRED_CLIPPED_END: bool,

    /// Skip first N flows, starting from the read's start, when considering duplicates. Useful for flow based reads where sometimes there
    /// is noise in the first flows
    /// (for this argument, \"read start\" means 5' end).
    #[arg(
        long = "FLOW_SKIP_FIRST_N_FLOWS",
        value_name = "i32",
        default_value_t = 0
    )]
    pub(crate) FLOW_SKIP_FIRST_N_FLOWS: i32,

    /// Treat position of read trimming based on quality as the known end (relevant for flow based reads). Default false - if the read
    /// is trimmed on quality its end is not defined and the read is duplicate of any read starting at the same place.
    #[arg(
        long = "FLOW_Q_IS_KNOWN_END",
        value_name = "bool",
        default_value_t = false, action=clap::ArgAction::Set
    )]
    pub(crate) FLOW_Q_IS_KNOWN_END: bool,

    /// Use specific quality summing strategy for flow based reads. The strategy ensures that the same
    /// (and correct) quality value is used for all bases of the same homopolymer.
    #[arg(
        long = "FLOW_QUALITY_SUM_STRATEGY",
        value_name = "bool",
        default_value_t = false, action=clap::ArgAction::Set
    )]
    pub(crate) FLOW_QUALITY_SUM_STRATEGY: bool,

    /// Threshold for considering a quality value high enough to be included when calculating FLOW_QUALITY_SUM_STRATEGY calculation.
    #[arg(
        long = "FLOW_EFFECTIVE_QUALITY_THRESHOLD",
        value_name = "i32",
        default_value_t = 15
    )]
    pub(crate) FLOW_EFFECTIVE_QUALITY_THRESHOLD: i32,

    /// Maximal difference of the read end position that counted as equal. Useful for flow based
    /// reads where the end position might vary due to sequencing errors.
    /// (for this argument, \"read end\" means 3' end)
    #[arg(
        long = "UNPAIRED_END_UNCERTAINTY",
        value_name = "i32",
        default_value_t = 0
    )]
    pub(crate) UNPAIRED_END_UNCERTAINTY: i32,

    /// Make the end location of single end read be significant when considering duplicates,
    /// in addition to the start location, which is always significant (i.e. require single-ended reads to start and
    /// end on the same position to be considered duplicate)
    /// (for this argument, \"read end\" means 3' end).
    #[arg(
        long = "USE_END_IN_UNPAIRED_READS",
        value_name = "bool",
        default_value_t = false, action=clap::ArgAction::Set
    )]
    pub(crate) USE_END_IN_UNPAIRED_READS: bool,

    /// This number, plus the maximum RAM available to the JVM, determine the memory footprint used by
    /// some of the sorting collections.  If you are running out of memory, try reducing this number.
    #[arg(
        long = "SORTING_COLLECTION_SIZE_RATIO",
        value_name = "f64",
        default_value_t = 0.25
    )]
    pub(crate) SORTING_COLLECTION_SIZE_RATIO: f64,

    /// One or more directories with space available to be used by this program for temporary storage of working files
    #[arg(long = "TMP_DIR", value_name = "Vec<PathBuf>", default_values  = [std::env::temp_dir().into_os_string()])]
    pub(crate) TMP_DIR: Vec<PathBuf>,

    /// Logging level.
    #[arg(long = "LOGGING_LEVEL", value_name = "LevelFilter", default_value_t=LevelFilter::Info)]
    pub(crate) LOGGING_LEVEL: LevelFilter,

    /// If a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG (DS),
    /// indicates the size of the duplicate set. The smallest possible DS value is 2 which occurs when two
    /// reads map to the same portion of the reference only one of which is marked as duplicate. The second
    /// tag, DUPLICATE_SET_INDEX_TAG (DI), represents a unique identifier for the duplicate set to which the
    /// record belongs. This identifier is the index-in-file of the representative read that was selected out
    /// of the duplicate set.
    #[arg(
        long = "TAG_DUPLICATE_SET_MEMBERS",
        value_name = "bool",
        default_value_t = false
        , action=clap::ArgAction::Set
    )]
    pub(crate) TAG_DUPLICATE_SET_MEMBERS: bool,

    /// Max memory to use.
    /// This option is not in Java MarkDuplicates but is here to handle java's `Runtime.getRuntime().maxMemory()` code.
    #[arg(long = "MAX_MEMORY", value_name = "usize", default_value_t=usize::MAX)]
    pub(crate) MAX_MEMORY: usize,

    /// Max memory to use.
    /// This option is not in Java MarkDuplicates but is here to handle java's `Runtime.getRuntime().maxMemory()` code.
    #[arg(long = "READ_NAME_REGEX", value_name = "String", default_value_t=OpticalDuplicateFinder::DEFAULT_READ_NAME_REGEX.to_string())]
    pub(crate) READ_NAME_REGEX: String,

    /// If true remove 'optical' duplicates and other duplicates that appear to have arisen from the
    /// sequencing process instead of the library preparation process, even if REMOVE_DUPLICATES is false.
    /// If REMOVE_DUPLICATES is true, all duplicates are removed and this option is ignored.
    #[arg(long = "REMOVE_SEQUENCING_DUPLICATES", value_name = "bool", default_value_t=false, action=clap::ArgAction::Set)]
    pub(crate) REMOVE_SEQUENCING_DUPLICATES: bool,

    /// If true remove 'optical' duplicates and other duplicates that appear to have arisen from the
    /// sequencing process instead of the library preparation process, even if REMOVE_DUPLICATES is false.
    /// If REMOVE_DUPLICATES is true, all duplicates are removed and this option is ignored.
    #[arg(long = "TAGGING_POLICY", value_name = "DuplicateTaggingPolicy", value_enum, default_value_t=DuplicateTaggingPolicy::DontTag)]
    pub(crate) TAGGING_POLICY: DuplicateTaggingPolicy,

    /// If not null, assume that the input file has this order even if the header says otherwise.
    #[arg(long = "ASSUME_SORTED", value_name = "SortOrder",default_value_t=false, action=clap::ArgAction::Set)]
    pub(crate) ASSUME_SORTED: bool,

    /// If true, assume that the input file is coordinate sorted even if the header says otherwise.
    /// Deprecated, used ASSUME_SORT_ORDER=coordinate instead.
    #[arg(long = "ASSUME_SORT_ORDER", value_name = "SortOrder")]
    pub(crate) ASSUME_SORT_ORDER: Option<SortOrder>,

    /// Comment(s) to include in the output file's header.
    #[arg(long = "COMMENT", value_name = "Vec<String>")]
    pub(crate) COMMENT: Vec<String>,

    /// Value of VN tag of PG record to be created. If not specified, the version will be detected automatically.
    #[arg(long = "PROGRAM_GROUP_VERSION", value_name = "String", default_value_t=String::new())]
    pub(crate) PROGRAM_GROUP_VERSION: String,

    /// Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically.
    #[arg(long = "PROGRAM_GROUP_COMMAND_LINE", value_name = "String", default_value_t=String::new())]
    pub(crate) PROGRAM_GROUP_COMMAND_LINE: String,

    /// Value of PN tag of PG record to be created.
    #[arg(long = "PROGRAM_GROUP_NAME", value_name = "String", default_value_t=String::new())]
    pub(crate) PROGRAM_GROUP_NAME: String,

    /// Clear DT tag from input SAM records. Should be set to false if input SAM doesn't have this tag.  Default true
    #[arg(long = "CLEAR_DT", value_name = "bool", default_value_t = true, action=clap::ArgAction::Set)]
    pub(crate) CLEAR_DT: bool,

    /// SAM tag to uniquely identify the molecule from which a read was derived.  Use of this option requires that
    /// the BARCODE_TAG option be set to a non null value.  Default null.
    #[arg(long = "MOLECULAR_IDENTIFIER_TAG", value_name = "String", default_value_t = String::new())]
    pub(crate) MOLECULAR_IDENTIFIER_TAG: String,

    /// If true do not write duplicates to the output file instead of writing them with appropriate flags set.
    #[arg(long = "REMOVE_DUPLICATES", value_name = "bool", default_value_t=false, action=clap::ArgAction::Set)]
    pub(crate) REMOVE_DUPLICATES: bool,

    /// Add PG tag to each read in a SAM or BAM
    #[arg(long = "ADD_PG_TAG_TO_READS", value_name = "bool", default_value_t=true, action=clap::ArgAction::Set)]
    pub(crate) ADD_PG_TAG_TO_READS: bool,

    /// File to write duplication metrics to
    #[arg(long = "METRICS_FILE", value_name = "PathBuf")]
    pub(crate) METRICS_FILE: PathBuf,

    /// Maximum thread to use.
    #[arg(long = "THREADS", value_name = "usize", default_value_t = 1)]
    pub(crate) THREADS: usize,

    #[arg(skip)]
    pub(crate) commandline: String,
}

impl Cli {
    pub(crate) fn parse_and_post() -> Cli {
        let mut cli = <Self as Parser>::parse();

        cli.after_action();

        cli
    }

    /// do some more works after parsing command line argument.
    pub(crate) fn after_action(&mut self) {
        if self.ASSUME_SORT_ORDER.is_some() || self.ASSUME_SORTED {
            if self.ASSUME_SORT_ORDER.is_none() {
                self.ASSUME_SORT_ORDER = Some(SortOrder::Coordinate);
                self.ASSUME_SORTED = false; // to maintain the "mutex" regarding these two arguments.
            }
        }

        if self.PROGRAM_GROUP_NAME.is_empty() {
            self.PROGRAM_GROUP_NAME = "MarkDuplicates".to_string();
        }

        // make command line string
        self.commandline = std::env::args()
            .into_iter()
            .reduce(|s, x| {
                write!(&mut s, " {}", x).unwrap();
                s
            })
            .unwrap();
    }
}
