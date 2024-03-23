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

}