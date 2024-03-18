use clap::{ArgAction, Args, Parser, Subcommand};

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
}