use clap::Parser;
use cmdline::cli::Cli;
use markdup::markduplicates::MarkDuplicates;
use utils::logging::init_global_logger;

#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

mod cmdline;
mod hts;
mod markdup;
mod utils;

#[cfg(test)]
mod tests;

fn main() {
    let mut cli = Cli::parse();
    cli.after_action();

    init_global_logger(cli.LOGGING_LEVEL);

    let mut md = MarkDuplicates::new(cli);

    md.do_work();
}
