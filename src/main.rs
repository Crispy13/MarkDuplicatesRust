use clap::Parser;
use cmdline::cli::Cli;
use markdup::markduplicates::MarkDuplicates;

#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

mod cmdline;
mod markdup;
mod hts;
mod utils;
mod tests;

fn main() {
    let cli = Cli::parse();

    let mut md = MarkDuplicates::new(cli);

    md.do_work();
}
