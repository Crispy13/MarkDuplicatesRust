use std::borrow::{Borrow, BorrowMut};
use std::fmt;
use std::process::exit;
use std::sync::{Mutex, OnceLock};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use log::{logger, LevelFilter};
use log4rs::append::console::{ConsoleAppender, Target};
use log4rs::append::file::FileAppender;
use log4rs::config::runtime::ConfigBuilder;
use log4rs::config::{Appender, Config, Logger, Root};
use log4rs::encode::pattern::PatternEncoder;
use rust_htslib::bam::Record;

use bio_types::genome::AbstractInterval;

struct RuntimeLogConfig {
    handle: log4rs::Handle,
    root_level: LevelFilter,
}

const LOG_DEST_STDERR: &'static str = "stderr";
pub(crate) fn add_logger_and_set_config(logger: Logger) {
    let mut logger_vec = logger_vec().lock().unwrap();
    logger_vec.push(logger);

    let runtime_log_config = RUNTIME_LOG_CONFIG.get().unwrap().lock().unwrap();
    let mut config_builder = default_configbuilder();

    for logger in logger_vec.iter() {
        config_builder = config_builder.logger(logger.clone())
    }

    // let stderr_appender = ConsoleAppender::builder().encoder(
    //     Box::new(PatternEncoder::new("[{d(%Y-%m-%d %H:%M:%S)}] stderr2 {t} {l}>> {m}{n}"))
    // ).target(Target::Stderr).build();
    // let file_appender = FileAppender::builder().append(false).build("debug.log").unwrap();

    drop(logger_vec); // release lock

    let config = config_builder
        // .appender(Appender::builder().build("stderr2", Box::new(stderr_appender)))
        .build(
            Root::builder()
                .appender(LOG_DEST_STDERR)
                .build(runtime_log_config.root_level.clone()),
        )
        .unwrap();

    runtime_log_config.handle.set_config(config);
}

fn logger_vec() -> &'static Mutex<Vec<Logger>> {
    static logger_vec: OnceLock<Mutex<Vec<Logger>>> = OnceLock::new();
    logger_vec.get_or_init(|| Mutex::new(vec![]))
}

pub(crate) fn default_configbuilder() -> ConfigBuilder {
    let stderr = ConsoleAppender::builder()
        .target(log4rs::append::console::Target::Stderr)
        .encoder(Box::new(log4rs::encode::pattern::PatternEncoder::new(
            "{l:<7} {d(%Y-%m-%d %H:%M:%S)}     {t}  {m}{n}",
        )))
        .build();

    let config = Config::builder()
        // .appender(Appender::builder().build("file", Box::new(file_appender)))
        .appender(Appender::builder().build(LOG_DEST_STDERR, Box::new(stderr)));
    // .logger(Logger::builder().build("app::backend::db", LevelFilter::Info))
    // .logger(Logger::builder()
    //     .appender("requests")
    //     .additive(false)
    //     .build("app::requests", LevelFilter::Info))

    config
}

static RUNTIME_LOG_CONFIG: OnceLock<Mutex<RuntimeLogConfig>> = OnceLock::new();

// pub(crate) fn get_global_log_handle() -> &'static log4rs::Handle {
//     &RUNTIME_LOG_CONFIG.get().unwrap().lock().unwrap().handle
// }

pub(crate) fn init_global_logger(level: LevelFilter) {
    let configbuilder = default_configbuilder();
    let config = configbuilder
        .build(Root::builder().appender(LOG_DEST_STDERR).build(level))
        .unwrap();

    let handle = log4rs::init_config(config).unwrap();

    RUNTIME_LOG_CONFIG.get_or_init(|| {
        Mutex::new(RuntimeLogConfig {
            handle,
            root_level: level,
        })
    });
}

pub(crate) struct ProgressLogger {
    log: &'static str,
    n: usize,
    verb: &'static str,
    noun: &'static str,

    last_chrom: String,
    last_pos: i64,
    last_read_name: String,

    start_time: Instant,
    last_start_time: Option<u64>,
    count_non_increasing: i64,

    processed: usize,
}

impl ProgressLogger {
    pub(crate) fn new(log: &'static str, n: usize, verb: &'static str, noun: &'static str) -> Self {
        Self {
            log,
            n,
            verb,
            noun,
            start_time: Instant::now(),
            processed: 0,
            last_start_time: None,
            last_chrom: String::new(),
            last_pos: 0,
            last_read_name: String::new(),
            count_non_increasing: 0,
        }
    }

    pub(crate) fn record(&mut self, rec: &Record) -> bool {
        let read_name = std::str::from_utf8(rec.qname()).unwrap_or_else(|err| {
            eprintln!("{}", err);
            panic!("Invalid UTF-8 encountered in read name.");
        });

        if "*".eq(read_name) {
            self.check_and_then_record("", 0, read_name)
        } else {
            self.check_and_then_record(rec.contig(), rec.pos(), read_name)
        }
    }

    fn check_and_then_record(&mut self, chrom: &str, pos: i64, rname: &str) -> bool {
        if !chrom.is_empty() && chrom.eq(&self.last_chrom) && pos < self.last_pos {
            self.count_non_increasing += 1;
        } else {
            self.last_chrom.clear();
            self.last_chrom.push_str(chrom);
        }

        self.last_pos = pos;
        self.last_read_name.clear();
        self.last_read_name.push_str(rname);

        if self.last_start_time.is_none() {
            self.last_start_time = Some(self.start_time.elapsed().as_secs());
        }

        self.processed += 1;
        if self.processed % self.n == 0 {
            self.__record();
            true
        } else {
            false
        }
    }

    fn __record(&mut self) {
        let seconds = self.start_time.elapsed().as_secs();

        let last_period_seconds = seconds - self.last_start_time.as_ref().unwrap();

        self.last_start_time = Some(seconds);

        let elapsed = format_elapsed_time(seconds);
        let period = pad(&last_period_seconds.to_string(), 4);
        let processed = pad(&self.processed.to_string(), 13);

        let read_info = if self.last_chrom.is_empty() {
            "*/*".to_string()
        } else {
            format!("{}:{}", self.last_chrom, self.last_pos)
        };

        let rn_info = if !self.last_read_name.is_empty() && self.count_non_increasing > 1000 {
            format!(".  Last read name: {}", self.last_read_name)
        } else {
            String::new()
        };

        let n = if self.processed % self.n == 0 {
            self.n
        } else {
            self.processed % self.n
        };

        log::info!(target: &self.log,
            "{} {} {} .  Elapsed time: {}s.  Time for last {}: {}s.  Last read position: {}{}",
            self.verb,
            processed,
            self.noun,
            elapsed,
            n,
            period,
            read_info,
            rn_info
        );
    }

    pub(crate) fn record_with_chrom_pos(&mut self, chrom: &str, pos: i64) {
        self.check_and_then_record(chrom, pos, "");
    }

    fn _record(&self, chrom: impl fmt::Display, pos: i32) {
        const level: &str = "info";
        macro_rules! do_log {
            ($level:tt) => {
                log::$level!(target: self.log, "chr={} pos={} {} {}.", chrom, pos, self.noun, self.verb)
            };
        }

        macro_rules! match_level_and_log {
            ($($level:tt),+) => {
                match level {
                    $(
                    stringify!($level) => {
                        do_log!($level)
                    }
                    ),+

                    _ => unreachable!(),
                }

            };
        }

        match_level_and_log!(debug, info, warn, error);
    }
}

fn format_elapsed_time(seconds: u64) -> String {
    let s = seconds % 60;
    let all_minutes = seconds / 60;
    let m = all_minutes % 60;
    let h = all_minutes / 60;

    format!("{:0>2}:{:0>2}:{:0>2}", h, m, s)
}

fn pad(s: &str, length: usize) -> String {
    let mut in_builder = String::with_capacity(length.max(s.len()));

    while in_builder.len() < length - s.len() {
        in_builder.push_str(" ");
    }

    in_builder.push_str(s);

    in_builder
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_add_logger() {
        init_global_logger(LevelFilter::Info);

        log::info!(target:"A", "Hey A1 info!");

        add_logger_and_set_config(
            Logger::builder()
                .additive(true)
                .build("A", LevelFilter::Debug),
        );

        log::debug!(target:"A", "Hey A2 debug!");
        log::info!(target:"A", "Hey A2 info!");
    }

    #[test]
    fn progress_logger() {
        init_global_logger(LevelFilter::Debug);

        let pr1 = ProgressLogger::new("A", 10000, "compared", "ReadEnds to Keeper");
        let pr2 = ProgressLogger::new("B", 1000, "compared", "ReadEnds to Keeper");

        pr1._record("chr1", 121212);
        pr2._record("chr2", 34343434);
    }
}
