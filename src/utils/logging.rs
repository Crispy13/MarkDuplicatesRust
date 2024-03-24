use std::borrow::{Borrow, BorrowMut};
use std::fmt;
use std::sync::{Mutex, OnceLock};

use log::{logger, LevelFilter};
use log4rs::append::console::{ConsoleAppender, Target};
use log4rs::append::file::FileAppender;
use log4rs::config::runtime::ConfigBuilder;
use log4rs::config::{Appender, Config, Logger, Root};
use log4rs::encode::pattern::PatternEncoder;

struct RuntimeLogConfig {
    handle: log4rs::Handle,
    root_level: LevelFilter,
}

const LOG_DEST_STDERR:&'static str = "stderr";
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
        .encoder(Box::new(PatternEncoder::new(
            "[{d(%Y-%m-%d %H:%M:%S)}] {T} {t} {l}>> {m}{n}",
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
    log: String,
    n: u32,
    verb: &'static str,
    noun: &'static str,
}

impl ProgressLogger {
    

    pub(crate) fn new(log: String, n: u32, verb: &'static str, noun: &'static str) -> Self {
        Self { log, n, verb, noun }
    }

    pub(crate) fn record(&self, chrom:impl fmt::Display, pos:i32) {
        const level: &str="info";
        macro_rules! do_log {
            ($level:tt) => {
                log::$level!(target: self.log.as_str(), "chr={} pos={} {} {}.", chrom, pos, self.noun, self.verb);
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



#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_add_logger() {
        init_global_logger(LevelFilter::Info);

        log::info!(target:"A", "Hey A1 info!");

        add_logger_and_set_config(
            Logger::builder().additive(true).build("A", LevelFilter::Debug)
        );

        log::debug!(target:"A", "Hey A2 debug!");
        log::info!(target:"A", "Hey A2 info!");

        
    }

    #[test]
    fn progress_logger() {
        init_global_logger(LevelFilter::Debug);

        let pr1 = ProgressLogger::new("A".to_owned(), 10000, "compared", "ReadEnds to Keeper");
        let pr2 = ProgressLogger::new("B".to_owned(), 1000, "compared", "ReadEnds to Keeper");

        pr1.record("chr1", 121212,);
        pr2.record("chr2", 34343434, );

    }
}