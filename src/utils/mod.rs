pub(crate) mod errors;
pub(crate) mod logging;

use tikv_jemalloc_ctl;

use std::{
    fmt::{Debug, Display},
    hash::{BuildHasher, BuildHasherDefault, DefaultHasher, Hash, Hasher},
    marker::PhantomData,
    path::PathBuf,
    sync::{atomic::AtomicUsize, OnceLock},
};

#[derive(Debug)]
pub(crate) struct Error<C: Display + Debug> {
    kind: ErrorKind,
    cxt: C,
}

#[derive(Debug)]
pub(crate) enum ErrorKind {
    InvalidUTF8,
}

impl<C: Display + Debug> std::error::Error for Error<C> {}

impl<C: Display + Debug> Display for Error<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(self, f)
    }
}

pub(crate) trait PathExt<C: Display + Debug = Self> {
    fn try_to_str(&self) -> Result<&str, Error<C>>;
}

impl PathExt<String> for PathBuf {
    fn try_to_str(&self) -> Result<&str, Error<String>> {
        match self.to_str() {
            Some(s) => Ok(s),
            None => Err(Error {
                kind: ErrorKind::InvalidUTF8,
                cxt: format!("{:?}", self),
            }),
        }
    }
}

type CommonHasherBuilder = DefaultHasher;
#[derive(Default)]
pub(crate) struct CommonHasher<H: Default + Hasher = DefaultHasher> {
    // dummy: PhantomData<H>,
    hasher_builder: BuildHasherDefault<H>,
}

pub(crate) fn hash_code<T: Hash>(value: T) -> u64 {
    static hasher_builder: OnceLock<CommonHasher> = OnceLock::new();
    let mut hasher = hasher_builder
        .get_or_init(|| CommonHasher::default())
        .hasher_builder
        .build_hasher();

    value.hash(&mut hasher);
    hasher.finish()
}

impl CommonHasher {
    pub(crate) fn get_hash_code<T: Hash>(value: T) -> u64 {
        let hasher_builder = Self::default().hasher_builder;

        let mut hasher = hasher_builder.build_hasher();

        value.hash(&mut hasher);
        hasher.finish()
    }
}

pub(crate) struct MemoryStatsRecord {
    pub(crate) allocated: usize,
    pub(crate) resident: usize,
}

pub(crate) struct MemoryStats {
    allocated: AtomicUsize,
    resident: AtomicUsize,
}

impl MemoryStats {
    /// Get memory usage, `(allocated, resident)`.
    pub(crate) fn record(&self) -> MemoryStatsRecord {
        let r = match Self::get_allocated_and_resident_mem_of_app() {
            Ok(r) => r,
            Err(err) => {
                // if failed to get mem stats, just load previous value.
                return MemoryStatsRecord {
                    allocated: self.allocated.load(std::sync::atomic::Ordering::Relaxed),
                    resident: self.resident.load(std::sync::atomic::Ordering::Relaxed),
                };
            }
        };

        // save mem stats
        self.allocated
            .store(r.allocated, std::sync::atomic::Ordering::Relaxed);
        self.resident
            .store(r.resident, std::sync::atomic::Ordering::Relaxed);

        r
    }

    fn get_allocated_and_resident_mem_of_app() -> Result<MemoryStatsRecord, anyhow::Error> {
        let e = tikv_jemalloc_ctl::epoch::mib().or_else(|err| Err(anyhow::anyhow!("{:?}", err)))?;

        let allocated = tikv_jemalloc_ctl::stats::allocated::mib()
            .or_else(|err| Err(anyhow::anyhow!("{:?}", err)))?;

        let resident = tikv_jemalloc_ctl::stats::resident::mib()
            .or_else(|err| Err(anyhow::anyhow!("{:?}", err)))?;

        e.advance()
            .or_else(|err| Err(anyhow::anyhow!("{:?}", err)))?;

        Ok(MemoryStatsRecord {
            allocated: allocated
                .read()
                .or_else(|err| Err(anyhow::anyhow!("{:?}", err)))?,
            resident: resident
                .read()
                .or_else(|err| Err(anyhow::anyhow!("{:?}", err)))?,
        })
    }
}

pub(crate) fn mem_stats() -> &'static MemoryStats {
    #[allow(non_upper_case_globals)]
    static mem_stats: OnceLock<MemoryStats> = OnceLock::new();

    mem_stats.get_or_init(|| MemoryStats {
        allocated: AtomicUsize::new(0),
        resident: AtomicUsize::new(0),
    })
}

pub(crate) fn get_resident_mem_of_app() -> usize {
    let e = tikv_jemalloc_ctl::epoch::mib().unwrap();
    let resident = tikv_jemalloc_ctl::stats::resident::mib().unwrap();

    e.advance().unwrap();
    resident.read().unwrap()
}

pub(crate) fn get_allocated_mem_of_app() -> usize {
    let e = tikv_jemalloc_ctl::epoch::mib().unwrap();
    let allocated = tikv_jemalloc_ctl::stats::allocated::mib().unwrap();

    e.advance().unwrap();
    allocated.read().unwrap()
}

pub(crate) fn get_allocated_and_resident_mem_of_app() -> (usize, usize) {
    let e = tikv_jemalloc_ctl::epoch::mib().unwrap();

    let allocated = tikv_jemalloc_ctl::stats::allocated::mib().unwrap();

    let resident = tikv_jemalloc_ctl::stats::resident::mib().unwrap();

    e.advance().unwrap();

    (allocated.read().unwrap(), resident.read().unwrap())
}

#[inline]
pub(crate) fn human_readable_byte_count(bytes: usize) -> String {
    if bytes < 1024 {
        return format!("{} B", bytes);
    }

    let exp = ((bytes as f64).ln() / 1024_f64.ln()) as usize;

    format!(
        "{:.1} {}B",
        bytes as f64 / 1024_usize.pow(exp as u32) as f64,
        "kMGTPE".get((exp - 1)..(exp)).unwrap()
    )
}

#[cfg(test)]
mod test {
    #[global_allocator]
    static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

    use std::hash::DefaultHasher;

    use crate::utils::{
        get_allocated_and_resident_mem_of_app, get_allocated_mem_of_app, human_readable_byte_count,
    };

    use super::{get_resident_mem_of_app, CommonHasher};

    #[test]
    fn get_mem_sizes() {
        let (a, r) = get_allocated_and_resident_mem_of_app();

        println!("allocated:{}", a);
        println!("resident:{}", r);

        let b = Vec::<i8>::with_capacity(1024_usize.pow(2));

        let c = Vec::<i16>::with_capacity(1024_usize.pow(2));

        let (a, r) = get_allocated_and_resident_mem_of_app();

        println!("allocated:{}", a);
        println!("resident:{}", r);

        drop(b);

        let (a, r) = get_allocated_and_resident_mem_of_app();

        println!("allocated:{}", a);
        println!("resident:{}", r);

        println!("{}", human_readable_byte_count(r));
    }

    #[test]
    fn hasher_check() {
        let ch1 = CommonHasher::get_hash_code(-112);
        let ch2 = CommonHasher::get_hash_code(-112);

        println!("{} {}", ch1, ch2);
    }
}
