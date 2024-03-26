pub(crate) mod logging;

use std::{
    fmt::{Debug, Display},
    hash::{BuildHasher, BuildHasherDefault, DefaultHasher, Hash, Hasher},
    marker::PhantomData,
    path::PathBuf,
    sync::OnceLock,
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

pub(crate) fn hash_code<T:Hash>(value: T) -> u64 {
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

#[cfg(test)]
mod test {
    use std::hash::DefaultHasher;

    use super::CommonHasher;

    #[test]
    fn hasher_check() {
        let ch1 = CommonHasher::get_hash_code(-112);
        let ch2 = CommonHasher::get_hash_code(-112);

        println!("{} {}", ch1, ch2);
    }
}
