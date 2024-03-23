use std::{
    fmt::{Debug, Display},
    path::PathBuf,
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
