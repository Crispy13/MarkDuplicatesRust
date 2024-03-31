use thiserror::Error as ErrorDerive;


#[derive(ErrorDerive, Debug)]
pub(crate) enum Error {
    #[error("Invalid UTF8 encountered when parsing read name.")]
    InvalidUTF8ReadName
}