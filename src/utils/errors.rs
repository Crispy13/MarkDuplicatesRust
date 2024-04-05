use thiserror::Error as ErrorDerive;

#[derive(ErrorDerive, Debug)]
pub(crate) enum Error {
    #[error("Invalid UTF8 encountered when parsing read name.")]
    InvalidUTF8ReadName,

    #[error("Expected to read {eb} bytes but {rb} bytes. input file={input_file}")]
    FailedDeserializeFromByteFile {
        eb: usize,
        rb: usize,
        input_file: String,
    },
}
