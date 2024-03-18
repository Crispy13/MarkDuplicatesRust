pub(crate) enum SAMTag {
    PG
}

impl SAMTag {
    pub(crate) fn name(&self) -> &'static str {
        match self {
            SAMTag::PG => "PG",
        }
    }
}