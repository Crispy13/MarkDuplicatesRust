use rust_htslib::bam::{extra_ext::RecordExt, Record};

use crate::utils::hash_code;

pub(crate) struct EstimateLibraryComplexity {}

impl EstimateLibraryComplexity {
    pub(crate) fn get_read_barcode_value(rec: &Record, tag: &str) -> u64 {
        if tag.is_empty() {
            return 0;
        }

        match rec.get_str_aux(tag.as_bytes()) {
            Ok(a) => hash_code(a),
            Err(_) => 0,
        }
    }

    fn get_read_one_barcode_value() {}
}
