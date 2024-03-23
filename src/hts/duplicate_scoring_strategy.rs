use std::{fmt::Display, sync::OnceLock};

use clap::ValueEnum;
use rust_htslib::bam::{ext::BamRecordExtensions, Record};

use super::utils::murmur3::Murmur3;

type Error = Box<dyn std::error::Error + Send + Sync>;
#[derive(ValueEnum, Clone, Debug, Copy)]
pub(crate) enum ScoringStrategy {
    SUM_OF_BASE_QUALITIES,
    TOTAL_MAPPED_REFERENCE_LENGTH,
    RANDOM,
}

impl Display for ScoringStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, f)
    }
}

/** An enum to use for storing temporary attributes on SAMRecords. */
pub(crate) enum Attr {
    DuplicateScore,
}

pub(crate) struct DuplicateScoringStrategy {}

impl DuplicateScoringStrategy {
    fn new() {}

    fn hasher() -> &'static Murmur3 {
        static murmur3: OnceLock<Murmur3> = OnceLock::new();
        murmur3.get_or_init(|| Murmur3::new(1))
    }

    /** Calculates a score for the read which is the sum of scores over Q15. */
    fn get_sum_of_base_qualities(rec: &Record) -> i32 {
        let mut score = 0_i32;
        rec.qual().iter().copied().for_each(|b| {
            if b >= 15 {
                score += b as i32;
            }
        });

        score
    }

    /// .
    ///
    /// # Errors
    ///
    /// This function will return an error if it fails to get mate cigar.
    pub(crate) fn compute_duplicate_score(
        rec: &Record,
        scoring_strategy: ScoringStrategy,
        assume_mate_cigar: bool,
    ) -> Result<i16, Error> {
        // the original java code can load previously saved score.
        // But yet, I have not found to do this in rust-htslib.

        let mut score = 0_i16;

        match scoring_strategy {
            ScoringStrategy::SUM_OF_BASE_QUALITIES => {
                // two (very) long reads worth of high-quality bases can go over Short.MAX_VALUE/2
                // and risk overflow.
                score += Self::get_sum_of_base_qualities(rec).min(i16::MAX as i32 / 2) as i16;
            }
            ScoringStrategy::TOTAL_MAPPED_REFERENCE_LENGTH => {
                if !rec.is_unmapped() {
                    // no need to remember the score since this scoring mechanism is symmetric
                    score = rec
                        .cigar()
                        .get_reference_length()
                        .min(i16::MAX as usize / 2) as i16;
                }

                if assume_mate_cigar && rec.is_paired() && !rec.is_mate_unmapped() {
                    score += rec
                        .mate_cigar()?
                        .get_reference_length()
                        .min(i16::MAX as usize / 2) as i16;
                }
            }
            ScoringStrategy::RANDOM => {
                // start with a random number between Short.MIN_VALUE/4 and Short.MAX_VALUE/4
                score += (Self::hasher().hash_bytes(rec.qname()) & 0b11_1111_1111_1111) as i16;
                // subtract Short.MIN_VALUE/4 from it to end up with a number between
                // 0 and Short.MAX_VALUE/2. This number can be then discounted in case the read is
                // not passing filters. We need to stay far from overflow so that when we add the two
                // scores from the two read mates we do not overflow since that could cause us to chose a
                // failing read-pair instead of a passing one.
                score -= i16::MAX / 4;
            }
        }

        // make sure that filter-failing records are heavily discounted. (the discount can happen twice, once
        // for each mate, so need to make sure we do not subtract more than Short.MIN_VALUE overall.)
        score += if rec.is_quality_check_failed() {
            i16::MAX / 2
        } else {
            0
        };

        // TODO: save `score` to record.

        Ok(score)
    }
}
