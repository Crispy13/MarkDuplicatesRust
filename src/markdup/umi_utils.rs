use bio_types::genome::AbstractInterval;
use std::sync::OnceLock;

use anyhow::{anyhow, Error};
use regex::Regex;
use rust_htslib::bam::{
    ext::BamRecordExtensions, extra_ext::RecordExt, record::Aux, Record
};
pub(crate) struct UmiUtil {}

impl UmiUtil {
    const DUPLEX_UMI_DELIMITER: &'static str = "-";
    const CONTIG_SEPARATOR: &'static str = ":";
    const UMI_NAME_SEPARATOR: &'static str = "/";
    const TOP_STRAND_DUPLEX: &'static str = "/A";
    const BOTTOM_STRAND_DUPLEX: &'static str = "/B";

    fn allowed_umi() -> &'static Regex {
        static rep: OnceLock<Regex> = OnceLock::new();
        rep.get_or_init(|| Regex::new(r"^[ATCGNatcgn-]*$").unwrap())
    }

    /**
     * Creates a top-strand normalized duplex UMI.
     * Single stranded UMIs are by definition already top-strand normalized, they require no transformation.
     * Duplex UMIs that come from a top strand read are also by definition, top-strand normalized. A duplex
     * UMI from a bottom strand can be normalized to be identical to the read from its corresponding top strand
     * by swapping the content of the UMI around the "-" found in duplex UMIs.  For example, a bottom strand
     * duplex UMI reading ATC-CGG when top-strand normalized will read CGG-ATC.
     *
     * @param record SAM record to retrieve UMI from.
     * @param umiTag The tag used in the bam file that designates the UMI, null returns null
     * @return Normalized Duplex UMI.  If the UMI isn't duplex, it returns the UMI unaltered.
     */
    pub(crate) fn get_top_strand_normalized_umi<'r>(
        rec: &'r Record,
        umi_tag: &str,
        duplex_umi: bool,
    ) -> Result<String, Error> {
        if umi_tag.is_empty() {
            return Ok("".to_string());
        }

        let umi = match rec.get_str_aux(umi_tag.as_bytes()) {
            Ok(s) => s,
            Err(err) => Err(err)?,
        };

        if Self::allowed_umi().captures(umi).is_none() {
            Err(anyhow!(
                "UMI found with illegal characters. \
            UMIs must match the regular expression ^[ATCGNatcgn-]*$."
            ))?
        }

        if !duplex_umi {
            return Ok(umi.to_string());
        }

        let split = umi.split(Self::DUPLEX_UMI_DELIMITER).collect::<Vec<_>>();
        if split.len() != 2 {
            Err(anyhow!(
                "Duplex UMIs must be of the form X-Y \
            where X and Y are equal length UMIs, for example AT-GA.  Found UMI, {}",
                umi
            ))?
        }

        match Self::get_strand(rec)? {
            ReadStrand::BOTTOM => Ok(format!(
                "{}{}{}",
                split.get(1).unwrap(),
                Self::DUPLEX_UMI_DELIMITER,
                split.get(0).unwrap()
            )),
            _ => Ok(umi.to_string()),
        }
    }

    /**
     * Determines if the read represented by a SAM record belongs to the top or bottom strand
     * or if it cannot determine strand position due to one of the reads being unmapped.
     * Top strand is defined as having a read 1 unclipped 5' coordinate
     * less than the read 2 unclipped 5' coordinate.  If a read is unmapped
     * we do not attempt to determine the strand to which the read or its mate belongs.
     * If the mate belongs to a different contig from the read, then the reference
     * index for the read and its mate is used in leu of the unclipped 5' coordinate.
     * @param rec Record to determine top or bottom strand
     * @return Top or bottom strand, unknown if it cannot be determined.
     */
    fn get_strand(rec: &Record) -> Result<ReadStrand, Error> {
        if rec.is_unmapped() || rec.is_mate_unmapped() {
            return Ok(ReadStrand::UNKNOWN);
        }

        // If the read pair are aligned to different contigs we use
        // the reference index to determine relative 5' coordinate ordering.
        // Both the read and its mate should not have their unmapped flag set to true.
        if rec.tid() != rec.mtid() {
            if rec.is_first_in_template() == (rec.tid() < rec.mtid()) {
                return Ok(ReadStrand::TOP);
            } else {
                return Ok(ReadStrand::BOTTOM);
            }
        }

        let read_5prime_start = if rec.is_reverse() {
            rec.get_unclipped_end()
        } else {
            rec.get_unclipped_start()
        };

        let mate_5prime_start = if rec.is_mate_reverse() {
            rec.get_mate_unclipped_end()
        } else {
            rec.get_mate_unclipped_start()
        }?;

        if rec.is_first_in_template() == (read_5prime_start <= mate_5prime_start) {
            Ok(ReadStrand::TOP)
        } else {
            Ok(ReadStrand::BOTTOM)
        }
    }

    /**
     * Set molecular identifier tag of record using UMI and alignment start position along with top and bottom strand information
     * @param rec SAMRecord to set molecular index of
     * @param assignedUmi Assigned or inferred UMI to use in the molecular index tag
     * @param molecularIdentifierTag SAM tag to use as molecular identifier, if null no modification to the record will be made
     * @param duplexUmis Treat UMI as duplex, if true /A and /B will be added to denote top and bottom strands respectively
     */
    pub(crate) fn set_molecular_identifier(
        rec: &mut Record,
        assigned_umi: &str,
        molecular_identifier_tag: &str,
        duplex_umis: bool,
    ) {
        if molecular_identifier_tag.is_empty() {
            return;
        }

        let pos_from_strand = if rec.is_reverse() {
            rec.reference_start()
        } else {
            rec.reference_end()
        };

        let mut molecular_identifier = format!(
            "{}{}{}{}{}",
            rec.contig(),
            Self::CONTIG_SEPARATOR,
            pos_from_strand,
            Self::UMI_NAME_SEPARATOR,
            assigned_umi
        );

        if duplex_umis {
            // Reads whose strand position can be determined will have their
            // molecularIdentifier set to include an identifier appended that
            // indicates top or bottom strand.
            let strand = Self::get_strand(rec).unwrap();
            match strand {
                ReadStrand::TOP => molecular_identifier.push_str(Self::TOP_STRAND_DUPLEX),
                ReadStrand::BOTTOM => molecular_identifier.push_str(Self::BOTTOM_STRAND_DUPLEX),
                ReadStrand::UNKNOWN => {
                    // If we can't determine strand position nothing
                    // is appended to the molecularIdentifier.
                }
            }
        }

        rec.push_aux(
            molecular_identifier_tag.as_bytes(),
            Aux::String(&molecular_identifier),
        )
        .unwrap();
    }
}

/**
 * An enum to hold the strand position (TOP or BOTTOM) of a read.
 */
pub(crate) enum ReadStrand {
    TOP,
    BOTTOM,
    UNKNOWN,
}
