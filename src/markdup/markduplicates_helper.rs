use macro_sup::set_mlog;
use rust_htslib::bam::{
    ext::BamRecordExtensions,
    record::{Aux, RecordExt},
    HeaderView, Record,
};

use crate::{
    hts::{
        duplicate_scoring_strategy::DuplicateScoringStrategy, utils::sorting_collection::CowForSC,
    },
    markdup::markduplicates::MarkDuplicatesExt,
    utils::hash_code,
};

use super::{
    markduplicates::MarkDuplicates,
    umi_utils::UmiUtil,
    utils::{
        read_ends::ReadEnds, read_ends_for_mark_duplicates::ReadEndsForMarkDuplicates,
        read_ends_for_mark_duplicates_with_barcodes::ReadEndsBarcodeData,
        read_name_parser::ReadNameParserExt,
    },
};
use anyhow::{anyhow, Error};

set_mlog!(stringify!(MarkDuplicatesHelper));

pub(crate) enum MarkDuplicatesHelper {
    MarkDuplicatesHelper,
    MarkDuplicatesForFlowHelper,
}

impl MarkDuplicatesHelper {
    const CLIPPING_TAG_NAME: &'static str = "tm";
    const CLIPPING_TAG_CONTAINS_A: [char; 1] = ['A'];
    const CLIPPING_TAG_CONTAINS_AQ: [char; 2] = ['A', 'Q'];
    const CLIPPING_TAG_CONTAINS_QZ: [char; 2] = ['Q', 'Z'];

    const END_INSIGNIFICANT_VALUE: i32 = 0;

    const ATTR_DUPLICATE_SCORE: &'static str = "ForFlowDuplicateScore";

    pub(super) fn validate_flow_parameters(md: &MarkDuplicates) {
        if md.cli.UNPAIRED_END_UNCERTAINTY != 0 && !md.cli.USE_END_IN_UNPAIRED_READS {
            panic!("Invalid parameter combination. UNPAIRED_END_UNCERTAINTY can not be specified when USE_END_IN_UNPAIRED_READS not specified");
        }
    }

    /**
     * Builds a read ends object that represents a single read.
     */
    pub(super) fn build_read_ends_normal(
        md: &mut MarkDuplicates,
        header: &HeaderView,
        index: i64,
        rec: &Record,
        use_barcode: bool,
    ) -> Result<ReadEndsForMarkDuplicates, Error> {
        let mut read_ends = ReadEnds::default();

        read_ends.read1_reference_index = rec.tid();
        read_ends.read1_coordinate = if rec.is_reverse() {
            rec.get_unclipped_end() as i32
        } else {
            rec.get_unclipped_start() as i32 + 1 // +1 to convert it into 1-based.
        };

        read_ends.orientation = if rec.is_reverse() {
            ReadEnds::R
        } else {
            ReadEnds::F
        };

        let read1_index_in_file = index;
        let score = DuplicateScoringStrategy::compute_duplicate_score(
            rec,
            md.cli.DUPLICATE_SCORING_STRATEGY,
            false,
        )
        .unwrap(); // we can assure that no error occurs when `assume_mate_cigar=false`.

        // Doing this lets the ends object know that it's part of a pair
        if rec.is_paired() && !rec.is_mate_unmapped() {
            read_ends.read2_reference_index = rec.mtid()
        }

        // Fill in the library ID
        read_ends.library_id = md.library_id_generator.get_library_id(rec)?;

        // Fill in the location information for optical duplicates
        if md.optical_duplicate_finder.add_location_information(
            std::str::from_utf8(rec.qname())?.to_string(),
            &mut read_ends,
        ) {
            // calculate the RG number (nth in list)
            read_ends.read_group = 0;

            let rg = rec
                .aux(b"RG")
                .or_else(|err| Err(Error::from(err)))
                .and_then(|aux| match aux {
                    Aux::String(v) => Ok(v),
                    _ => Err(anyhow!("Invalid type for RG tag")),
                });
            let read_groups = header.get_read_groups();

            if let (Ok(rg), Ok(read_groups)) = (rg, read_groups) {
                for read_group in read_groups {
                    if read_group.get_read_group_id().eq(rg) {
                        break;
                    } else {
                        read_ends.read_group += 1;
                    }
                }
            }
        }

        let mut read_ends_for_md = ReadEndsForMarkDuplicates {
            score,
            read1_index_in_file,
            read_ends,
            ..Default::default()
        };

        if use_barcode {
            let top_strand_normalized_umi = UmiUtil::get_top_strand_normalized_umi(
                rec,
                &md.cli.BARCODE_TAG,
                md.cli.DUPLEX_UMI,
            )?;

            let barcode = hash_code(top_strand_normalized_umi);

            let mut barcode_data = ReadEndsBarcodeData {
                barcode,
                ..Default::default()
            };

            if (!rec.is_paired() || rec.is_first_in_template()) {
                barcode_data.read_one_barcode = md.get_read_one_barcode_value(rec);
            } else {
                barcode_data.read_two_barcode = md.get_read_two_barcode_value(rec);
            }

            read_ends_for_md.barcode_data = Some(barcode_data);
        }

        read_ends_for_md.read1_qname = std::str::from_utf8(rec.qname()).unwrap().to_string();

        Ok(read_ends_for_md)
    }

    /**
     * Builds a read ends object that represents a single read - for flow based read
     */
    pub(super) fn build_read_ends_for_flow(
        md: &mut MarkDuplicates,
        header: &HeaderView,
        index: i64,
        rec: &Record,
        use_barcode: bool,
    ) -> Result<ReadEndsForMarkDuplicates, Error> {
        let mut ends = Self::build_read_ends_normal(md, header, index, rec, use_barcode)?;

        // this code only supported unpaired reads
        if rec.is_paired() && !rec.is_mate_unmapped() {
            panic!(
                "FLOW_MODE does not support paired reads. offending read: {}",
                std::str::from_utf8(rec.qname())?
            );
        }

        // adjust start/end coordinates
        ends.read_ends.read1_coordinate =
            Self::get_read_end_coordinate(rec, !rec.is_reverse(), true, &md.cli);

        if md.cli.USE_UNPAIRED_CLIPPED_END {
            ends.read_ends.read2_coordinate =
                Self::get_read_end_coordinate(rec, rec.is_reverse(), false, &md.cli);
        }

        // adjust score
        if md.cli.FLOW_QUALITY_SUM_STRATEGY {
            ends.score =
                Self::compute_flow_duplicate_score(md, rec, ends.read_ends.read2_coordinate);
        }

        Ok(ends)
    }

    fn get_read_end_coordinate(
        rec: &Record,
        start_end: bool,
        certain: bool,
        cli: &crate::cmdline::cli::Cli,
    ) -> i32 {
        let mut flow_order = FlowOrder::new(rec);
        let unclipped_coor = if start_end {
            rec.get_unclipped_start()
        } else {
            rec.get_unclipped_end()
        } as i32;

        let alignment_coor = if start_end {
            rec.reference_start()
        } else {
            rec.reference_end()
        } as i32;

        // this code requires a valid flow order
        if flow_order.is_valid() {
            // simple case
            if cli.USE_UNPAIRED_CLIPPED_END {
                return alignment_coor;
            }

            // "skipping" case
            if certain && cli.FLOW_SKIP_FIRST_N_FLOWS != 0 {
                let bases = rec.seq().encoded;
                let mut hmer_base = if start_end {
                    bases[0]
                } else {
                    bases[bases.len() - 1]
                };

                let mut hmers_left = cli.FLOW_SKIP_FIRST_N_FLOWS; // number of hmer left to trim

                // advance flow order to base
                while flow_order.current() != hmer_base {
                    flow_order.advance();
                    hmers_left -= 1;
                }

                let mut hmer_size = 0;
                while {
                    hmer_size += 1;
                    hmer_size < bases.len()
                } {
                    let query = if start_end {
                        bases[hmer_size]
                    } else {
                        bases[bases.len() - 1 - hmer_size]
                    };

                    if query != hmer_base {
                        hmers_left -= 1;
                        if hmers_left <= 0 {
                            break;
                        } else {
                            hmer_base = if start_end {
                                bases[hmer_size]
                            } else {
                                bases[bases.len() - 1 - hmer_size]
                            };

                            flow_order.advance();

                            while flow_order.current() != hmer_base {
                                flow_order.advance();
                                hmers_left -= 1;
                            }

                            if hmers_left <= 0 {
                                break;
                            }
                        }
                    }
                }

                let coor = unclipped_coor
                    + if start_end {
                        hmer_size as i32
                    } else {
                        -(hmer_size as i32)
                    };

                return if cli.USE_UNPAIRED_CLIPPED_END {
                    if start_end {
                        coor.max(alignment_coor)
                    } else {
                        coor.min(alignment_coor)
                    }
                } else {
                    coor
                };
            }

            // "known end" case
            if {
                if cli.FLOW_Q_IS_KNOWN_END {
                    Self::is_adapter_clipped(rec)
                } else {
                    Self::is_adapter_clipped_with_q(rec)
                }
            } {
                return unclipped_coor;
            }

            // "uncertain quality clipped" case
            if !certain && Self::is_quality_clipped(rec) {
                return Self::END_INSIGNIFICANT_VALUE;
            }
        }

        // if here, return a default
        unclipped_coor
    }

    fn clipping_tag_contains_any(rec: &Record, chars: &[char]) -> bool {
        if let Ok(clipping_tag_value) = rec.get_str_aux(Self::CLIPPING_TAG_NAME.as_bytes()) {
            for ch in chars {
                if clipping_tag_value.find(*ch).is_some() {
                    return true;
                }
            }
            return false;
        } else {
            return false;
        };
    }

    fn is_adapter_clipped(rec: &Record) -> bool {
        Self::clipping_tag_contains_any(rec, &Self::CLIPPING_TAG_CONTAINS_A)
    }

    fn is_adapter_clipped_with_q(rec: &Record) -> bool {
        Self::clipping_tag_contains_any(rec, &Self::CLIPPING_TAG_CONTAINS_AQ)
    }

    fn is_quality_clipped(rec: &Record) -> bool {
        Self::clipping_tag_contains_any(rec, &Self::CLIPPING_TAG_CONTAINS_QZ)
    }

    fn compute_flow_duplicate_score(md: &MarkDuplicates, rec: &Record, end: i32) -> i16 {
        if end == Self::END_INSIGNIFICANT_VALUE {
            return -1;
        }

        // cached attribute not implemented yet.
        // java code search there is cached attribute first if then get it.
        // To do it in Rust, each record object must not removed or replaced on every read() function, but it is.
        // For now, ignore the part loading cached attribute. Just make it every time this method is called.

        let mut score = 0_i16;
        score +=
            (Self::get_flow_sum_of_base_qualities(rec, md.cli.FLOW_EFFECTIVE_QUALITY_THRESHOLD)
                as i16)
                .min(i16::MAX / 2);

        score += if rec.is_quality_check_failed() {
            i16::MIN / 2
        } else {
            0
        };

        score
    }

    fn get_flow_sum_of_base_qualities(rec: &Record, threshold: i32) -> i32 {
        let mut score = 0;

        // access qualities and bases
        let quals = rec.qual();
        let bases = rec.seq().encoded;

        // create iteration range and direction
        let starting_offset = if !rec.is_reverse() { 0 } else { bases.len() };

        let end_offset = if !rec.is_reverse() { bases.len() } else { 0 };

        let iter_incr = if !rec.is_reverse() {
            |x: &mut usize| *x += 1
        } else {
            |x: &mut usize| *x -= 1
        };

        // loop on bases, extract qual related to homopolymer from start of homopolymer
        let mut last_base = 0;
        let mut effective_qual = 0;

        let mut i = starting_offset;

        let threshold = threshold as u8;
        while {
            iter_incr(&mut i);
            i != end_offset
        } {
            let base = bases[i];
            if base != last_base {
                effective_qual = quals[i];
            }

            if effective_qual >= threshold {
                score += effective_qual;
            }

            last_base = base;
        }

        score as i32
    }

    pub(crate) fn get_read_duplicate_score(
        &self,
        md: &MarkDuplicates,
        rec: &Record,
        paired_ends: &ReadEndsForMarkDuplicates,
    ) -> i16 {
        match self {
            MarkDuplicatesHelper::MarkDuplicatesHelper => {
                DuplicateScoringStrategy::compute_duplicate_score(
                    rec,
                    md.cli.DUPLICATE_SCORING_STRATEGY,
                    false,
                )
                .unwrap()
            }
            MarkDuplicatesHelper::MarkDuplicatesForFlowHelper => {
                if md.cli.FLOW_QUALITY_SUM_STRATEGY {
                    Self::compute_flow_duplicate_score(
                        md,
                        rec,
                        paired_ends.read_ends.read2_coordinate,
                    )
                } else {
                    md.get_read_duplicate_score(rec, paired_ends)
                }
            }
        }
    }

    fn generate_duplicate_indexes(
        &self,
        md: &mut MarkDuplicates,
        use_barcodes: bool,
        index_optical_duplicates: bool,
    ) {
        md.sort_indices_for_duplicates(index_optical_duplicates);

        let mut next_chunk: Vec<CowForSC<ReadEndsForMarkDuplicates>> = Vec::with_capacity(200);

        // First just do the pairs
        mlog::info!("Traversing read pair information and detecting duplicates.");

        let mut first_of_next_chunk;

        let mut pair_sort_iter = md.pair_sort.iter().unwrap();

        if let Some(v) = pair_sort_iter.next() {
            // first_of_next_chunk = v;

            next_chunk.push(v);
            first_of_next_chunk = next_chunk.last().unwrap();

            // pair_sort_iter.for_each(|next| {
            //     if MarkDuplicates::are_comparable_for_duplicates(
            //         &first_of_next_chunk,
            //         &next,
            //         true,
            //         use_barcodes,
            //     ) {
            //         next_chunk.push(next);
            //     } else {
            //         md.handle_chunk(next_chunk);
            //         next_chunk.clear();
            //         next_chunk.push(next);
            //         first_of_next_chunk = next_chunk.last().unwrap();
            //     }
            // })
        }
    }
}

struct FlowOrder {
    flow_order: String, // the flow order string
    flow_index: i32,    // the current position on the flow order
}

impl FlowOrder {
    fn new(rec: &Record) -> Self {
        let mut s = Self {
            flow_order: String::new(),
            flow_index: 0,
        };

        // access flow order from record's read group
        if let Ok(Some(fo)) = rec.get_read_group().and_then(|rg| Ok(rg.get_flow_order())) {
            s.flow_order.push_str(fo);
            return s;
        }

        // fallback on finding a flow order elsewhere
        let header = rec.header().unwrap();
        let header_map = header.header_map();
        if let Some(rgs) = header_map.get_read_groups() {
            for rg in rgs {
                if let Some(fo) = rg.get_flow_order() {
                    s.flow_order = fo.to_string();
                    return s;
                }
            }
        }

        s
    }

    fn is_valid(&self) -> bool {
        !self.flow_order.is_empty()
    }

    fn current(&self) -> u8 {
        self.flow_order.as_bytes()[self.flow_index as usize]
    }

    fn advance(&mut self) {
        self.flow_index += 1;
        if self.flow_index >= self.flow_order.len() as i32 {
            self.flow_index = 0;
        }
    }
}
