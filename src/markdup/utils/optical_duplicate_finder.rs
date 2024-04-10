use std::{collections::HashMap, fmt};

use macro_sup::set_mlog;
use serde::{Deserialize, Serialize};

use crate::utils::{graph_utils::Graph, logging::ProgressLogger};

use super::{
    physical_location::PhysicalLocation,
    read_name_parser::{impl_read_name_parser_ext, ReadNameParser, ReadNameParserExt},
};
use regex::Regex;

type Error = Box<dyn std::error::Error + Send + Sync>;

#[derive(Serialize, Deserialize)]
pub(crate) struct OpticalDuplicateFinder {
    pub(crate) optical_duplicate_pixel_distance: i32,
    big_duplicate_set_size: i32,
    max_duplicate_set_size: i64,

    rnp: ReadNameParser,
}

impl Default for OpticalDuplicateFinder {
    fn default() -> Self {
        Self {
            optical_duplicate_pixel_distance: Self::DEFAULT_OPTICAL_DUPLICATE_DISTANCE,
            big_duplicate_set_size: Self::DEFAULT_BIG_DUPLICATE_SET_SIZE,
            max_duplicate_set_size: Self::DEFAULT_MAX_DUPLICATE_SET_SIZE,
            rnp: ReadNameParser::new(),
        }
    }
}

set_mlog!(stringify!(OpticalDuplicateFinder));

impl OpticalDuplicateFinder {
    pub(crate) const DEFAULT_OPTICAL_DUPLICATE_DISTANCE: i32 = 100;
    pub(crate) const DEFAULT_MAX_DUPLICATE_SET_SIZE: i64 = 300000;
    pub(crate) const DEFAULT_BIG_DUPLICATE_SET_SIZE: i32 = 1000;

    #[allow(non_upper_case_globals)]
    const log: &'static str = stringify!(OpticalDuplicateFinder);

    pub(crate) fn find_optical_duplicates<P>(&self, loc_vec: &[&P], keeper: Option<&P>) -> Vec<bool>
    where
        P: PhysicalLocation + PartialEq + fmt::Debug,
    {
        let length = loc_vec.len() as i32;

        let optical_duplicate_flags = vec![false; length as usize];

        // If there is only one or zero reads passed in (so there are obviously no optical duplicates),
        // or if there are too many reads (so we don't want to try to run this expensive n^2 algorithm),
        // then just return an array of all false
        if self.rnp.read_name_regex.is_empty()
            || length < 2
            || length as i64 > self.max_duplicate_set_size
        {
            return optical_duplicate_flags;
        }

        let actual_keeper = Self::keeper_or_null(loc_vec, keeper);

        let log_progress = length > self.big_duplicate_set_size;

        let (progress_logger_for_keeper, progress_logger_for_rest);
        if log_progress {
            progress_logger_for_keeper = Some(ProgressLogger::new(
                Self::log,
                10000,
                "compared",
                "ReadEnds to keeper",
            ));
            progress_logger_for_rest = Some(ProgressLogger::new(
                Self::log,
                1000,
                "compared",
                "ReadEnds to others",
            ));

            mlog::info!("Large duplicate set. size = {}", length);
            mlog::debug!("About to compare to keeper: {:?}", actual_keeper);
        } else {
            progress_logger_for_keeper = None;
            progress_logger_for_rest = None;
        }

        if length >= (if keeper.is_none() { 3 } else { 4 }) {
            self.get_optical_duplicates_flag_with_graph(
                loc_vec,
                actual_keeper,
                optical_duplicate_flags,
                progress_logger_for_keeper,
                progress_logger_for_rest,
                log_progress,
            )
        } else {
            self.get_optical_duplicates_flag_fast(
                loc_vec,
                actual_keeper,
                optical_duplicate_flags,
                progress_logger_for_keeper,
                progress_logger_for_rest,
                log_progress,
            )
        }
    }
    /** Returns the keeper if it is contained within the list and has location information, otherwise null. */
    fn keeper_or_null<'k, P: PhysicalLocation + PartialEq>(
        loc_vec: &[&P],
        keeper: Option<&'k P>,
    ) -> Option<&'k P> {
        if let Some((k, true)) = keeper.and_then(|v| Some((v, v.has_location()))) {
            for loc in loc_vec {
                if (*loc).eq(k) {
                    return Some(k);
                }
            }
        }

        None
    }

    fn get_optical_duplicates_flag_with_graph<P>(
        &self,
        loc_vec: &[&P],
        keeper: Option<&P>,
        mut optical_duplicate_flags: Vec<bool>,
        mut progress_logger_for_keeper: Option<ProgressLogger>,
        mut progress_logger_for_rest: Option<ProgressLogger>,
        log_progress: bool,
    ) -> Vec<bool>
    where
        P: PhysicalLocation + PartialEq + fmt::Debug,
    {
        // Make a graph where the edges are reads that lie within the optical duplicate pixel distance from each other,
        // we will then use the union-find algorithm to cluster the graph and find optical duplicate groups
        let mut optical_distance_relation_graph: Graph<i32> = Graph::new();
        if log_progress {
            mlog::debug!("Building adjacency graph for duplicate group");
        }

        let mut tile_r_map: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut keeper_index = -1_i32;

        for (i, currect_loc) in (0..loc_vec.len() as i32).zip(loc_vec.iter().copied()) {
            if let Some(true) = keeper.and_then(|ac| Some(currect_loc == ac)) {
                keeper_index = i;
            }

            if currect_loc.has_location() {
                let key = currect_loc.get_read_group().wrapping_shl(16) as i32
                    + currect_loc.get_tile() as i32;

                match tile_r_map.entry(key) {
                    std::collections::hash_map::Entry::Occupied(mut ent) => {
                        ent.get_mut().push(key);
                    }
                    std::collections::hash_map::Entry::Vacant(ent) => {
                        ent.insert(vec![i]);
                    }
                }
            }

            optical_distance_relation_graph.add_node(i);
        }

        // Since because finding adjacent optical duplicates is an O(n^2) operation, we can subdivide the input into its
        // readgroups in order to reduce the amount of redundant checks across readgroups between reads.
        for tile_group in tile_r_map.values() {
            if tile_group.len() > 1 {
                Self::fill_graph_from_a_group(
                    loc_vec,
                    tile_group,
                    log_progress,
                    progress_logger_for_keeper.as_mut(),
                    self.optical_duplicate_pixel_distance,
                    &mut optical_distance_relation_graph,
                );
            }
        }

        if log_progress {
            mlog::debug!(
                "Finished building adjacency graph for duplicate group, moving onto clustering"
            );
        }

        // Keep a map of the reads and their cluster assignments
        let optical_duplicate_cluster_map = optical_distance_relation_graph.cluster();
        let mut cluster_to_representative_read: HashMap<i32, i32> = HashMap::new();

        let mut keeper_cluster: Option<i32> = None;

        // Specially mark the keeper as specifically not a duplicate if it exists
        if keeper_index >= 0 {
            let cluster_map_value = optical_duplicate_cluster_map
                .get(&keeper_index)
                .unwrap()
                .clone();
            cluster_to_representative_read.insert(cluster_map_value, keeper_index);
            keeper_cluster = Some(cluster_map_value);
        }

        for (record_index, record_assigned_cluster) in
            optical_duplicate_cluster_map.iter().map(|(k, v)| (**k, *v))
        {
            // logging here for same reason as above
            if log_progress {
                let loc = loc_vec.get(record_index as usize).unwrap();
                progress_logger_for_rest
                    .as_mut()
                    .unwrap()
                    .record_with_chrom_pos(&loc.get_read_group().to_string(), loc.get_x() as i64);
            }

            // If its not the first read we've seen for this cluster, mark it as an optical duplicate
            if let Some((rrv, true)) = cluster_to_representative_read
                .get(&record_assigned_cluster)
                .and_then(|v| Some((*v, record_index.ne(&keeper_index))))
            {
                let representative_loc = loc_vec.get(rrv as usize).unwrap();
                let current_record_loc = loc_vec.get(record_index as usize).unwrap();

                // If not in the keeper cluster, then keep the minX -> minY valued duplicate (note the tile must be equal for reads to cluster together)
                if !(keeper_index >= 0
                    && keeper_cluster
                        .as_ref()
                        .map_or_else(|| false, |v| record_assigned_cluster.eq(v)))
                    && (current_record_loc.get_x() < representative_loc.get_x()
                        || (current_record_loc.get_x() == representative_loc.get_x()
                            && current_record_loc.get_y() < representative_loc.get_y()))
                {
                    // Mark the old min as an optical duplicate, and save the new min
                    *optical_duplicate_flags.get_mut(rrv as usize).unwrap() = true;
                    cluster_to_representative_read.insert(record_assigned_cluster, record_index);
                } else {
                    // If a smaller read has already been visited, mark the test read as an optical duplicate
                    *optical_duplicate_flags
                        .get_mut(record_index as usize)
                        .unwrap() = true;
                }
            } else {
                cluster_to_representative_read.insert(record_assigned_cluster, record_index);
            }
        }

        optical_duplicate_flags
    }

    fn fill_graph_from_a_group<P>(
        whole_vec: &[&P],
        group_vec: &[i32],
        log_progress: bool,
        mut progress_logger_for_keeper: Option<&mut ProgressLogger>,
        distance: i32,
        optical_distance_relation_graph: &mut Graph<i32>,
    ) where
        P: PhysicalLocation + PartialEq + fmt::Debug,
    {
        for i in (0..group_vec.len()) {
            let i_index = group_vec.get(i).unwrap().clone();

            let current_loc = whole_vec.get(i_index as usize).unwrap().clone();

            // The main point of adding this log and if statement (also below) is a workaround a bug in the JVM
            // which causes a deep exception (https://github.com/broadinstitute/picard/issues/472).
            // It seems that this is related to https://bugs.openjdk.java.net/browse/JDK-8033717 which
            // was closed due to non-reproducibility. We came across a bam file that evoked this error
            // every time we tried to duplicate-mark it. The problem seemed to be a duplicate-set of size 500,000,
            // and this loop seemed to kill the JVM for some reason. This logging statement (and the one in the
            // loop below) solved the problem.
            if log_progress {
                progress_logger_for_keeper
                    .as_mut()
                    .unwrap()
                    .record_with_chrom_pos(
                        &current_loc.get_read_group().to_string(),
                        current_loc.get_x() as i64,
                    );
            }

            for j in ((i + 1)..group_vec.len()) {
                let j_index = group_vec.get(j).unwrap().clone();
                let other = whole_vec.get(j_index as usize).unwrap();

                if Self::close_enough_short(current_loc, other, distance) {
                    optical_distance_relation_graph.add_edge(i_index, j_index);
                }
            }
        }
    }

    fn close_enough_short<L>(lhs: &L, rhs: &L, distance: i32) -> bool
    where
        L: PhysicalLocation + PartialEq + fmt::Debug,
    {
        lhs != rhs
            && (lhs.get_x() - rhs.get_x()).abs() <= distance
            && (lhs.get_y() - rhs.get_y()).abs() <= distance
    }

    fn get_optical_duplicates_flag_fast<P>(
        &self,
        loc_vec: &[&P],
        actual_keeper: Option<&P>,
        mut optical_duplicate_flags: Vec<bool>,
        progress_logger_for_keeper: Option<ProgressLogger>,
        mut progress_logger_for_rest: Option<ProgressLogger>,
        log_progress: bool,
    ) -> Vec<bool>
    where
        P: PhysicalLocation + PartialEq + fmt::Debug,
    {
        let length = loc_vec.len();

        // First go through and compare all the reads to the keeper
        if let Some(ac) = actual_keeper {
            for (other, opt_dup_flag) in loc_vec.iter().zip(optical_duplicate_flags.iter_mut()) {
                *opt_dup_flag =
                    Self::close_enough(ac, other, self.optical_duplicate_pixel_distance);
                // The main point of adding this log and if statement (also below) is a workaround a bug in the JVM
                // which causes a deep exception (https://github.com/broadinstitute/picard/issues/472).
                // It seems that this is related to https://bugs.openjdk.java.net/browse/JDK-8033717 which
                // was closed due to non-reproducibility. We came across a bam file that evoked this error
                // every time we tried to duplicate-mark it. The problem seemed to be a duplicate-set of size 500,000,
                // and this loop seemed to kill the JVM for some reason. This logging statement (and the one in the
                // loop below) solved the problem.
            }
        }

        if log_progress {
            mlog::debug!("Done with comparing to keeper, now the rest.");
        }

        // Now go through and do each pairwise comparison not involving the actualKeeper
        for (i, lhs) in (0..(length)).zip(loc_vec.iter().copied()) {
            // no comparisons to actualKeeper since those are all handled above
            if let Some(true) = actual_keeper.and_then(|ac| Some(lhs == ac)) {
                continue;
            }

            // logging here for same reason as above
            if log_progress {
                progress_logger_for_rest
                    .as_mut()
                    .unwrap() // can assure it is Some when `log_progress` is true
                    .record_with_chrom_pos(&lhs.get_read_group().to_string(), lhs.get_x() as i64);
            }

            for j in ((i + 1)..length) {
                let rhs = loc_vec.get(j as usize).unwrap().clone();
                if let Some(true) = actual_keeper.and_then(|ac| Some(rhs == ac)) {
                    continue; // no comparisons to actualKeeper since those are all handled above
                }

                if optical_duplicate_flags[i] && optical_duplicate_flags[j] {
                    continue; // both already marked, no need to check
                }

                if Self::close_enough(lhs, rhs, self.optical_duplicate_pixel_distance) {
                    // At this point we want to mark either lhs or rhs as duplicate. Either could have been marked
                    // as a duplicate of the keeper (but not both - that's checked above), so be careful about which
                    // one to now mark as a duplicate.

                    let index = if optical_duplicate_flags[j] { i } else { j };
                    optical_duplicate_flags[index] = true;
                }
            }
        }

        optical_duplicate_flags
    }

    /** Simple method to test whether two physical locations are close enough to each other to be deemed optical dupes. */
    fn close_enough<P>(lhs: &P, rhs: &P, distance: i32) -> bool
    where
        P: PhysicalLocation + PartialEq + fmt::Debug,
    {
        lhs != rhs // no comparing an object to itself (checked using object identity)!
            && lhs.has_location() // no comparing objects without locations
            && rhs.has_location()
            && lhs.get_read_group() == rhs.get_read_group() // must be in the same RG to be optical duplicates
            && lhs.get_tile() == rhs.get_tile() // and the same tile
            && (lhs.get_x() - rhs.get_x()).abs() <= distance
            && (lhs.get_y() - rhs.get_y()).abs() <= distance
    }
}

impl ReadNameParserExt for OpticalDuplicateFinder {
    impl_read_name_parser_ext!(self, self.rnp);
}
