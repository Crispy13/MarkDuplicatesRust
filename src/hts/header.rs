use rust_htslib::bam::Header;
use std::fmt::Write;
use std::{collections::HashSet, sync::OnceLock};

pub(crate) struct PgIdGenerator {
    record_counter: usize,
    ids_that_are_already_taken: HashSet<String>,
}

impl PgIdGenerator {
    fn new(header: &Header) -> Self {
        let mut ids_that_are_already_taken = HashSet::new();
        header
            .to_hashmap()
            .remove("PG")
            .unwrap()
            .into_iter()
            .for_each(|mut pgr| {
                ids_that_are_already_taken.insert(pgr.remove("ID").unwrap());
            });

        let record_counter = ids_that_are_already_taken.len();

        Self {
            ids_that_are_already_taken,
            record_counter,
        }
    }

    fn get_non_colliding_id(&mut self, record_id: String) -> String {
        if !self.ids_that_are_already_taken.contains(&record_id) {
            // don't remap 1st record. If there are more records
            // with this id, they will be remapped in the 'else'.
            self.ids_that_are_already_taken.insert(record_id.clone());
            self.record_counter += 1;

            record_id
        } else {
            let mut new_id = String::new();
            // Below we tack on one of roughly 1.7 million possible 4 digit base36 at random. We do this because
            // our old process of just counting from 0 upward and adding that to the previous id led to 1000s of
            // calls idsThatAreAlreadyTaken.contains() just to resolve 1 collision when merging 1000s of similarly
            // processed bams.
            while {
                write!(
                    &mut new_id,
                    "{}.{}",
                    record_id,
                    SamFileHeaderMerger::positive_four_digit_base36_str(self.record_counter)
                )
                .unwrap();
                self.record_counter += 1;

                self.ids_that_are_already_taken.contains(&new_id)
            } {}

            new_id
        }
    }
}

struct SamFileHeaderMerger {
    // int_to_base36: [char; 36],
}

impl SamFileHeaderMerger {
    fn new() -> Self {
        Self {}
    }

    pub(crate) fn int_to_base36() -> &'static [char; 36] {
        #[allow(non_upper_case_globals)]
        static int_to_base36_once_lock: OnceLock<[char; 36]> = OnceLock::new();
        int_to_base36_once_lock.get_or_init(|| {
            let mut int_to_base36 = [Default::default(); 36];
            let a_val = 'A' as u32;
            let zero_val = '0' as u32;

            for i in (0..10) {
                int_to_base36[i] = char::from_u32(zero_val + i as u32).unwrap();
            }

            for i in (0..26) {
                int_to_base36[i + 10] = char::from_u32(a_val + i as u32).unwrap();
            }

            int_to_base36
        })
    }

    fn positive_four_digit_base36_str(mut left_over: usize) -> String {
        if left_over == 0 {
            return "0".into();
        }

        let mut char_vec = Vec::with_capacity(10);

        let int_to_base_36 = Self::int_to_base36();

        while left_over > 0 {
            let value_index = left_over % 36;
            char_vec.push(int_to_base_36[value_index]);
            left_over /= 36;
        }

        char_vec.into_iter().rev().collect::<String>()
    }
}
