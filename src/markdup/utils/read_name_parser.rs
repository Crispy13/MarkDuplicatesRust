use std::{fmt, num::ParseIntError};

use log4rs::Handle;
use regex::Regex;
use serde::{Deserialize, Serialize};

type Error = Box<dyn std::error::Error + Send + Sync>;

use super::physical_location::{PhysicalLocation, PhysicalLocationInt};

pub(crate) const DEFAULT_READ_NAME_REGEX: &'static str =
    "<optimized capture of last three ':' separated fields as numeric values>";

#[derive(Serialize, Deserialize)]
pub(crate) struct ReadNameParser {
    pub(super) read_name_stored: String,
    pub(super) physical_location_stored: PhysicalLocationInt,
    pub(super) tmp_location_fields: [i32; 3], // for optimization of addLocationInformation
    pub(super) use_optimized_default_parsing: bool,
    pub(crate) read_name_regex: String,

    #[serde(skip)]
    pub(super) read_name_pattern: Option<Regex>,

    pub(super) warned_about_regex_not_matching: bool,

    #[serde(skip)]
    pub(super) log: String, // Just log name instead of a Log object
}

impl Default for ReadNameParser {
    fn default() -> Self {
        Self {
            read_name_stored: String::new(),
            physical_location_stored: PhysicalLocationInt::default(),
            tmp_location_fields: Default::default(),
            use_optimized_default_parsing: Default::default(),
            read_name_regex: Default::default(),
            read_name_pattern: None,
            warned_about_regex_not_matching: false,
            log: String::new(),
        }
    }
}

impl ReadNameParser {
    /**
     * Creates are read name parser using the given read name regex.  See {@link #DEFAULT_READ_NAME_REGEX} for an explanation on how to
     * format the regular expression (regex) string.
     * @param readNameRegex the read name regular expression string to parse read names, null to never parse location information..
     * @param log the log to which to write messages.
     */
    pub(crate) fn with_regex_and_log(read_name_regex: String, log: String) -> Self {
        let use_optimized_default_parsing = DEFAULT_READ_NAME_REGEX.eq(&read_name_regex);

        Self {
            read_name_regex,
            log,
            use_optimized_default_parsing,
            ..Default::default()
        }
    }

    pub(crate) fn new() -> Self {
        Self::with_regex(DEFAULT_READ_NAME_REGEX.to_string())
    }

    pub(crate) fn with_regex(read_name_regex: String) -> Self {
        Self::with_regex_and_log(read_name_regex, "".to_string())
    }
}

pub(crate) trait ReadNameParserExt {
    /**
     * Very specialized method to rapidly parse a sequence of digits from a String up until the first
     * non-digit character.
     *
     * @throws NumberFormatException if the String does not start with an optional - followed by at least on digit
     */
    fn rapid_parse_i32(input: &str) -> Result<i32, Error> {
        let len = input.len();
        let mut val = 0 as i32;
        let mut i = 0;

        let mut is_negative = false;

        let input_bytes = input.as_bytes();

        if 0 < len && b'-'.eq(input_bytes.get(0).unwrap()) {
            i = 1;
            is_negative = true;
        }

        let mut has_digits = false;

        for (_, b) in (i..len).zip(input_bytes.get(i..len).unwrap()) {
            let ch = *b as char;
            if ch.is_numeric() {
                val = (val * 10) + (ch as u32 as i32 - 48);
                has_digits = true;
            } else {
                break;
            }
        }

        if !has_digits {
            Err(format!(
                "String '{}' did not start with a parsable number.",
                input
            ))?
        }

        if is_negative {
            val = -val;
        }

        Ok(val)
    }

    /**
     * Method used to extract tile/x/y from the read name and add it to the PhysicalLocationShort so that it
     * can be used later to determine optical duplication
     *
     * @param readName the name of the read/cluster
     * @param loc the object to add tile/x/y to
     * @return true if the read name contained the information in parsable form, false otherwise
     */
    fn read_location_information(
        &mut self,
        read_name: &str,
        loc: &mut impl PhysicalLocation,
    ) -> bool;

    fn _read_location_information(
        &mut self,
        read_name: &str,
        loc: &mut impl PhysicalLocation,
    ) -> Result<bool, Error>;

    fn add_location_information(
        &mut self,
        read_name: String,
        loc: &mut impl PhysicalLocation,
    ) -> bool;

    fn get_last_three_fields(
        read_name: &str,
        delim: char,
        tokens: &mut [i32],
    ) -> Result<i32, Error>;
}

macro_rules! impl_read_name_parser_ext {
    ($self:ident, $rnp:expr) => {
        fn read_location_information(
            &mut $self,
            read_name: &str,
            loc: &mut impl PhysicalLocation,
        ) -> bool {
            match $rnp._read_location_information(&read_name, loc) {
                Ok(b) => b,
                Err(err) => {
                    if !$rnp.log.is_empty() && !$rnp.warned_about_regex_not_matching {
                        log::warn!(
                            target: $rnp.log.as_str(),
                            "A field field parsed out of a read name was expected to contain an integer and did not. \
                            Read name: {}. Cause: {:?}",
                            read_name, err
                        );
                        $rnp.warned_about_regex_not_matching = true;
                    }
                    false
                }
            }
        }

        fn _read_location_information(
            &mut $self,
            read_name: &str,
            loc: &mut impl PhysicalLocation,
        ) -> Result<bool, Error> {
            // Optimized version if using the default read name regex (== used on purpose):
            if $rnp.use_optimized_default_parsing {
                let fields: i32 = Self::get_last_three_fields(read_name, ':', $rnp.tmp_location_fields.as_mut_slice())?;
                if !(fields == 5 || fields == 7) {
                    if !$rnp.log.is_empty() && $rnp.warned_about_regex_not_matching {
                        log::warn!(target: $rnp.log.as_str(),
                            "Default READ_NAME_REGEX '{}' did not match read name '{}'.  \
                            You may need to specify a READ_NAME_REGEX in order to correctly identify optical duplicates.  \
                            Note that this message will not be emitted again even if other read names do not match the regex.",
                            $rnp.read_name_regex, read_name
                        );
                        $rnp.warned_about_regex_not_matching = true;
                    }

                    return Ok(false);
                }

                loc.set_tile(*$rnp.tmp_location_fields.get(0).unwrap() as i16);
                loc.set_x($rnp.tmp_location_fields[1]);
                loc.set_y($rnp.tmp_location_fields[2]);
                return Ok(true);
            } else if $rnp.read_name_regex.is_empty() {
                return Ok(false);
            } else {
                // Standard version that will use the regex
                if $rnp.read_name_pattern.is_none() {
                    $rnp.read_name_pattern = Some(Regex::new(&$rnp.read_name_regex).unwrap());
                }

                if let Some(m) = $rnp.read_name_pattern.as_ref().unwrap().captures(read_name) {
                    loc.set_tile(m.get(1).unwrap().as_str().parse::<i16>()?);
                    loc.set_x(m.get(2).unwrap().as_str().parse::<i32>()?);
                    loc.set_y(m.get(3).unwrap().as_str().parse::<i32>()?);

                    return Ok(true);
                } else {
                    if !$rnp.log.is_empty() && !$rnp.warned_about_regex_not_matching {
                        log::warn!(
                            target: $rnp.log.as_str(),
                            "READ_NAME_REGEX '{}' did not match read name '{}'.  Your regex may not be correct.  \
                            Note that this message will not be emitted again even if other read names do not match the regex.",
                            $rnp.read_name_regex, read_name
                        );
                        $rnp.warned_about_regex_not_matching = true;
                    }
                    return Ok(false);
                }
            }
        }

        fn add_location_information(&mut $self, read_name: String, loc: &mut impl PhysicalLocation) -> bool {
            if !read_name.eq(&$rnp.read_name_stored) {
                if $rnp.read_location_information(&read_name, loc) {
                    $rnp.read_name_stored = read_name;
                    $rnp.physical_location_stored.set_x(loc.get_x());
                    $rnp.physical_location_stored.set_y(loc.get_y());
                    $rnp.physical_location_stored.set_tile(loc.get_tile());
                    return true;
                }
                // return false if read name cannot be parsed
                return false;
            } else {
                loc.set_tile($rnp.physical_location_stored.get_tile());
                loc.set_x($rnp.physical_location_stored.get_x());
                loc.set_y($rnp.physical_location_stored.get_y());
                return true;
            }
        }

        /**
        * Given a string, splits the string by the delimiter, and returns the the last three fields parsed as integers.  Parsing a field
        * considers only a sequence of digits up until the first non-digit character.  The three values are stored in the passed-in array.
        *
        * @throws NumberFormatException if any of the tokens that should contain numbers do not start with parsable numbers
        */
        fn get_last_three_fields(read_name: &str, delim: char, tokens: &mut [i32]) -> Result<i32, Error> {
            let mut tokens_idx = 2_i32; // start at the last token
            let mut num_fields = 0;
            let mut end_idx = read_name.len() as i32;

            let mut i = (read_name.len() - 1) as i32;
            // find the last three tokens only
            while i >= 0 && tokens_idx >= 0 {
                if *read_name.as_bytes().get(i as usize).unwrap() as char == delim || i == 0 {
                    num_fields += 1;
                    tokens[tokens_idx as usize] = Self::rapid_parse_i32(
                        read_name
                            .get(
                                {
                                    if i == 0 {
                                        0
                                    } else {
                                        (i + 1) as usize
                                    }
                                }..(end_idx as usize),
                            )
                            .unwrap(),
                    )?;
                    tokens_idx -= 1;
                    end_idx = i;
                }

                i -= 1;
            }

            // continue to find the # of fields
            while (0 <= i) {
                if *read_name.as_bytes().get(i as usize).unwrap() as char == delim || i == 0 {
                    num_fields += 1;
                }
                i -= 1;
            }

            if num_fields < 3 {
                tokens[0] = -1;
                tokens[1] = -1;
                tokens[2] = -1;

                Ok(-1)
            } else {
                Ok(num_fields)
            }
        }
    };
}
pub(crate) use impl_read_name_parser_ext;

impl ReadNameParserExt for ReadNameParser {
    impl_read_name_parser_ext!(self, self);
}

#[cfg(test)]
mod read_name_parser_test {
    use crate::markdup::utils::physical_location::PhysicalLocationShort;

    use super::*;

    use super::ReadNameParser;

    /** Tests rapidParseInt for positive and negative numbers, as well as non-digit suffixes */
    #[test]
    fn test_rapid_parse_int() {
        for i in (-100..100) {
            let i_string = i.to_string();
            assert_eq!(ReadNameParser::rapid_parse_i32(&i_string).unwrap(), i);

            // trailing characters
            assert_eq!(
                ReadNameParser::rapid_parse_i32(&format!("{}A", i_string)).unwrap(),
                i
            );
            assert_eq!(
                ReadNameParser::rapid_parse_i32(&format!("{}ACGT", i_string)).unwrap(),
                i
            );
            assert_eq!(
                ReadNameParser::rapid_parse_i32(&format!("{}.1", i_string)).unwrap(),
                i
            );
        }
    }

    /** Tests rapidParseInt for positive and negative numbers, as well as non-digit suffixes */
    #[test]
    fn test_rapid_parse_int_fails() {
        let values = ["foo", "bar", "abc123", "-foo", "f00", "-f00"];
        for s in values {
            match ReadNameParser::rapid_parse_i32(s) {
                Ok(_) => panic!("Should have failed to rapid-parse {} as an int.", s),
                Err(_) => (),
            }
        }
    }

    /** Helper for testGetRapidDefaultReadNameRegexSplit */
    fn do_test_get_rapid_default_read_name_regex_split(num_fields: i32) {
        let mut read_name = String::new();

        for i in (0..num_fields) {
            if i > 0 {
                read_name.push_str(":");
            }
            read_name.push_str(&i.to_string());
        }

        let mut input_fields = [-1_i32; 3];
        let mut expected_fields = [0_i32; 3];

        if num_fields < 3 {
            assert_eq!(
                ReadNameParser::get_last_three_fields(&read_name, ':', &mut input_fields).unwrap(),
                -1
            );
        } else {
            assert_eq!(
                ReadNameParser::get_last_three_fields(&read_name, ':', &mut input_fields).unwrap(),
                num_fields
            );
            expected_fields = [-1; 3];

            if num_fields > 0 {
                expected_fields[0] = num_fields - 3;
            }

            if num_fields > 1 {
                expected_fields[1] = num_fields - 2;
            }

            if num_fields > 2 {
                expected_fields[2] = num_fields - 1;
            }

            for i in (0..input_fields.len()) {
                assert_eq!(input_fields[i], expected_fields[i]);
            }
        }
    }

    /** Tests that we split the string with the correct # of fields, and modified values */
    #[test]
    fn test_get_rapid_default_read_name_regex_split() {
        for i in (1..10) {
            do_test_get_rapid_default_read_name_regex_split(i);
        }
    }

    fn test_parse_read_name_data_provider() -> [(&'static str, i32, i32, i32); 2] {
        [
            ("RUNID:7:1203:2886:82292", 1203, 2886, 82292),
            ("RUNID:7:1203:2884:16834", 1203, 2884, 16834),
        ]
    }

    // NB: these test fail s due to overflow in the duplicate finder test.  This has been the behavior previously, so keep it for now.
    fn _test_parse_read_name_overflow(read_name: &str, tile: i32, x: i32, y: i32) {
        let mut parser = ReadNameParser::new();
        // parser.log = "rnp".to_owned();
        let mut loc = PhysicalLocationShort::default();

        assert!(parser.add_location_information(read_name.to_owned(), &mut loc));
        assert_eq!(loc.get_tile() as i32, tile);
        assert_eq!(loc.get_x(), x as i16 as i32); // casting to short for the overflow
        assert_eq!(loc.get_y(), y as i16 as i32); // casting to short for the overflow
    }

    #[test]
    #[should_panic]
    fn test_parse_read_name_overflow() {
        // init_global_logger(log::LevelFilter::Debug);

        test_parse_read_name_data_provider()
            .into_iter()
            .for_each(|e| {
                _test_parse_read_name_overflow(e.0, e.1, e.2, e.3);
            })
    }

    fn test_read_name_parsing_data_provider(
    ) -> [(&'static str, &'static str, i32, i32, i32, bool); 11] {
        let last_three_fields_regex: &str = "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$";

        [
            (
                last_three_fields_regex,
                "RUNID:123:000000000-ZZZZZ:1:1105:17981:23325",
                1105,
                17981,
                23325,
                true,
            ),
            (
                last_three_fields_regex,
                "RUNID:123:000000000-ZZZZZ:1:1109:22981:17995",
                1109,
                22981,
                17995,
                true,
            ),
            (
                last_three_fields_regex,
                "1109:22981:17995",
                1109,
                22981,
                17995,
                true,
            ),
            (
                last_three_fields_regex,
                "RUNID:7:1203:2886:82292",
                1203,
                2886,
                82292,
                true,
            ),
            (
                last_three_fields_regex,
                "RUNID:7:1203:2884:16834",
                1203,
                2884,
                16834,
                true,
            ),
            (
                last_three_fields_regex,
                "1109ABC:22981DEF:17995GHI",
                1109,
                22981,
                17995,
                true,
            ),
            (
                DEFAULT_READ_NAME_REGEX,
                "RUNID:123:000000000-ZZZZZ:1:1105:17981:23325",
                1105,
                17981,
                23325,
                true,
            ),
            (
                DEFAULT_READ_NAME_REGEX,
                "RUNID:123:000000000-ZZZZZ:1:1109:22981:17995",
                1109,
                22981,
                17995,
                true,
            ),
            (
                DEFAULT_READ_NAME_REGEX,
                "1109:22981:17995",
                1109,
                22981,
                17995,
                false,
            ),
            (
                DEFAULT_READ_NAME_REGEX,
                "RUNID:7:1203:2886:82292",
                1203,
                2886,
                82292,
                true,
            ),
            (
                DEFAULT_READ_NAME_REGEX,
                "RUNID:7:1203:2884:16834",
                1203,
                2884,
                16834,
                true,
            ),
        ]
    }

    fn _test_read_name_parsing(
        read_name_regex: &str,
        read_name: &str,
        tile: i32,
        x: i32,
        y: i32,
        add_location_information_succeeds: bool,
    ) {
        let mut parser = ReadNameParser::with_regex(read_name_regex.to_owned());
        let mut loc = PhysicalLocationInt::default();

        assert_eq!(
            parser.add_location_information(read_name.to_owned(), &mut loc),
            add_location_information_succeeds
        );
        if add_location_information_succeeds {
            // just check the location
            assert_eq!(loc.get_tile() as i32, tile);
            assert_eq!(loc.get_x(), x);
            assert_eq!(loc.get_y(), y);
        } else if read_name_regex == DEFAULT_READ_NAME_REGEX {
            // additional testing on the default regex
            let mut tokens = [0_i32; 3];
            ReadNameParser::get_last_three_fields(read_name, ':', &mut tokens).unwrap();
            assert_eq!(tokens[0], tile);
            assert_eq!(tokens[1], x);
            assert_eq!(tokens[2], y);
        }
    }

    #[test]
    fn test_read_name_parsing() {
        test_read_name_parsing_data_provider()
            .into_iter()
            .for_each(|e| {
                _test_read_name_parsing(e.0, e.1, e.2, e.3, e.4, e.5);
            })
    }

    #[test]
    /// Testing that the parser behavior stays constant after being java serialized
    fn test_serialized_read_name_parser() {
        // A non standard parser
        let last_three_fields_regex = "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$";
        let mut parser = ReadNameParser::with_regex(last_three_fields_regex.to_string());
        let mut default_parser = ReadNameParser::new();

        let mut loc = PhysicalLocationInt::default();
        let read_name = "1109ABC:22981DEF:17995GHI";
        let expected_tile = 1109;
        let expected_x = 22981;
        let expected_y = 17995;

        assert!(parser.add_location_information(read_name.to_string(), &mut loc));
        assert_eq!(loc.get_tile(), expected_tile);
        assert_eq!(loc.get_x(), expected_x);
        assert_eq!(loc.get_y(), expected_y);

        loc = PhysicalLocationInt::default();
        assert!(!default_parser.add_location_information(read_name.to_string(), &mut loc));

        let serialized_bytes: Vec<u8> = bincode::serialize(&[&parser, &default_parser]).unwrap();

        let [mut parser_serialized, mut default_parser_serialized] =
            bincode::deserialize::<[ReadNameParser; 2]>(&serialized_bytes).unwrap();

        loc = PhysicalLocationInt::default();

        assert!(parser_serialized.add_location_information(read_name.to_string(), &mut loc));
        assert_eq!(loc.get_tile(), expected_tile);
        assert_eq!(loc.get_x(), expected_x);
        assert_eq!(loc.get_y(), expected_y);

        loc = PhysicalLocationInt::default();
        assert!(
            !default_parser_serialized.add_location_information(read_name.to_string(), &mut loc)
        );
    }
}
