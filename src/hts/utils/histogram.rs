use std::{borrow::Borrow, collections::BTreeMap, fmt::Display, hash::Hash, ops::Deref};

/**
 * Class for computing and accessing histogram type data.  Stored internally in
 * a sorted Map so that keys can be iterated in order.
 *
 * @author Tim Fennell
 */
#[derive(Clone, Debug)]
pub(crate) struct Histogram<N> {
    map: BTreeMap<N, Bin<N>>,
    bin_label: String,
    value_label: String,
}

impl<N> Default for Histogram<N>
where
    N: PartialEq,
{
    fn default() -> Self {
        Self {
            map: Default::default(),
            bin_label: "BIN".to_string(),
            value_label: "VALUE".to_string(),
        }
    }
}

impl<N> Histogram<N>
where
    N: PartialEq + Ord + Copy + Into<f64>,
{
    const SERIAL_VERSION_UID: i64 = 1;

    /** Constructs a new Histogram with default bin and value labels. */
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /** Constructs a new Histogram with supplied bin and value labels. */
    pub(crate) fn from_labels<S: ToString>(bin_label: S, value_label: S) -> Self {
        Self {
            bin_label: bin_label.to_string(),
            value_label: value_label.to_string(),
            ..Default::default()
        }
    }

    /** Constructs a new Histogram that'll use the supplied comparator to sort keys. */
    pub(crate) fn with_comparator(comparator: fn()) -> Self {
        unimplemented!("There's no way to use Comparator with BTreeMap for now.")
    }

    /** Constructor that takes labels for the bin and values and a comparator to sort the bins. */
    pub(crate) fn with_options(bin_label: String, value_label: String, comparator: fn()) -> Self {
        unimplemented!("There's no way to use Comparator with BTreeMap for now.")
    }

    /** Prefill the histogram with the supplied set of bins. */
    pub(crate) fn prefill_bins(&mut self, ids: impl Iterator<Item = N>) {
        ids.for_each(|id| {
            self.map.insert(id, Bin::new(id));
        })
    }

    /** Increments the value in the designated bin by the supplied increment. */
    pub(crate) fn increment(&mut self, id: N, increment: f64) {
        match self.map.entry(id) {
            std::collections::btree_map::Entry::Vacant(ent) => {
                ent.insert(Bin::new(id));
            }
            std::collections::btree_map::Entry::Occupied(mut ent) => {
                ent.get_mut().value += increment
            }
        }
    }

    /** Increments the value in the designated bin by the supplied increment. */
    pub(crate) fn increment1(&mut self, id: impl Into<N>) {
        self.increment(id.into(), 1.0);
    }

    pub(crate) fn get_bin_label(&self) -> &str {
        self.bin_label.as_str()
    }

    pub(crate) fn set_bin_label(&mut self, bin_label: String) {
        self.bin_label = bin_label;
    }

    pub(crate) fn get_value_label(&self) -> &str {
        self.value_label.as_str()
    }

    pub(crate) fn set_value_label(&mut self, value_label: String) {
        self.value_label = value_label;
    }

    /**
     * Assuming that the key type for the histogram is a Number type, returns the mean of
     * all the items added to the histogram.
     */
    pub(crate) fn get_mean(&self) -> f64 {
        // Could use simply getSum() / getCount(), but that would require iterating over the
        // values() set twice, which seems inefficient given how simply the computation is.
        let (mut product, mut total_count) = (0.0, 0.0);
        for bin in self.map.values() {
            let id_value = bin.get_id_value();
            let count = bin.get_value();

            product += id_value * count;
            total_count += count;
        }

        product / total_count
    }

    /**
     * Returns the sum of the products of the histgram bin ids and the number of entries in each bin.
     * Note: This is only supported if this histogram stores instances of Number.
     */
    pub(crate) fn get_sum(&self) -> f64 {
        self.map.values().fold(0.0, |mut a, b| {
            a += b.get_value() * b.get_id_value();
            a
        })
    }

    /**
     * Returns the sum of the number of entries in each bin.
     */
    pub(crate) fn get_sum_of_values(&self) -> f64 {
        self.map.values().fold(0.0, |mut a, b| {
            a += b.get_value();
            a
        })
    }

    pub(crate) fn get_standard_deviation(&self) -> f64 {
        let mean = self.get_mean();

        let mut count = 0.0;
        let mut total = 0.0;

        self.map.values().for_each(|bin| {
            let local_count = bin.get_value();
            let value = bin.get_id_value();

            count += local_count;
            total += local_count * (value - mean).powf(2.0);
        });

        (total / (count - 1.0)).sqrt()
    }

    /**
     * Calculates the mean bin size
     */
    pub(crate) fn get_mean_bin_size(&self) -> f64 {
        self.get_sum_of_values() / self.size() as f64
    }

    /**
     * Returns the size of this histogram.
     */
    pub(crate) fn size(&self) -> usize {
        self.map.len()
    }

    /**
     * Returns the comparator used to order the keys in this histogram, or
     * {@code null} if this histogram uses the {@linkplain Comparable
     * natural ordering} of its keys.
     *
     * @return the comparator used to order the keys in this histogram,
     *         or {@code null} if this histogram uses the natural ordering
     *         of its keys
     */
    pub(crate) fn comparator(&self) {
        unimplemented!("There's no way to get comparator from `BTreeMap` in Rust for now.")
    }

    /**
     * Calculates the median bin size
     */
    pub(crate) fn get_median_bin_size(&self) -> f64 {
        if self.size() == 0 {
            return 0.0;
        }

        self.get_sum_of_values() / self.size() as f64
    }

    /**
     * Retrieves the bin associated with the given key.
     */
    pub(crate) fn get<K>(&self, key: &K) -> Option<&Bin<N>>
    where
        N: Borrow<K>,
        K: Ord,
    {
        self.map.get(key)
    }
}

impl<N> PartialEq for Histogram<N>
where
    N: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.bin_label == other.bin_label
            && self.value_label == other.value_label
            && self.map == other.map
    }
}

impl<N> Display for Histogram<N>
where
    N: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#?}", self.map)
    }
}

impl<N> Hash for Histogram<N>
where
    N: Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.bin_label.hash(state);
        self.value_label.hash(state);
        self.map.hash(state);
    }
}

/** Represents a bin in the Histogram. */
#[derive(Clone, Debug)]
pub(crate) struct Bin<N> {
    id: N,
    value: f64,
}

impl<N> Bin<N>
where
    N: PartialEq + Copy,
{
    const SERIAL_VERSION_UID: i64 = 1;

    pub(crate) fn new(id: N) -> Self {
        Self { id, value: 0.0 }
    }

    /** Gets the ID of this bin. */
    pub(crate) fn get_id(&self) -> N {
        self.id
    }

    /** Gets the value in the bin. */
    pub(crate) fn get_value(&self) -> f64 {
        self.value
    }
}

// pub(crate) trait IsNumber {}

// macro_rules! impl_is_number {
//     ($ty:ty) => {
//         impl IsNumber for $ty {}
//     };
// }

// impl_is_number!(i16);
// impl_is_number!(i32);
// impl_is_number!(i64);
// impl_is_number!(u16);
// impl_is_number!(u32);
// impl_is_number!(u64);

pub(crate) trait BinExt {
    fn get_id_value(&self) -> f64;
}

impl<N> BinExt for Bin<N>
where
    N: Copy + Into<f64>,
{
    fn get_id_value(&self) -> f64 {
        self.id.into()
    }
}

impl<N> Display for Bin<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl<N> PartialEq for Bin<N>
where
    N: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value && self.id == other.id
    }
}

impl<N> Hash for Bin<N>
where
    N: Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
        self.value.to_bits().hash(state);
    }

    // fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
    //     let result = self.id.hash(state);
    //     let temp:u64 = if self.value != 0.0 {
    //         self.value.to_bits()
    //     } else {
    //         0
    //     };
    //     let result = 31 * result + (temp ^ temp.wrapping_shr(32))

    // }
}

/// A f64 wrapper `Ord` is implemented for.
/// if `PartialOrd::partial_cmp` is `None`, this program will stop running.
#[derive(Copy, Clone)]
#[allow(non_camel_case_types)]
pub(crate) struct f64H(f64);

impl Deref for f64H {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl PartialEq for f64H {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl PartialOrd for f64H {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Eq for f64H {}

impl Ord for f64H {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.0.partial_cmp(&other.0) {
            Some(ord) => ord,
            None => panic!(
                "Code failed. Failed to compare f64 values. lhs={} != rhs={}. \
                At least one of the values must not be used in `Histogram`.",
                self.0, other.0
            ),
        }
    }
}

impl Into<f64> for f64H {
    fn into(self) -> f64 {
        self.0
    }
}

impl From<f64> for f64H {
    fn from(value: f64) -> Self {
        Self(value)
    }
}
