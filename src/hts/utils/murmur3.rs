pub(crate) struct Murmur3 {
    seed: u32,
}

impl Murmur3 {
    const serialVersionUID: i64 = 1;

    /** Constructs a Murmur3 hash with the given seed. */
    pub(crate) fn new(seed: u32) -> Self {
        Self { seed }
    }

    pub(crate) fn hash_bytes(&self, bytes: &[u8]) -> u32 {
        let mut h1 = self.seed;

        let len = bytes.len();

        for i in (0..len - 1).step_by(2) {
            let mut k1 = bytes[i] as u32 | ((bytes[i + 1] as u32) << 16);
            k1 = Self::mix_k1(k1);
            h1 = Self::mix_h1(h1, k1);
        }

        // deal with any remaining bytes
        if (len & 1) == 1 {
            let mut k1 = *bytes.last().unwrap() as u32;
            k1 = Self::mix_k1(k1);
            h1 ^= k1;
        }

        Self::fmix(h1, 2 * len as u32)
    }

    #[allow(non_upper_case_globals)]
    fn mix_k1(mut k1: u32) -> u32 {
        const c1: u32 = 0xcc9e2d51;
        const c2: u32 = 0x1b873593;

        k1 *= c1;
        k1 = k1.rotate_left(15);
        k1 *= c2;

        k1
    }

    fn mix_h1(mut h1: u32, k1: u32) -> u32 {
        h1 ^= k1;
        h1 = h1.rotate_left(13);
        h1 = h1 * 5 + 0xe6546b64;

        h1
    }

    fn fmix(mut h1: u32, length: u32) -> u32 {
        h1 ^= length;
        h1 ^= h1 >> 16;
        h1 *= 0x85ebca6b;
        h1 ^= h1 >> 13;
        h1 *= 0xc2b2ae35;
        h1 ^= h1 >> 16;

        h1
    }
}
