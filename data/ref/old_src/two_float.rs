use std::fmt;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u64);

impl HashValue {
    pub fn from_u64(hashvalue: u64) -> Self {
        HashValue(hashvalue)
    }

    pub fn perfect_hash(val1: f32, val2: f32) -> Self {
        // Save two f32 as is in a u64
        let hashvalue = (val1.round().to_bits() as u64) << 32 | (val2.round().to_bits() as u64);
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> [f32; 2] {
        let val1 = f32::from_bits((self.0 >> 32) as u32);
        let val2 = f32::from_bits((self.0 & 0x00000000FFFFFFFF) as u32);
        [val1, val2]
    }
}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let values = self.reverse_hash();
        write!(f, "HashValue({}), values={:?}", self.0, values)
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let val = self.reverse_hash();
        write!(f, "{}\t{}\t{}", self.0, val[0], val[1])
    }
}

pub type HashCollection = Vec<HashValue>;

pub fn normalize_angle_degree(val: f32, min: f32, max: f32) -> f32 {
    (val - min) / (max - min)
}
