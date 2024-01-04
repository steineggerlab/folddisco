use std::fmt;
use crate::geometry::core::GeometricHash;
use crate::geometry::core::HashType;
// use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u64);

impl HashValue {
    // Constructor
    // pub fn new(dist: u16, angle: u16) -> Self {
    //     let hashvalue = dist << 8 | angle ;
    //     HashValue(hashvalue)
    // }

    pub fn from_u64(hashvalue: u64) -> Self {
        HashValue(hashvalue)
    }

    pub fn perfect_hash(dist: f32, angle: f32) -> Self {
        let hashvalue = (dist.to_bits() as u64) << 32 | angle.to_bits() as u64;
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> (f32, f32) {
        let dist_bits = (self.0 >> 32) as u32;
        let angle_bits = (self.0 & 0x00000000FFFFFFFF) as u32;
        let dist = f32::from_bits(dist_bits);
        let angle = f32::from_bits(angle_bits);
        (dist, angle)
    }

    // pub fn loose_hash(&self) -> Self {
    //     let (dist, angle) = self.reverse_hash();
    //     let dist = dist / 2; // 20 seems to be too loose
    //     let angle = angle / 5; // 20 seems to be too loose
    //     let hashvalue = dist << 8 | angle;
    //     HashValue(hashvalue)
    // }

    // pub fn
}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (dist, angle) = self.reverse_hash();
        write!(f, "HashValue({}), dist={}, angle={}", self.0, dist, angle)
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (dist, angle) = self.reverse_hash();
        write!(f, "{}\t{}\t{}", self.0, dist, angle)
        // write!(f, "{}", self.0)
    }
}

pub type HashCollection = Vec<HashValue>;

pub fn discretize_angle(val: f32) -> u16 {
    // min = -180, max = 180, bins = 2^16,
    let cont_f = 360.0_f32 / (2.0_f32.powi(16) - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    (val - (-180.0) * (disc_f) + 0.5) as u16
}

pub fn continuize_angle(val: u16) -> f32 {
    // min = -180, max = 180, bins = 2^16, disc_f = 2^16 / 360, cont_f = 360 / 2^16
    let cont_f = 360.0_f32 / (2.0_f32.powi(16) - 1.0_f32);

    (val as f32) * (cont_f) - 180.0
}

// #[derive(Debug)]
// pub struct GeometricHasher {
//     angles: Vec<f32>,
//     distances: Vec<f32>,
//     hashvalue: Vec<u64>, // TODO: Figure out the data type for hash
// }

// impl GeometricHasher {
//     pub fn new(angles: Vec<f32>, distances: Vec<f32>) -> GeometricHasher {
//         GeometricHasher {
//             angles: angles,
//             distances: distances,
//             hashvalue: Vec::new(),
//         }
//     }
// }

// impl Hasher for HashValue {
//     fn finish(&self) -> u64 {
//         self.0 as u64
//     }

//     fn write(&mut self, _bytes: &[u8]) {

//     }

//     fn write_u64(&mut self, i: u64) {
//         self.0 = i as u16;
//     }
// }


// pub trait GeometricHash {
//     fn from_u64(hash: u64) -> Self;
//     fn to_u64(&self) -> u64;
//     fn perfect_hash(&self, feature: Vec<f32>) -> u64;
//     fn reverse_hash(&self, hash: u64) -> Vec<f32>;
//     fn hash_type(&self) -> HashType;
// }

impl GeometricHash for HashValue {
    fn from_u64(hash: u64) -> Self {
        HashValue(hash)
    }

    fn to_u64(&self) -> u64 {
        self.0
    }

    fn perfect_hash(feature: Vec<f32>) -> Self {
        let dist = feature[0];
        let angle = feature[1];
        let hashvalue = (dist.to_bits() as u64) << 32 | angle.to_bits() as u64;
        HashValue(hashvalue)
    }

    fn reverse_hash(&self) -> Vec<f32> {
        let dist_bits = (self.0 >> 32) as u32;
        let angle_bits = (self.0 & 0x00000000FFFFFFFF) as u32;
        let dist = f32::from_bits(dist_bits);
        let angle = f32::from_bits(angle_bits);
        vec![dist, angle]
    }

    fn hash_type(&self) -> HashType {
        HashType::SimpleHash
    }
}