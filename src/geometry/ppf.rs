// Point Pair Features (PPF):

use std::fmt;
use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u32);

impl HashValue {
    pub fn from_u32(hashvalue: u32) -> Self {
        HashValue(hashvalue)
    }

    pub fn perfect_hash(ppf: &[f32; 4]) -> Self {
        let h0 = ppf[0].round() as u32 / 2;
        let h1 = discretize_angle(ppf[1]) as u32;
        let h2 = discretize_angle(ppf[2]) as u32;
        let h3 = discretize_angle(ppf[3]) as u32;
        assert!(h0 < 256);
        assert!(h1 < 256);
        assert!(h2 < 256);
        assert!(h3 < 256);
        let hashvalue = h0 << 24 | h1 << 16 | h2 << 8 | h3;
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> [f32; 4] {
        let h0 = (self.0 >> 24) as u8;
        let h1 = (self.0 >> 16) as u8;
        let h2 = (self.0 >> 8) as u8;
        let h3 = (self.0 & 0x000000FF) as u8;
        [(h0 * 2) as f32, h1 as f32, h2 as f32, h3 as f32]
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
        let ppf = self.reverse_hash();
        write!(f, "HashValue({}), dist={}, angle={},{},{}", self.0, ppf[0], ppf[1], ppf[2], ppf[3])
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ppf = self.reverse_hash();
        write!(f, "{}\t{}\t{}\t{}\t{}", self.0, ppf[0], ppf[1], ppf[2], ppf[3])
        // write!(f, "{}", self.0)
    }
}

pub type HashCollection = Vec<HashValue>;

pub fn discretize_angle(val: f32) -> u8 {
    // min = 0, max = 180, bins = 2^16,
    let cont_f = 180.0_f32 / (6.0_f32 - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    (val * (disc_f) + 0.5) as u8
}

pub fn continuize_angle(val: u8) -> f32 {
    // min = 0, max = 180, bins = 2^16, disc_f = 2^16 / 360, cont_f = 360 / 2^16
    let cont_f = 180.0_f32 / (6.0_f32 - 1.0_f32);
    (val as f32) * (cont_f)
}

// pub fn discretize_angle(val: f32) -> u16 {
//     // min = -180, max = 180, bins = 2^16,
//     let cont_f = 360.0_f32 / (2.0_f32.powi(16) - 1.0_f32);
//     let disc_f = 1.0_f32 / cont_f;
//     (val - (-180.0) * (disc_f) + 0.5) as u16
// }

// pub fn continuize_angle(val: u16) -> f32 {
//     // min = -180, max = 180, bins = 2^16, disc_f = 2^16 / 360, cont_f = 360 / 2^16
//     let cont_f = 360.0_f32 / (2.0_f32.powi(16) - 1.0_f32);

//     (val as f32) * (cont_f) - 180.0
// }


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
