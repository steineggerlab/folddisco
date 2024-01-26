use std::fmt;
use crate::geometry::core::GeometricHash;
use crate::geometry::core::HashType;
// use std::hash::Hasher;

// Residue 1: 5 bits
// Residue 2: 5 bits
// Distance1: 16 bins 4 bits
// Distance delta: 16 bins 4 bits
// Angle: 32 bins 5 bits
// total: 23 bits

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u32);

impl HashValue {
    pub fn from_u32(hashvalue: u32) -> Self {
        HashValue(hashvalue)
    }

    pub fn hash_type(&self) -> HashType {
        HashType::PDBMotif
    }

    pub fn perfect_hash(aa1: u8, aa2: u8, ca_dist: f32, delta_dist: f32, angle: f32) -> Self {
        let aa1_hash = aa1 as u32;
        let aa2_hash = aa2 as u32;
        assert!(i != j);

        let h4 = discretize_value(ca_dist, 2.0, 20.0, 16.0);
        
        
        assert!(aa1_hash < 256);
        assert!(aa2_hash < 256);

        assert!(h4 < 256);
        let hashvalue = aa1_hash << 24 | aa2_hash << 16 | h3 << 8 | h4;
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



// Point Pair Features (PPF):

use std::fmt;
use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u32);

impl HashValue {
    pub fn from_u32(hashvalue: u32) -> Self {
        HashValue(hashvalue)
    }

    pub fn perfect_hash(aa1: u8, aa2: u8, i: usize, j: usize, dist: f32) -> Self {
        let h1 = aa1 as u32;
        let h2 = aa2 as u32;
        assert!(i != j);
        let h3 = ((i as i32 - j as i32).abs() + 1).ilog2();
        let h4 = discretize_value(dist, 2.0, 20.0, 4.0);

        assert!(h1 < 256);
        assert!(h2 < 256);
        assert!(h3 < 256);
        assert!(h4 < 256);
        let hashvalue = h1 << 24 | h2 << 16 | h3 << 8 | h4;
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> [f32; 4] {
        let h1 = (self.0 >> 24) as u8;
        let h2 = (self.0 >> 16) as u8;
        let h3 = (self.0 >> 8) as u8;
        let h4 = (self.0 & 0x000000FF) as u8;
        [h1 as f32, h2 as f32, h3 as f32, h4 as f32]
    }

    pub fn neighbors(&self) -> Vec<HashValue> {
        let mut neighbors = Vec::new();
        let mut hash = self.0;
        // Append hashes with 1 bit difference for h3 & h4
        // h3-1, h3, h3+1; h4-1, h4, h4+1 (if possible) total 9 hashes
        // Explain code with comments
        for i in 0..9 {
            let mut h3 = (hash >> 8) as u8;
            let mut h4 = (hash & 0x000000FF) as u8;
            // h3-1, h3, h3+1; h4-1, h4, h4+1 (if possible) total 9 hashes
            // Explain code with comments
            if i < 3 {
                h3 = h3.wrapping_sub(1);
            } else if i > 5 {
                h3 = h3.wrapping_add(1);
            }
            if i % 3 == 0 {
                h4 = h4.wrapping_sub(1);
            } else if i % 3 == 2 {
                h4 = h4.wrapping_add(1);
            }
            // Check if h3 and h4 are valid

            let hashvalue = (h3 as u32) << 8 | h4 as u32;
            neighbors.push(HashValue(hashvalue));

        }
        neighbors
    }
    
}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let feature = self.reverse_hash();
        write!(
            f,
            "HashValue({}), aa1={}, aa2={}, logdist={}, dist={}",
            self.0, feature[0], feature[1], feature[2], feature[3]
        )
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let feature  = self.reverse_hash();
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.0, feature[0], feature[1], feature[2], feature[3]
        )
    }
}

pub fn discretize_value(val: f32, min: f32, max: f32, num_bin: f32) -> u32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    (val * (disc_f) + 0.5) as u32
}


pub type HashCollection = Vec<HashValue>;

pub fn map_aa_to_u8(aa: &[u8; 3]) -> u8 {
    match aa {
        b"ALA" => 0,
        b"ARG" => 1,
        b"ASN" => 2,
        b"ASP" => 3,
        b"CYS" => 4,
        b"GLN" => 5,
        b"GLU" => 6,
        b"GLY" => 7,
        b"HIS" => 8,
        b"ILE" => 9,
        b"LEU" => 10,
        b"LYS" => 11,
        b"MET" => 12,
        b"PHE" => 13,
        b"PRO" => 14,
        b"SER" => 15,
        b"THR" => 16,
        b"TRP" => 17,
        b"TYR" => 18,
        b"VAL" => 19,
        _ => panic!("Invalid AA"),
    }
}
pub fn map_u8_to_aa(aa: u8) -> &'static str {
    match aa {
        0 => "ALA",
        1 => "ARG",
        2 => "ASN",
        3 => "ASP",
        4 => "CYS",
        5 => "GLN",
        6 => "GLU",
        7 => "GLY",
        8 => "HIS",
        9 => "ILE",
        10 => "LEU",
        11 => "LYS",
        12 => "MET",
        13 => "PHE",
        14 => "PRO",
        15 => "SER",
        16 => "THR",
        17 => "TRP",
        18 => "TYR",
        19 => "VAL",
        _ => panic!("Invalid AA"),
    }
}

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