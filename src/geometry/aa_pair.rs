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