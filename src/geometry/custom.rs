// Geometric features from trRosetta paper
// https://doi.org/10.1073/pnas.1914677117
// Paper bin size: 0.5A, 15 degree
// Original Features
// 1. Cb-Cb distance (2 - 20A)            36 bins
// 2. 3 dihedrals
// - omega: between ca1-cb1-cb2-ca2 (-180 ~ 180) 24 bins
// - theta1: between n1-ca1-cb1-cb2 (-180 ~ 180)  24 bins
// - theta2: between cb1-cb2-ca2-n2 (-180 ~ 180)  24 bins
// 3. 2 planar angles
// - phi1: between ca1-cb1-cb2 (0 ~ 180)        12 bins
// - phi2: between cb1-cb2-ca2 (0 ~ 180)        12 bins
// 36 + 24 + 24 + 24 + 12 + 12 = 132 bins

// Implementation
// Dist - 1 bin = 1A; u8
// Dihedral 16 bins = 22.5 degree; u8 * 3
// Planar 8 bins = 22.5 degree; u8 * 2
// 16^5 * 20;

        // let h_cb_dist = discretize_value(cbd, 2.0, 20.0, 8.0); 3
        // // Torsion angles
        // let h_omega = discretize_value(omega, -1.0, 1.0, 4.0);
        // let h_theta1 = discretize_value(theta1, -1.0, 1.0, 4.0);
        // let h_theta2 = discretize_value(theta2, -1.0, 1.0, 4.0);
        // // Planar angles
        // let h_phi1 = discretize_value(phi1, 0.0, 180.0, 8.0);
        // let h_phi2 = discretize_value(phi2, 0.0, 180.0, 8.0);

use std::fmt;
use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u64);

impl HashValue {
    pub fn from_u64(hashvalue: u64) -> Self {
        HashValue(hashvalue)
    }
    pub fn as_u64(&self) -> u64 {
        self.0
    }

    pub fn perfect_hash(
        cb_dist: f32,
        omega: f32,
        theta1: f32,
        theta2: f32,
        phi1: f32,
        phi2: f32,
    ) -> Self {
        let mut cbd = cb_dist;
        if cb_dist > 20.0 {
            cbd = 20.0;
        }
        let h_cb_dist = discretize_value(cbd, 2.0, 20.0, 8.0);
        // Torsion angles
        let h_omega = discretize_value(omega, -1.0, 1.0, 4.0);
        let h_theta1 = discretize_value(theta1, -1.0, 1.0, 4.0);
        let h_theta2 = discretize_value(theta2, -1.0, 1.0, 4.0);
        // Planar angles
        let h_phi1 = discretize_value(phi1, 0.0, 180.0, 8.0);
        let h_phi2 = discretize_value(phi2, 0.0, 180.0, 8.0);

        assert!(h_cb_dist < 256);
        assert!(h_omega < 256);
        assert!(h_theta1 < 256);
        assert!(h_theta2 < 256);
        assert!(h_phi1 < 256);
        assert!(h_phi2 < 256);
        let hashvalue = h_cb_dist << 40
            | h_omega << 32
            | h_theta1 << 24
            | h_theta2 << 16
            | h_phi1 << 8
            | h_phi2;
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> [f32; 6] {
        let h_cb_dist = (self.0 >> 40) as u8;
        let h_omega = (self.0 >> 32) as u8;
        let h_theta1 = (self.0 >> 24) as u8;
        let h_theta2 = (self.0 >> 16) as u8;
        let h_phi1 = (self.0 >> 8) as u8;
        let h_phi2 = (self.0 & 0x000000FF) as u8;
        [
            h_cb_dist as f32,
            h_omega as f32,
            h_theta1 as f32,
            h_theta2 as f32,
            h_phi1 as f32,
            h_phi2 as f32,
        ]
    }

    pub fn neighbors(&self, include_self: bool) -> Vec<HashValue> {
        let mut neighbors = Vec::new();
        // Get neighbors for each feature
        // Just add 1 or subtract 1 from each feature
        // If the feature is at the boundary, don't add or subtract
        let hash = self.0;
        // for i in 0..1 { // Neighbor only for cb distance
        for i in 0..6 {
            let mut h = (hash >> (40 - i * 8)) as u8;
            // Add 1 or subtract 1 checking boundary
            if h < 255u8 {
                h = h.wrapping_add(1);
                // Update hash using h
                let hashvalue = (h as u64) << (40 - i * 8);
                let hashvalue = hashvalue | (hash & 0x00FFFFFFFFFFFFFF);
                neighbors.push(HashValue(hashvalue));
            }
            if h > 0u8 {
                h = h.wrapping_sub(1);
                let hashvalue = (h as u64) << (40 - i * 8);
                let hashvalue = hashvalue | (hash & 0x00FFFFFFFFFFFFFF);
                neighbors.push(HashValue(hashvalue));
            }
        }
        if include_self {
            neighbors.push(*self);
        }
        neighbors
    }
    
}

// impl Hasher for HashValue {
//     fn finish(&self) -> u64 {
//         self.0
//     }

//     fn write(&mut self, _bytes: &[u8]) {
//         unimplemented!()
//     }
// }

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let values = self.reverse_hash();
        write!(f, "HashValue({}), values={:?}", self.0, values)
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let val = self.reverse_hash();
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.0, val[0], val[1], val[2], val[3], val[4], val[5]
        )
        // write!(f, "{}", self.0)
    }
}

pub type HashCollection = Vec<HashValue>;

pub fn discretize_value(val: f32, min: f32, max: f32, num_bin: f32) -> u64 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    (val * (disc_f) + 0.5) as u64
}

pub fn continuize_value(val: u64, min: f32, max: f32, num_bin: f32) -> f32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    (val as f32) * (cont_f)
}

pub fn normalize_angle_degree(val: f32, min: f32, max: f32) -> f32 {
    (val - min) / (max - min)
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