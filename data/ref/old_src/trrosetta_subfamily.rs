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

use std::fmt;
use crate::geometry::core::GeometricHash;
use crate::geometry::core::HashType;
use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u64);


impl GeometricHash for HashValue {
    fn perfect_hash(&self, feature: Vec<f32>) -> Self {
        let cb_dist = feature[0];
        let omega = feature[1];
        let theta1 = feature[2];
        let theta2 = feature[3];
        let phi1 = feature[4];
        let phi2 = feature[5];
        HashValue::perfect_hash(cb_dist, omega, theta1, theta2, phi1, phi2)
    }
    fn reverse_hash(&self, hash: u64) -> Vec<f32> {
        let values = self.reverse_hash();
        vec![
            continuize_value(values[0] as u64, 2.0, 20.0, 16.0),
            continuize_value(values[1] as u64, -1.0, 1.0, 4.0),
            continuize_value(values[2] as u64, -1.0, 1.0, 8.0),
            continuize_value(values[3] as u64, -1.0, 1.0, 8.0),
            continuize_value(values[4] as u64, 0.0, 180.0, 16.0),
            continuize_value(values[5] as u64, 0.0, 180.0, 16.0),
        ]
    }
    fn hash_type(&self) -> HashType {
        HashType::TRRosettaHash
    }
}

impl HashValue {
    pub fn from_u64(hashvalue: u64) -> Self {
        HashValue(hashvalue)
    }
    pub fn as_u64(&self) -> u64 {
        self.0
    }
    pub fn as_usize(&self) -> usize {
        self.0 as usize
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
        let h_cb_dist = discretize_value(cbd, 2.0, 20.0, 16.0);
        // Torsion angles
        let h_omega = discretize_value(omega, -1.0, 1.0, 4.0);
        let h_theta1 = discretize_value(theta1, -1.0, 1.0, 8.0);
        let h_theta2 = discretize_value(theta2, -1.0, 1.0, 8.0);
        // Planar angles
        let h_phi1 = discretize_value(phi1, 0.0, 180.0, 16.0);
        let h_phi2 = discretize_value(phi2, 0.0, 180.0, 16.0);

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
