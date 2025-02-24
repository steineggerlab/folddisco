// File: pdb_tr.rs
// Created: 2024-03-27 17:35:35
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved
// PDB motif + 2 torsion angles

use std::fmt;
use crate::geometry::core::HashType;
use crate::utils::convert::discretize_f32_value_into_u32 as discretize_value;
use crate::utils::convert::continuize_u32_value_into_f32 as continuize_value;
use crate::utils::convert::*;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(pub u32);

pub const HYBRID_NBIN_DIST: f32 = 16.0;
pub const HYBRID_NBIN_SIN_COS: f32 = 4.0;

impl HashValue {
    #[inline]
    pub fn perfect_hash(feature: &Vec<f32>, nbin_dist: usize, nbin_angle: usize) -> u32 {
        let nbin_dist = if nbin_dist > 16 {
            16.0
        } else if nbin_dist == 0 {
            HYBRID_NBIN_DIST
        } else {
            nbin_dist as f32
        };
        let nbin_angle = if nbin_angle > 4 {
            4.0
        } else if nbin_angle == 0 {
            HYBRID_NBIN_SIN_COS
        } else {
            nbin_angle as f32
        };
        let res1_group = feature[0] as u32;
        let res2_group = feature[1] as u32;
        let ca_dist = discretize_value(
            feature[2], MIN_DIST, MAX_DIST, nbin_dist
        );
        let cb_dist = discretize_value(
            feature[3], MIN_DIST, MAX_DIST, nbin_dist
        );
        // Angle is expected to be in radians
        let sin_ca_cb_angle = feature[4].sin();
        let cos_ca_cb_angle = feature[4].cos();
        let sin_ca_cb_angle = discretize_value(
            sin_ca_cb_angle, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_ca_cb_angle = discretize_value(
            cos_ca_cb_angle, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        // Two torsion angles: 
        let sin_phi1 = feature[5].sin();
        let cos_phi1 = feature[5].cos();
        let sin_phi2 = feature[6].sin();
        let cos_phi2 = feature[6].cos();
        let sin_phi1 = discretize_value(
            sin_phi1, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_phi1 = discretize_value(
            cos_phi1, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let sin_phi2 = discretize_value(
            sin_phi2, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_phi2 = discretize_value(
            cos_phi2, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        
        let sin_bb_phi1 = feature[7].sin();
        let cos_bb_phi1 = feature[7].cos();
        let sin_bb_phi2 = feature[8].sin();
        let cos_bb_phi2 = feature[8].cos();
        let sin_bb_phi1 = discretize_value(
            sin_bb_phi1, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_bb_phi1 = discretize_value(
            cos_bb_phi1, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let sin_bb_phi2 = discretize_value(
            sin_bb_phi2, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_bb_phi2 = discretize_value(
            cos_bb_phi2, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        
        let hashvalue = res1_group << 30 | res2_group << 28 | ca_dist << 24
            | cb_dist << 20 | sin_ca_cb_angle << 18 | cos_ca_cb_angle << 16
            | sin_phi1 << 14 | cos_phi1 << 12 | sin_phi2 << 10 | cos_phi2 << 8 
            | sin_bb_phi1 << 6 | cos_bb_phi1 << 4 | sin_bb_phi2 << 2 | cos_bb_phi2;
        hashvalue
    }
    
    pub fn perfect_hash_default(feature: &Vec<f32>) -> u32 {
        HashValue::perfect_hash(feature, HYBRID_NBIN_DIST as usize, HYBRID_NBIN_SIN_COS as usize)
    }
    
    pub fn reverse_hash_default(&self) -> [f32; 9] {
        self.reverse_hash(HYBRID_NBIN_DIST as usize, HYBRID_NBIN_SIN_COS as usize)
    }
    
    pub fn reverse_hash(&self, nbin_dist: usize, nbin_angle: usize) -> [f32; 9] {
        let res1_group = ((self.0 >> 30) & BITMASK32_2BIT)as f32;
        let res2_group = ((self.0 >> 28) & BITMASK32_2BIT) as f32;
        let ca_dist = continuize_value(
            (self.0 >> 24) & BITMASK32_4BIT as u32, 
            MIN_DIST, MAX_DIST, nbin_dist as f32
        );
        let cb_dist = continuize_value(
            (self.0 >> 20) & BITMASK32_4BIT as u32,
            MIN_DIST, MAX_DIST, nbin_dist as f32
        );
        let sin_ca_cb_angle = continuize_value(
            (self.0 >> 18) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_ca_cb_angle = continuize_value(
            (self.0 >> 16) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let sin_phi1 = continuize_value(
            (self.0 >> 14) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_phi1 = continuize_value(
            (self.0 >> 12) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let sin_phi2 = continuize_value(
            (self.0 >> 10) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_phi2 = continuize_value(
            (self.0 >> 8) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let sin_bb_phi1 = continuize_value(
            (self.0 >> 6) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_bb_phi1 = continuize_value(
            (self.0 >> 4) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let sin_bb_phi2 = continuize_value(
            (self.0 >> 2) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_bb_phi2 = continuize_value(
            self.0 & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        

        let ca_cb_angle = sin_ca_cb_angle.atan2(cos_ca_cb_angle).to_degrees();
        let phi1 = sin_phi1.atan2(cos_phi1).to_degrees();
        let phi2 = sin_phi2.atan2(cos_phi2).to_degrees();
        let bb_phi1 = sin_bb_phi1.atan2(cos_bb_phi1).to_degrees();
        let bb_phi2 = sin_bb_phi2.atan2(cos_bb_phi2).to_degrees();
        
        [res1_group, res2_group, ca_dist, cb_dist, ca_cb_angle, phi1, phi2, bb_phi1, bb_phi2]
    }
    
    pub fn hash_type(&self) -> HashType {
        HashType::Hybrid
    }
    
    pub fn from_u32(hashvalue: u32) -> Self {
        HashValue(hashvalue)
    }
    
    pub fn as_u32(&self) -> u32 {
        self.0
    }
    
    pub fn from_u64(hashvalue: u64) -> Self {
        HashValue(hashvalue as u32)
    }
    
    pub fn as_u64(&self) -> u64 {
        self.0 as u64
    }
    
    pub fn is_symmetric(&self) -> bool {
        let values = self.reverse_hash_default();
        // Residue pair is symmetric and phi is symmetric
        (values[0] == values[1]) && (values[5] == values[6])
    }
}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let values = self.reverse_hash_default();
        write!(f, "HashValue({}), values={:?}", self.0, values)
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let values = self.reverse_hash_default();
        write!(f, "{}\t{:?}", self.0, values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::core::GeometricHash;
    use crate::utils::convert::map_aa_to_u8_group;
    #[test]
    fn test_geometrichash_works() {
        // Test perfect hash
        let raw_feature = (
            b"PHE", b"VAL", 14.0_f32, 15.9_f32, 116.0_f32, 80.0_f32, -100.0_f32, 40.0_f32, 30.0_f32
        );
        let raw_feature2 = (
            b"VAL", b"PHE", 15.9_f32, 14.0_f32, 116.0_f32, 80.0_f32, -100.0_f32, -120.0_f32, -30.0_f32
        );
        let raw_feature = vec![
            map_aa_to_u8_group(raw_feature.0) as f32, map_aa_to_u8_group(raw_feature.1) as f32,
            raw_feature.2, raw_feature.3, raw_feature.4.to_radians(),
            raw_feature.5.to_radians(), raw_feature.6.to_radians(),
            raw_feature.7.to_radians(), raw_feature.8.to_radians()
        ];
        let raw_feature2 = vec![
            map_aa_to_u8_group(raw_feature2.0) as f32, map_aa_to_u8_group(raw_feature2.1) as f32,
            raw_feature2.2, raw_feature2.3, raw_feature2.4.to_radians(),
            raw_feature2.5.to_radians(), raw_feature2.6.to_radians(),
            raw_feature2.7.to_radians(), raw_feature2.8.to_radians()
        ];
        let start = std::time::Instant::now();
        for _ in 0..10000 {
            let _ = HashValue::perfect_hash_default(&raw_feature);
            let _ = HashValue::perfect_hash_default(&raw_feature2);
        }
        let duration = start.elapsed();
        println!("Time elapsed in perfect_hash_default() is: {:?}", duration);
        let hash = HashValue::perfect_hash_default(&raw_feature);
        let hash = GeometricHash::from_u32(hash, HashType::Hybrid);
        match hash {
            GeometricHash::Hybrid(hash) => {
                println!("{:?}", hash);
            },
            _ => panic!("Invalid hash type"),
        }
    }
}