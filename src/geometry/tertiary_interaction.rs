// Hash based on feature to build 3Di character in Foldseek
// Main features are angles from neighboring C-alpha atoms

// TODO: Implement this module




// File: pdb_tr.rs
// Created: 2024-03-27 17:35:35
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved
// PDB motif + 2 torsion angles

use std::fmt;
use crate::geometry::core::HashType;
use crate::geometry::util::discretize_f32_value_into_u32 as discretize_value;
use crate::geometry::util::continuize_u32_value_into_f32 as continuize_value;
use crate::geometry::util::*;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(pub u32);

impl HashValue {
    pub fn perfect_hash(feature: Vec<f32>, nbin_dist: usize, nbin_angle: usize) -> Self {
        // Added one more quantization for distance
        let nbin_dist = if nbin_dist > 16 { 16.0 } else { nbin_dist as f32 };
        let nbin_angle = if nbin_angle > 8 { 8.0 } else { nbin_angle as f32 };
        
        let cos_phi_12 = feature[0].cos();
        let cos_phi_12 = discretize_value(
            cos_phi_12, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_phi_34 = feature[1].cos();
        let cos_phi_34 = discretize_value(
            cos_phi_34, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_phi_15 = feature[2].cos();
        let cos_phi_15 = discretize_value(
            cos_phi_15, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_phi_35 = feature[3].cos();
        let cos_phi_35 = discretize_value(
            cos_phi_35, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_phi_14 = feature[4].cos();
        let cos_phi_14 = discretize_value(
            cos_phi_14, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_phi_23 = feature[5].cos();
        let cos_phi_23 = discretize_value(
            cos_phi_23, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_phi_13 = feature[6].cos();
        let cos_phi_13 = discretize_value(
            cos_phi_13, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let ca_dist = discretize_value(
            feature[7], MIN_DIST, MAX_DIST, nbin_dist
        );
        // let seq_dist = if feature[8] < -4.0 {
        //     0_u32
        // } else if feature[8] > 4.0 {
        //     8_u32
        // } else {
        //     feature[8] as u32 + 4
        // };
        let seq_dist = 0u32;

        let hashvalue = (cos_phi_12 << 26) | (cos_phi_34 << 23) | (cos_phi_15 << 20) |
                        (cos_phi_35 << 17) | (cos_phi_14 << 14) | (cos_phi_23 << 11) |
                        (cos_phi_13 << 8) | (ca_dist << 4) | (seq_dist);
        HashValue(hashvalue)
    }

    pub fn perfect_hash_default(feature: Vec<f32>) -> Self {
       Self::perfect_hash(feature, NBIN_DIST as usize, NBIN_SIN_COS as usize)
    }
    
    pub fn reverse_hash_default(&self) -> Vec<f32> {
        self.reverse_hash(NBIN_DIST as usize, NBIN_SIN_COS as usize)
    }
    
    pub fn reverse_hash(&self, nbin_dist: usize, nbin_angle: usize) -> Vec<f32> {
        let nbin_dist = if nbin_dist > 16 { 16.0 } else { nbin_dist as f32 };
        let nbin_angle = if nbin_angle > 8 { 8.0 } else { nbin_angle as f32 };
        let cos_phi_12 = continuize_value(
            (self.0 >> 26) & BITMASK32_3BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        ).acos().to_degrees();
        let cos_phi_34 = continuize_value(
            (self.0 >> 23) & BITMASK32_3BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        ).acos().to_degrees();
        let cos_phi_15 = continuize_value(
            (self.0 >> 20) & BITMASK32_3BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        ).acos().to_degrees();
        let cos_phi_35 = continuize_value(
            (self.0 >> 17) & BITMASK32_3BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        ).acos().to_degrees();
        let cos_phi_14 = continuize_value(
            (self.0 >> 14) & BITMASK32_3BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        ).acos().to_degrees();
        let cos_phi_23 = continuize_value(
            (self.0 >> 11) & BITMASK32_3BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        ).acos().to_degrees();
        let cos_phi_13 = continuize_value(
            (self.0 >> 8) & BITMASK32_3BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        ).acos().to_degrees();
        let ca_dist = continuize_value(
            (self.0 >> 4) & BITMASK32_4BIT, MIN_DIST, MAX_DIST, nbin_dist
        );
        let seq_dist = (self.0 & BITMASK32_4BIT) as f32 - 4.0;

        vec![
            cos_phi_12, cos_phi_34, cos_phi_15, cos_phi_35, cos_phi_14, cos_phi_23,
            cos_phi_13, ca_dist, seq_dist
        ]

    }
    
    pub fn hash_type(&self) -> HashType {
        HashType::TertiaryInteraction
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
    use crate::geometry::{core::GeometricHash, util::map_aa_to_u8};
    #[test]
    fn test_geometrichash_works() {
        // Test perfect hash
        let raw_feature = (
            10.0_f32, 20.0_f32, 30.0_f32, 40.0_f32, 50.0_f32, 60.0_f32, 70.0_f32, 3.0_f32, 4.0_f32
        );
        let raw_feature = vec![
            raw_feature.0.to_radians(), raw_feature.1.to_radians(), raw_feature.2.to_radians(),
            raw_feature.3.to_radians(), raw_feature.4.to_radians(), raw_feature.5.to_radians(),
            raw_feature.6.to_radians(), raw_feature.7, raw_feature.8
        ];
        let hash: GeometricHash = GeometricHash::TertiaryInteraction(
            HashValue::perfect_hash_default(raw_feature)
        );
        match hash {
            GeometricHash::TertiaryInteraction(hash) => {
                println!("{:?}", hash);
            },
            _ => panic!("Invalid hash type"),
        }
    }
}