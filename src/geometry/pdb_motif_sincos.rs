// File: pdb_motif_sincos.rs
// Created: 2024-01-18 15:48:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use std::fmt;
use crate::geometry::core::HashType;
use crate::utils::convert::discretize_f32_value_into_u32 as discretize_value;
use crate::utils::convert::continuize_u32_value_into_f32 as continuize_value;
use crate::utils::convert::*;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(pub u32);

impl HashValue {
    pub fn perfect_hash(feature: &Vec<f32>, nbin_dist: usize, nbin_angle: usize) -> u32 {
        // Added one more quantization for distance
        let nbin_dist = if nbin_dist > 16 { 
            16.0
        } else if nbin_dist == 0 {
            NBIN_DIST
        } else { 
            nbin_dist as f32 
        };
        let nbin_angle = if nbin_angle > 16 {
            16.0
        } else if nbin_angle == 0 {
            NBIN_SIN_COS
        } else {
            nbin_angle as f32
        };
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        let ca_dist = discretize_value(
            feature[2], MIN_DIST, MAX_DIST, nbin_dist
        );
        let cb_dist = discretize_value(
            feature[3], MIN_DIST, MAX_DIST, nbin_dist
        );
        // Angle is expected to be in radians
        let sin_angle = feature[4].sin();
        let cos_angle = feature[4].cos();
        let sin_angle = discretize_value(
            sin_angle, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let cos_angle = discretize_value(
            cos_angle, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
        );
        let hashvalue = res1 << 21 | res2 << 16 | ca_dist << 12 
            | cb_dist << 8 | sin_angle << 4 | cos_angle;
        hashvalue
    }
    
    pub fn perfect_hash_default(feature: &Vec<f32>) -> u32 {
        HashValue::perfect_hash(feature, NBIN_DIST as usize, NBIN_SIN_COS as usize)
    }
    
    pub fn reverse_hash(&self, nbin_dist: usize, nbin_angle: usize) -> [f32; 5] {
        let res1 = ((self.0 >> 21) & BITMASK32_5BIT)as f32;
        let res2 = ((self.0 >> 16) & BITMASK32_5BIT) as f32;
        let ca_dist = continuize_value(
            (self.0 >> 12) & BITMASK32_4BIT as u32, 
            MIN_DIST, MAX_DIST, nbin_dist as f32
        );
        let cb_dist = continuize_value(
            (self.0 >> 8) & BITMASK32_4BIT as u32,
            MIN_DIST, MAX_DIST, nbin_dist as f32
        );
        let sin_angle = continuize_value(
            (self.0 >> 4) & BITMASK32_4BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_angle = continuize_value(
            self.0 & BITMASK32_4BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let angle = sin_angle.atan2(cos_angle).to_degrees();
        [res1, res2, ca_dist, cb_dist, angle]
    }
    
    pub fn reverse_hash_default(&self) -> [f32; 5] {
        self.reverse_hash(NBIN_DIST as usize, NBIN_SIN_COS as usize)
    }
    
    pub fn hash_type(&self) -> HashType {
        HashType::PDBMotifSinCos
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
    use crate::geometry::core::GeometricHash;
    use crate::utils::convert::map_aa_to_u8;
    #[test]
    fn test_geometrichash_works() {
        // Test perfect hash
        let raw_feature = (
            b"PHE", b"VAL", 14.0_f32, 15.9_f32, 116.0_f32
        );
        let raw_feature = vec![
            map_aa_to_u8(raw_feature.0) as f32, map_aa_to_u8(raw_feature.1) as f32,
            raw_feature.2, raw_feature.3, raw_feature.4.to_radians()
        ];
        let hash: GeometricHash = GeometricHash::PDBMotifSinCos(
            HashValue(HashValue::perfect_hash(&raw_feature, 8, 3))
        );
        match hash {
            GeometricHash::PDBMotifSinCos(hash) => {
                println!("{:?}", hash);
            },
            _ => panic!("Invalid hash type"),
        }
    }
}