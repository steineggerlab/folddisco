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

impl HashValue {
    pub fn perfect_hash(feature: Vec<f32>, nbin_dist: usize, nbin_angle: usize) -> Self {
        // Added one more quantization for distance
        let nbin_dist = if nbin_dist > 16 { 16.0 } else { nbin_dist as f32 };
        let nbin_angle = if nbin_angle > 4 { 4.0 } else { nbin_angle as f32 };
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
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

        let hashvalue = res1 << 25 | res2 << 20 | ca_dist << 16 
            | cb_dist << 12 | sin_ca_cb_angle << 10 | cos_ca_cb_angle << 8
            | sin_phi1 << 6 | cos_phi1 << 4 | sin_phi2 << 2 | cos_phi2;
        HashValue(hashvalue)
    }

    pub fn perfect_hash_default(feature: Vec<f32>) -> Self {
       Self::perfect_hash(feature, NBIN_DIST as usize, NBIN_SIN_COS as usize)
    }
    
    pub fn reverse_hash_default(&self) -> Vec<f32> {
        self.reverse_hash(NBIN_DIST as usize, NBIN_SIN_COS as usize)
    }
    
    pub fn reverse_hash(&self, nbin_dist: usize, nbin_angle: usize) -> Vec<f32> {
        let res1 = ((self.0 >> 25) & BITMASK32_5BIT)as f32;
        let res2 = ((self.0 >> 20) & BITMASK32_5BIT) as f32;
        let ca_dist = continuize_value(
            (self.0 >> 16) & BITMASK32_4BIT as u32, 
            MIN_DIST, MAX_DIST, nbin_dist as f32
        );
        let cb_dist = continuize_value(
            (self.0 >> 12) & BITMASK32_4BIT as u32,
            MIN_DIST, MAX_DIST, nbin_dist as f32
        );
        let sin_ca_cb_angle = continuize_value(
            (self.0 >> 10) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_ca_cb_angle = continuize_value(
            (self.0 >> 8) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let sin_phi1 = continuize_value(
            (self.0 >> 6) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_phi1 = continuize_value(
            (self.0 >> 4) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let sin_phi2 = continuize_value(
            (self.0 >> 2) & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_phi2 = continuize_value(
            self.0 & BITMASK32_2BIT as u32,
            MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );

        let ca_cb_angle = sin_ca_cb_angle.atan2(cos_ca_cb_angle).to_degrees();
        let phi1 = sin_phi1.atan2(cos_phi1).to_degrees();
        let phi2 = sin_phi2.atan2(cos_phi2).to_degrees();
        
        vec![res1, res2, ca_dist, cb_dist, ca_cb_angle, phi1, phi2]
    }
    
    pub fn hash_type(&self) -> HashType {
        HashType::PDBTrRosetta
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
            b"PHE", b"VAL", 14.0_f32, 15.9_f32, 116.0_f32, 80.0_f32, -100.0_f32
        );
        let raw_feature = vec![
            map_aa_to_u8(raw_feature.0) as f32, map_aa_to_u8(raw_feature.1) as f32,
            raw_feature.2, raw_feature.3, raw_feature.4.to_radians(),
            raw_feature.5.to_radians(), raw_feature.6.to_radians()
        ];
        let hash: GeometricHash = GeometricHash::PDBTrRosetta(
            HashValue::perfect_hash(raw_feature, 8, 3)
        );
        match hash {
            GeometricHash::PDBTrRosetta(hash) => {
                println!("{:?}", hash);
            },
            _ => panic!("Invalid hash type"),
        }
    }
}