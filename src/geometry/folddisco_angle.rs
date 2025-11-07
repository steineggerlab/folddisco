// File: folddisco_angle.rs
// Created: 2025-10-28 20:42:01
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2025 Hyunbin Kim, All rights reserved

use std::fmt;
use crate::geometry::core::HashType;
use crate::utils::convert::discretize_f32_value_into_u32 as discretize_value;
use crate::utils::convert::continuize_u32_value_into_f32 as continuize_value;
use crate::utils::convert::*;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(pub u32);

pub const NBIN_DIST: f32 = 8.0; // Limited to 8 because 3 bits are assigned
pub const NBIN_ANGLE_180: f32 = 32.0; // 5 bits
pub const NBIN_ANGLE_360: f32 = 32.0; // 5 bits

// 
pub const MIN_ANGLE_RAD: f32 = -std::f32::consts::PI;
pub const MAX_ANGLE_RAD: f32 = std::f32::consts::PI;


impl HashValue {
    #[inline]
    pub fn perfect_hash(feature: &Vec<f32>, nbin_dist: usize, nbin_angle: usize) -> u32 {
        let nbin_dist = if nbin_dist > NBIN_DIST as usize {
            NBIN_DIST
        } else if nbin_dist == 0 {
            NBIN_DIST
        } else {
            nbin_dist as f32
        };
        let nbin_angle = if nbin_angle > NBIN_ANGLE_360 as usize {
            NBIN_ANGLE_360
        } else if nbin_angle == 0 {
            NBIN_ANGLE_360
        } else {
            nbin_angle as f32
        };
        
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        let res_pair = map_aa_u32_pair_to_u32(res1, res2);

        let ca_dist = discretize_value(
            feature[2], MIN_DIST, MAX_DIST, nbin_dist
        );
        let cb_dist = discretize_value(
            feature[3], MIN_DIST, MAX_DIST, nbin_dist
        );
        // Angle is expected to be in radians but here just discretize radian directly
        
        // For ca-cb angle, max bin is 32 (5 bits)
        let ca_cb_angle = discretize_value(
            feature[4], 0.0, MAX_ANGLE_RAD, nbin_angle.min(NBIN_ANGLE_180)
        );

        // Two torsion angles: 
        let phi1 = discretize_value(
            feature[5], MIN_ANGLE_RAD, MAX_ANGLE_RAD, nbin_angle
        );
        let phi2 = discretize_value(
            feature[6], MIN_ANGLE_RAD, MAX_ANGLE_RAD, nbin_angle
        );

        // Bit map: 9 for residue pairs, 3 for ca_dist, 3 for cb_dist, 5 for ca-cb angle, 5 for phi1, 5 for phi2        
        let hashvalue = res_pair << 21 | ca_dist << 18 | cb_dist << 15 
            | ca_cb_angle << 10 | phi1 << 5 | phi2;
        hashvalue
    }
    
    pub fn perfect_hash_default(feature: &Vec<f32>) -> u32 {
        HashValue::perfect_hash(feature, NBIN_DIST as usize, NBIN_ANGLE_360 as usize)
    }
    
    pub fn reverse_hash_default(&self) -> [f32; 7] {
        self.reverse_hash(NBIN_DIST as usize, NBIN_ANGLE_360 as usize)
    }
    
    pub fn reverse_hash(&self, nbin_dist: usize, nbin_angle: usize) -> [f32; 7] {
        let res_pair = ((self.0 >> 21) & BITMASK32_9BIT) as u32;
        let (res1, res2) = map_u32_to_aa_u32_pair(res_pair);

        let ca_dist = continuize_value(
            (self.0 >> 18) & BITMASK32_3BIT as u32, 
            MIN_DIST, MAX_DIST, nbin_dist as f32
        );
        let cb_dist = continuize_value(
            (self.0 >> 15) & BITMASK32_3BIT as u32,
            MIN_DIST, MAX_DIST, nbin_dist as f32
        );
    
        let ca_cb_angle = continuize_value(
            (self.0 >> 10) & BITMASK32_5BIT as u32,
            0.0, MAX_ANGLE_RAD, (nbin_angle as f32).min(NBIN_ANGLE_180)
        );
        let phi1 = continuize_value(
            (self.0 >> 5) & BITMASK32_5BIT as u32,
            MIN_ANGLE_RAD, MAX_ANGLE_RAD, nbin_angle as f32
        );
        let phi2 = continuize_value(
            self.0 & BITMASK32_5BIT as u32,
            MIN_ANGLE_RAD, MAX_ANGLE_RAD, nbin_angle as f32
        );
        let ca_cb_angle = ca_cb_angle.to_degrees();
        let phi1 = phi1.to_degrees();
        let phi2 = phi2.to_degrees();

        [res1 as f32, res2 as f32, ca_dist, cb_dist, ca_cb_angle, phi1, phi2]
    }
    
    pub fn hash_type(&self) -> HashType {
        HashType::FolddiscoAngle
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
    use crate::utils::convert::map_aa_to_u8;
    
    #[test]
    fn test_geometrichash_works() {
        // Test perfect hash
        let raw_feature = (
            b"PHE", b"VAL", 14.0_f32, 15.9_f32, 116.0_f32, 80.0_f32, -100.0_f32
        );
        let raw_feature2 = (
            b"VAL", b"PHE", 15.9_f32, 14.0_f32, 116.0_f32, 80.0_f32, -100.0_f32
        );
        let raw_feature = vec![
            map_aa_to_u8(raw_feature.0) as f32, map_aa_to_u8(raw_feature.1) as f32,
            raw_feature.2, raw_feature.3, raw_feature.4.to_radians(),
            raw_feature.5.to_radians(), raw_feature.6.to_radians()
        ];
        let raw_feature2 = vec![
            map_aa_to_u8(raw_feature2.0) as f32, map_aa_to_u8(raw_feature2.1) as f32,
            raw_feature2.2, raw_feature2.3, raw_feature2.4.to_radians(),
            raw_feature2.5.to_radians(), raw_feature2.6.to_radians()
        ];
        let start = std::time::Instant::now();
        for _ in 0..10000 {
            let _ = HashValue::perfect_hash_default(&raw_feature);
            let _ = HashValue::perfect_hash_default(&raw_feature2);
        }
        let duration = start.elapsed();
        println!("Time elapsed in perfect_hash_default() is: {:?}", duration);
        let hash = HashValue::perfect_hash_default(&raw_feature);
        let hash = GeometricHash::FolddiscoAngle(HashValue::from_u32(hash));
        match hash {
            GeometricHash::PDBTrRosetta(hash) => {
                println!("{:?}", hash);
            },
            _ => panic!("Invalid hash type"),
        }
    }

}