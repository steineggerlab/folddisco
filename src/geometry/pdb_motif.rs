// File: pdb_motif.rs
// Created: 2024-01-18 15:48:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use std::fmt;
use crate::geometry::core::HashType;
use crate::geometry::util::discretize_f32_value_into_u32 as discretize_value;
use crate::geometry::util::continuize_u32_value_into_f32 as continuize_value;

// Residue 1: 5 bits; Residue 2: 5 bits; Distances: 16 bins 4 bits; 
// Angle: 32 bins 5 bits; total: 23 bits
const MIN_DIST: f32 = 2.0;
const MAX_DIST: f32 = 20.0;
const NBIN_DIST: f32 = 16.0;
const MIN_ANGLE: f32 = 0.0;
const MAX_ANGLE: f32 = 180.0;
const NBIN_ANGLE: f32 = 16.0;
// Bitmasks
const BITMASK32_4BIT: u32 = 0x0000000F;
const BITMASK32_5BIT: u32 = 0x0000001F;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(pub u32);

impl HashValue {
    pub fn perfect_hash(feature: Vec<f32>) -> Self {
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        let ca_dist = discretize_value(feature[2], MIN_DIST, MAX_DIST, NBIN_DIST);
        let cb_dist = discretize_value(feature[3], MIN_DIST, MAX_DIST, NBIN_DIST);
        let angle = discretize_value(feature[4], MIN_ANGLE, MAX_ANGLE, NBIN_ANGLE);
        let hashvalue = res1 << 18 | res2 << 13 | ca_dist << 9 | cb_dist << 5 | angle;
        HashValue(hashvalue)
    }
    pub fn reverse_hash(&self) -> Vec<f32> {
        let res1 = ((self.0 >> 18) & BITMASK32_5BIT)as f32;
        let res2 = ((self.0 >> 13) & BITMASK32_5BIT) as f32;
        let ca_dist = continuize_value(
            (self.0 >> 9) & BITMASK32_4BIT as u32, 
            MIN_DIST, MAX_DIST, NBIN_DIST
        );
        let cb_dist = continuize_value(
            (self.0 >> 5) & BITMASK32_4BIT as u32, 
            MIN_DIST, MAX_DIST, NBIN_DIST
        );
        let angle = continuize_value(
            self.0 & BITMASK32_5BIT as u32, 
            MIN_ANGLE, MAX_ANGLE, NBIN_ANGLE
        );
        vec![res1, res2, ca_dist, cb_dist, angle]
    }
    pub fn hash_type(&self) -> HashType {
        HashType::PDBMotif
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
        let values = self.reverse_hash();
        write!(f, "HashValue({}), values={:?}", self.0, values)
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let values = self.reverse_hash();
        write!(f, "{}\t{:?}", self.0, values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::util::map_aa_to_u8;
    #[test]
    fn test_hash_works() {
        // Test perfect hash
        let raw_feature = (
            b"PHE", b"VAL", 14.0_f32, 15.9_f32, 116.0_f32
        );
        let raw_feature = vec![
            map_aa_to_u8(raw_feature.0) as f32, map_aa_to_u8(raw_feature.1) as f32,
            raw_feature.2, raw_feature.3, raw_feature.4
        ];

        let hash = HashValue::perfect_hash(raw_feature.clone());
        let reversed = hash.reverse_hash();
        println!("{:?}", hash);
    }
}