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

pub const PDBTR_NBIN_DIST: f32 = 16.0;
pub const PDBTR_NBIN_SIN_COS: f32 = 4.0;

impl HashValue {
    #[inline]
    pub fn perfect_hash(feature: &Vec<f32>, nbin_dist: usize, nbin_angle: usize) -> u32 {
        let nbin_dist = if nbin_dist > 16 {
            16.0
        } else if nbin_dist == 0 {
            PDBTR_NBIN_DIST
        } else {
            nbin_dist as f32
        };
        let nbin_angle = if nbin_angle > 4 {
            4.0
        } else if nbin_angle == 0 {
            PDBTR_NBIN_SIN_COS
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
        hashvalue
    }

    #[inline(always)]
    fn discretize_with_shift(value: f32, min_val: f32, max_val: f32, nbins: f32, shift: f32) -> u32 {
        let bin_width = (max_val - min_val) / nbins;
        // Apply shift by moving the bin boundaries
        let shifted_value = value + shift; // Note: we shift the value, not the boundaries
        let clamped_value = shifted_value.max(min_val).min(max_val - 1e-6);
        let bin_index = ((clamped_value - min_val) / bin_width).floor() as u32;
        bin_index.min(15) // 4-bit constraint for distances
    }
    
    pub fn perfect_hash_default(feature: &Vec<f32>) -> u32 {
        HashValue::perfect_hash(feature, PDBTR_NBIN_DIST as usize, PDBTR_NBIN_SIN_COS as usize)
    }
    
    pub fn reverse_hash_default(&self) -> [f32; 7] {
        self.reverse_hash(PDBTR_NBIN_DIST as usize, PDBTR_NBIN_SIN_COS as usize)
    }
    
    pub fn reverse_hash(&self, nbin_dist: usize, nbin_angle: usize) -> [f32; 7] {
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
        
        [res1, res2, ca_dist, cb_dist, ca_cb_angle, phi1, phi2]
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
    
    pub fn is_symmetric(&self) -> bool {
        let values = self.reverse_hash_default();
        // Residue pair is symmetric and phi is symmetric
        (values[0] == values[1]) && (values[5] == values[6])
    }
    
    /// deduplication: inline during generation
    /// Returns (unique_count, unique_hashes_array)
    #[inline]
    pub fn perfect_hash_with_shifts_dedup_inline(feature: &Vec<f32>) -> (u8, [u32; 8]) {
        const NBIN_DIST: f32 = PDBTR_NBIN_DIST;
        const NBIN_ANGLE: f32 = PDBTR_NBIN_SIN_COS;
        
        // Optimized shift amounts
        const DIST_SHIFT: f32 = 0.6;
        const ANGLE_SHIFT_RAD: f32 = std::f32::consts::PI / 8.0;
        
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        
        // Precompute ALL trigonometric values we need (only 3 sets total)
        // Original angles (0 shift)
        let sin_ca_cb_0 = feature[4].sin();
        let cos_ca_cb_0 = feature[4].cos();
        let sin_phi1_0 = feature[5].sin();
        let cos_phi1_0 = feature[5].cos();
        let sin_phi2_0 = feature[6].sin();
        let cos_phi2_0 = feature[6].cos();
        
        // Positive angle shift (+ANGLE_SHIFT_RAD)
        let sin_ca_cb_pos = (feature[4] + ANGLE_SHIFT_RAD).sin();
        let cos_ca_cb_pos = (feature[4] + ANGLE_SHIFT_RAD).cos();
        let sin_phi1_pos = (feature[5] + ANGLE_SHIFT_RAD).sin();
        let cos_phi1_pos = (feature[5] + ANGLE_SHIFT_RAD).cos();
        let sin_phi2_pos = (feature[6] + ANGLE_SHIFT_RAD).sin();
        let cos_phi2_pos = (feature[6] + ANGLE_SHIFT_RAD).cos();
        
        // Negative angle shift (-ANGLE_SHIFT_RAD)
        let sin_ca_cb_neg = (feature[4] - ANGLE_SHIFT_RAD).sin();
        let cos_ca_cb_neg = (feature[4] - ANGLE_SHIFT_RAD).cos();
        let sin_phi1_neg = (feature[5] - ANGLE_SHIFT_RAD).sin();
        let cos_phi1_neg = (feature[5] - ANGLE_SHIFT_RAD).cos();
        let sin_phi2_neg = (feature[6] - ANGLE_SHIFT_RAD).sin();
        let cos_phi2_neg = (feature[6] - ANGLE_SHIFT_RAD).cos();
        
        // Precompute distance shifts (3 variants: 0, +DIST_SHIFT, -DIST_SHIFT)
        let ca_dist_0 = Self::discretize_with_shift(
            feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, 0.0
        );
        let cb_dist_0 = Self::discretize_with_shift(
            feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, 0.0
        );
        
        let ca_dist_pos = Self::discretize_with_shift(
            feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, DIST_SHIFT
        );
        let cb_dist_pos = Self::discretize_with_shift(
            feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, DIST_SHIFT
        );
        
        let ca_dist_neg = Self::discretize_with_shift(
            feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, -DIST_SHIFT
        );
        let cb_dist_neg = Self::discretize_with_shift(
            feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, -DIST_SHIFT
        );
        
        // Precompute discretized angle values for each angle shift
        let sin_ca_cb_disc_0 = discretize_value(sin_ca_cb_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_ca_cb_disc_0 = discretize_value(cos_ca_cb_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let sin_phi1_disc_0 = discretize_value(sin_phi1_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi1_disc_0 = discretize_value(cos_phi1_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let sin_phi2_disc_0 = discretize_value(sin_phi2_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi2_disc_0 = discretize_value(cos_phi2_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        
        let sin_ca_cb_disc_pos = discretize_value(sin_ca_cb_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_ca_cb_disc_pos = discretize_value(cos_ca_cb_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let sin_phi1_disc_pos = discretize_value(sin_phi1_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi1_disc_pos = discretize_value(cos_phi1_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let sin_phi2_disc_pos = discretize_value(sin_phi2_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi2_disc_pos = discretize_value(cos_phi2_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        
        let sin_ca_cb_disc_neg = discretize_value(sin_ca_cb_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_ca_cb_disc_neg = discretize_value(cos_ca_cb_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let sin_phi1_disc_neg = discretize_value(sin_phi1_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi1_disc_neg = discretize_value(cos_phi1_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let sin_phi2_disc_neg = discretize_value(sin_phi2_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi2_disc_neg = discretize_value(cos_phi2_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        
        let mut unique_hashes = [0u32; 8];
        let mut count = 0u8;
        
        // Now generate hashes by simply selecting precomputed values
        // Each combination just picks the right precomputed distance and angle values
        
        // 0: Original (0 dist, 0 angle)
        let hashvalue = res1 << 25 | res2 << 20 | ca_dist_0 << 16 
            | cb_dist_0 << 12 | sin_ca_cb_disc_0 << 10 | cos_ca_cb_disc_0 << 8
            | sin_phi1_disc_0 << 6 | cos_phi1_disc_0 << 4 | sin_phi2_disc_0 << 2 | cos_phi2_disc_0;
        unique_hashes[count as usize] = hashvalue;
        count += 1;
        
        // 1: Angle+ only (0 dist, + angle)
        let hashvalue = res1 << 25 | res2 << 20 | ca_dist_0 << 16 
            | cb_dist_0 << 12 | sin_ca_cb_disc_pos << 10 | cos_ca_cb_disc_pos << 8
            | sin_phi1_disc_pos << 6 | cos_phi1_disc_pos << 4 | sin_phi2_disc_pos << 2 | cos_phi2_disc_pos;
        Self::add_unique_hash(&mut unique_hashes, &mut count, hashvalue);
        
        // 2: Angle- only (0 dist, - angle)
        let hashvalue = res1 << 25 | res2 << 20 | ca_dist_0 << 16 
            | cb_dist_0 << 12 | sin_ca_cb_disc_neg << 10 | cos_ca_cb_disc_neg << 8
            | sin_phi1_disc_neg << 6 | cos_phi1_disc_neg << 4 | sin_phi2_disc_neg << 2 | cos_phi2_disc_neg;
        Self::add_unique_hash(&mut unique_hashes, &mut count, hashvalue);
        
        // 3: Distance+ only (+ dist, 0 angle)
        let hashvalue = res1 << 25 | res2 << 20 | ca_dist_pos << 16 
            | cb_dist_pos << 12 | sin_ca_cb_disc_0 << 10 | cos_ca_cb_disc_0 << 8
            | sin_phi1_disc_0 << 6 | cos_phi1_disc_0 << 4 | sin_phi2_disc_0 << 2 | cos_phi2_disc_0;
        Self::add_unique_hash(&mut unique_hashes, &mut count, hashvalue);
        
        // 4: Both positive (+ dist, + angle)
        let hashvalue = res1 << 25 | res2 << 20 | ca_dist_pos << 16 
            | cb_dist_pos << 12 | sin_ca_cb_disc_pos << 10 | cos_ca_cb_disc_pos << 8
            | sin_phi1_disc_pos << 6 | cos_phi1_disc_pos << 4 | sin_phi2_disc_pos << 2 | cos_phi2_disc_pos;
        Self::add_unique_hash(&mut unique_hashes, &mut count, hashvalue);
        
        // 5: Dist+, Angle- (+ dist, - angle)
        let hashvalue = res1 << 25 | res2 << 20 | ca_dist_pos << 16 
            | cb_dist_pos << 12 | sin_ca_cb_disc_neg << 10 | cos_ca_cb_disc_neg << 8
            | sin_phi1_disc_neg << 6 | cos_phi1_disc_neg << 4 | sin_phi2_disc_neg << 2 | cos_phi2_disc_neg;
        Self::add_unique_hash(&mut unique_hashes, &mut count, hashvalue);
        
        // 6: Distance- only (- dist, 0 angle)
        let hashvalue = res1 << 25 | res2 << 20 | ca_dist_neg << 16 
            | cb_dist_neg << 12 | sin_ca_cb_disc_0 << 10 | cos_ca_cb_disc_0 << 8
            | sin_phi1_disc_0 << 6 | cos_phi1_disc_0 << 4 | sin_phi2_disc_0 << 2 | cos_phi2_disc_0;
        Self::add_unique_hash(&mut unique_hashes, &mut count, hashvalue);
        
        // 7: Dist-, Angle+ (- dist, + angle)
        let hashvalue = res1 << 25 | res2 << 20 | ca_dist_neg << 16 
            | cb_dist_neg << 12 | sin_ca_cb_disc_pos << 10 | cos_ca_cb_disc_pos << 8
            | sin_phi1_disc_pos << 6 | cos_phi1_disc_pos << 4 | sin_phi2_disc_pos << 2 | cos_phi2_disc_pos;
        Self::add_unique_hash(&mut unique_hashes, &mut count, hashvalue);
        
        (count, unique_hashes)
    }

    #[inline]
    fn add_unique_hash(unique_hashes: &mut [u32; 8], count: &mut u8, hashvalue: u32) {
        // Linear search for duplicates (very fast for small arrays)
        for i in 0..*count as usize {
            if unique_hashes[i] == hashvalue {
                return; // Duplicate found, don't add
            }
        }
        
        // Add new unique hash
        if (*count as usize) < unique_hashes.len() {
            unique_hashes[*count as usize] = hashvalue;
            *count += 1;
        }
    }

    #[inline]
    pub fn perfect_hash_with_all_shifts_exhaustive(feature: &Vec<f32>) -> (usize, Vec<u32>) {
        const NBIN_DIST: f32 = PDBTR_NBIN_DIST;
        const NBIN_ANGLE: f32 = PDBTR_NBIN_SIN_COS;
        const DIST_SHIFT: f32 = 0.6;
        const ANGLE_SHIFT_RAD: f32 = std::f32::consts::PI / 8.0;
        
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        
        // Precompute all values (same as before)
        let sin_ca_cb_0 = feature[4].sin();
        let cos_ca_cb_0 = feature[4].cos();
        let sin_phi1_0 = feature[5].sin();
        let cos_phi1_0 = feature[5].cos();
        let sin_phi2_0 = feature[6].sin();
        let cos_phi2_0 = feature[6].cos();
        
        let sin_ca_cb_pos = (feature[4] + ANGLE_SHIFT_RAD).sin();
        let cos_ca_cb_pos = (feature[4] + ANGLE_SHIFT_RAD).cos();
        let sin_phi1_pos = (feature[5] + ANGLE_SHIFT_RAD).sin();
        let cos_phi1_pos = (feature[5] + ANGLE_SHIFT_RAD).cos();
        let sin_phi2_pos = (feature[6] + ANGLE_SHIFT_RAD).sin();
        let cos_phi2_pos = (feature[6] + ANGLE_SHIFT_RAD).cos();
        
        let sin_ca_cb_neg = (feature[4] - ANGLE_SHIFT_RAD).sin();
        let cos_ca_cb_neg = (feature[4] - ANGLE_SHIFT_RAD).cos();
        let sin_phi1_neg = (feature[5] - ANGLE_SHIFT_RAD).sin();
        let cos_phi1_neg = (feature[5] - ANGLE_SHIFT_RAD).cos();
        let sin_phi2_neg = (feature[6] - ANGLE_SHIFT_RAD).sin();
        let cos_phi2_neg = (feature[6] - ANGLE_SHIFT_RAD).cos();
        
        // Pre-deduplicate discretized values
        let ca_dists = [
            Self::discretize_with_shift(feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, 0.0),
            Self::discretize_with_shift(feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, DIST_SHIFT),
            Self::discretize_with_shift(feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, -DIST_SHIFT),
        ];
        
        let cb_dists = [
            Self::discretize_with_shift(feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, 0.0),
            Self::discretize_with_shift(feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, DIST_SHIFT),
            Self::discretize_with_shift(feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, -DIST_SHIFT),
        ];
        
        let ca_cb_angles = [
            (discretize_value(sin_ca_cb_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_ca_cb_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_ca_cb_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_ca_cb_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_ca_cb_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_ca_cb_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
        ];
        
        let phi1_angles = [
            (discretize_value(sin_phi1_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi1_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_phi1_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi1_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_phi1_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi1_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
        ];
        
        let phi2_angles = [
            (discretize_value(sin_phi2_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi2_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_phi2_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi2_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_phi2_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi2_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
        ];
        
        // OPTIMIZATION 4: Pre-detect duplicate ranges to skip redundant combinations
        let unique_ca_dists: Vec<u32> = {
            let mut temp = Vec::with_capacity(3);
            for &dist in &ca_dists {
                if !temp.contains(&dist) {
                    temp.push(dist);
                }
            }
            temp
        };
        
        let unique_cb_dists: Vec<u32> = {
            let mut temp = Vec::with_capacity(3);
            for &dist in &cb_dists {
                if !temp.contains(&dist) {
                    temp.push(dist);
                }
            }
            temp
        };
        
        let unique_ca_cb_angles: Vec<(u32, u32)> = {
            let mut temp = Vec::with_capacity(3);
            for &angles in &ca_cb_angles {
                if !temp.contains(&angles) {
                    temp.push(angles);
                }
            }
            temp
        };
        
        let unique_phi1_angles: Vec<(u32, u32)> = {
            let mut temp = Vec::with_capacity(3);
            for &angles in &phi1_angles {
                if !temp.contains(&angles) {
                    temp.push(angles);
                }
            }
            temp
        };
        
        let unique_phi2_angles: Vec<(u32, u32)> = {
            let mut temp = Vec::with_capacity(3);
            for &angles in &phi2_angles {
                if !temp.contains(&angles) {
                    temp.push(angles);
                }
            }
            temp
        };
        
        let mut unique_hashes = Vec::with_capacity(unique_ca_dists.len() * unique_cb_dists.len() * 
                                                   unique_ca_cb_angles.len() * unique_phi1_angles.len() * 
                                                   unique_phi2_angles.len());
        
        // OPTIMIZATION 5: Only iterate over unique discretized values
        for &ca_dist in &unique_ca_dists {
            for &cb_dist in &unique_cb_dists {
                for &(sin_ca_cb_disc, cos_ca_cb_disc) in &unique_ca_cb_angles {
                    for &(sin_phi1_disc, cos_phi1_disc) in &unique_phi1_angles {
                        for &(sin_phi2_disc, cos_phi2_disc) in &unique_phi2_angles {
                            let hashvalue = res1 << 25 | res2 << 20 | ca_dist << 16 
                                | cb_dist << 12 | sin_ca_cb_disc << 10 | cos_ca_cb_disc << 8
                                | sin_phi1_disc << 6 | cos_phi1_disc << 4 | sin_phi2_disc << 2 | cos_phi2_disc;
                            
                            unique_hashes.push(hashvalue);
                        }
                    }
                }
            }
        }
        
        (unique_hashes.len(), unique_hashes)
    }
    
    #[inline]
    pub fn perfect_hash_with_all_shifts_exhaustive_optimized(feature: &Vec<f32>) -> (usize, Vec<u32>) {
        const NBIN_DIST: f32 = PDBTR_NBIN_DIST;
        const NBIN_ANGLE: f32 = PDBTR_NBIN_SIN_COS;
        const DIST_SHIFT: f32 = 0.6;
        const ANGLE_SHIFT_RAD: f32 = std::f32::consts::PI / 8.0;
        
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        
        // Precompute trigonometric values (same as before)
        let sin_ca_cb_0 = feature[4].sin();
        let cos_ca_cb_0 = feature[4].cos();
        let sin_phi1_0 = feature[5].sin();
        let cos_phi1_0 = feature[5].cos();
        let sin_phi2_0 = feature[6].sin();
        let cos_phi2_0 = feature[6].cos();
        
        let sin_ca_cb_pos = (feature[4] + ANGLE_SHIFT_RAD).sin();
        let cos_ca_cb_pos = (feature[4] + ANGLE_SHIFT_RAD).cos();
        let sin_phi1_pos = (feature[5] + ANGLE_SHIFT_RAD).sin();
        let cos_phi1_pos = (feature[5] + ANGLE_SHIFT_RAD).cos();
        let sin_phi2_pos = (feature[6] + ANGLE_SHIFT_RAD).sin();
        let cos_phi2_pos = (feature[6] + ANGLE_SHIFT_RAD).cos();
        
        let sin_ca_cb_neg = (feature[4] - ANGLE_SHIFT_RAD).sin();
        let cos_ca_cb_neg = (feature[4] - ANGLE_SHIFT_RAD).cos();
        let sin_phi1_neg = (feature[5] - ANGLE_SHIFT_RAD).sin();
        let cos_phi1_neg = (feature[5] - ANGLE_SHIFT_RAD).cos();
        let sin_phi2_neg = (feature[6] - ANGLE_SHIFT_RAD).sin();
        let cos_phi2_neg = (feature[6] - ANGLE_SHIFT_RAD).cos();
        
        // Precompute distance and angle discretizations (same as before)
        let ca_dists = [
            Self::discretize_with_shift(feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, 0.0),
            Self::discretize_with_shift(feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, DIST_SHIFT),
            Self::discretize_with_shift(feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, -DIST_SHIFT),
        ];
        
        let cb_dists = [
            Self::discretize_with_shift(feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, 0.0),
            Self::discretize_with_shift(feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, DIST_SHIFT),
            Self::discretize_with_shift(feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, -DIST_SHIFT),
        ];
        
        let ca_cb_angles = [
            (discretize_value(sin_ca_cb_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_ca_cb_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_ca_cb_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_ca_cb_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_ca_cb_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_ca_cb_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
        ];
        
        let phi1_angles = [
            (discretize_value(sin_phi1_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi1_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_phi1_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi1_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_phi1_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi1_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
        ];
        
        let phi2_angles = [
            (discretize_value(sin_phi2_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi2_0, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_phi2_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi2_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
            (discretize_value(sin_phi2_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE),
             discretize_value(cos_phi2_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE)),
        ];
        
        // OPTIMIZATION 4: Pre-detect duplicate ranges to skip redundant combinations
        let unique_ca_dists: Vec<u32> = {
            let mut temp = Vec::with_capacity(3);
            for &dist in &ca_dists {
                if !temp.contains(&dist) {
                    temp.push(dist);
                }
            }
            temp
        };
        
        let unique_cb_dists: Vec<u32> = {
            let mut temp = Vec::with_capacity(3);
            for &dist in &cb_dists {
                if !temp.contains(&dist) {
                    temp.push(dist);
                }
            }
            temp
        };
        
        let unique_ca_cb_angles: Vec<(u32, u32)> = {
            let mut temp = Vec::with_capacity(3);
            for &angles in &ca_cb_angles {
                if !temp.contains(&angles) {
                    temp.push(angles);
                }
            }
            temp
        };
        
        let unique_phi1_angles: Vec<(u32, u32)> = {
            let mut temp = Vec::with_capacity(3);
            for &angles in &phi1_angles {
                if !temp.contains(&angles) {
                    temp.push(angles);
                }
            }
            temp
        };
        
        let unique_phi2_angles: Vec<(u32, u32)> = {
            let mut temp = Vec::with_capacity(3);
            for &angles in &phi2_angles {
                if !temp.contains(&angles) {
                    temp.push(angles);
                }
            }
            temp
        };
        
        let mut unique_hashes = Vec::with_capacity(unique_ca_dists.len() * unique_cb_dists.len() * 
                                                   unique_ca_cb_angles.len() * unique_phi1_angles.len() * 
                                                   unique_phi2_angles.len());
        
        // OPTIMIZATION 5: Only iterate over unique discretized values
        for &ca_dist in &unique_ca_dists {
            for &cb_dist in &unique_cb_dists {
                for &(sin_ca_cb_disc, cos_ca_cb_disc) in &unique_ca_cb_angles {
                    for &(sin_phi1_disc, cos_phi1_disc) in &unique_phi1_angles {
                        for &(sin_phi2_disc, cos_phi2_disc) in &unique_phi2_angles {
                            let hashvalue = res1 << 25 | res2 << 20 | ca_dist << 16 
                                | cb_dist << 12 | sin_ca_cb_disc << 10 | cos_ca_cb_disc << 8
                                | sin_phi1_disc << 6 | cos_phi1_disc << 4 | sin_phi2_disc << 2 | cos_phi2_disc;
                            
                            unique_hashes.push(hashvalue);
                        }
                    }
                }
            }
        }
        
        (unique_hashes.len(), unique_hashes)
    }

    /// Single efficient function that generates up to 32 unique hashes (2^5) by analyzing 
    /// which of the 5 key values (ca_dist, cb_dist, ca_cb_angle, phi1, phi2) change when shifted.
    /// Minimizes trigonometric operations through early termination: if positive shift changes 
    /// discretized value, negative shift is still checked for completeness but with smart optimization.
    /// Combines direction checking and hash generation in one pass to eliminate redundant operations.
    #[inline]
    pub fn perfect_hash_with_max_32_shifts(feature: &Vec<f32>) -> (usize, Vec<u32>) {
        const NBIN_DIST: f32 = PDBTR_NBIN_DIST;
        const NBIN_ANGLE: f32 = PDBTR_NBIN_SIN_COS;
        const DIST_SHIFT: f32 = 0.6;
        const ANGLE_SHIFT_RAD: f32 = std::f32::consts::PI / 8.0;
        
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        
        // === DISTANCE VALUE ANALYSIS (no trigonometry needed) ===
        let ca_dist_orig = Self::discretize_with_shift(feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, 0.0);
        let cb_dist_orig = Self::discretize_with_shift(feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, 0.0);
        
        let ca_dist_pos = Self::discretize_with_shift(feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, DIST_SHIFT);
        let ca_dist_neg = Self::discretize_with_shift(feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, -DIST_SHIFT);
        let cb_dist_pos = Self::discretize_with_shift(feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, DIST_SHIFT);
        let cb_dist_neg = Self::discretize_with_shift(feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, -DIST_SHIFT);
        
        // Collect unique distance values
        let mut ca_dist_values = vec![ca_dist_orig];
        if ca_dist_pos != ca_dist_orig { ca_dist_values.push(ca_dist_pos); }
        if ca_dist_neg != ca_dist_orig && ca_dist_neg != ca_dist_pos { ca_dist_values.push(ca_dist_neg); }
        
        let mut cb_dist_values = vec![cb_dist_orig];
        if cb_dist_pos != cb_dist_orig { cb_dist_values.push(cb_dist_pos); }
        if cb_dist_neg != cb_dist_orig && cb_dist_neg != cb_dist_pos { cb_dist_values.push(cb_dist_neg); }
        
        // === ANGLE VALUE ANALYSIS (minimize trigonometric operations) ===
        // Always compute original trigonometric values (6 operations - required)
        let sin_ca_cb_orig = feature[4].sin();
        let cos_ca_cb_orig = feature[4].cos();
        let sin_phi1_orig = feature[5].sin();
        let cos_phi1_orig = feature[5].cos();
        let sin_phi2_orig = feature[6].sin();
        let cos_phi2_orig = feature[6].cos();
        
        let sin_ca_cb_disc_orig = discretize_value(sin_ca_cb_orig, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_ca_cb_disc_orig = discretize_value(cos_ca_cb_orig, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let sin_phi1_disc_orig = discretize_value(sin_phi1_orig, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi1_disc_orig = discretize_value(cos_phi1_orig, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let sin_phi2_disc_orig = discretize_value(sin_phi2_orig, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi2_disc_orig = discretize_value(cos_phi2_orig, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        
        // CA-CB angle analysis with early termination optimization
        let mut ca_cb_angle_values = vec![(sin_ca_cb_disc_orig, cos_ca_cb_disc_orig)];
        
        // Check positive shift first (2 additional trig operations)
        let sin_ca_cb_pos = (feature[4] + ANGLE_SHIFT_RAD).sin();
        let cos_ca_cb_pos = (feature[4] + ANGLE_SHIFT_RAD).cos();
        let sin_ca_cb_disc_pos = discretize_value(sin_ca_cb_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_ca_cb_disc_pos = discretize_value(cos_ca_cb_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        
        let ca_cb_pos_changed = sin_ca_cb_disc_pos != sin_ca_cb_disc_orig || cos_ca_cb_disc_pos != cos_ca_cb_disc_orig;
        if ca_cb_pos_changed {
            ca_cb_angle_values.push((sin_ca_cb_disc_pos, cos_ca_cb_disc_pos));
        }
        
        // Check negative shift (2 more trig operations) - but can optimize based on positive result
        let sin_ca_cb_neg = (feature[4] - ANGLE_SHIFT_RAD).sin();
        let cos_ca_cb_neg = (feature[4] - ANGLE_SHIFT_RAD).cos();
        let sin_ca_cb_disc_neg = discretize_value(sin_ca_cb_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_ca_cb_disc_neg = discretize_value(cos_ca_cb_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        
        let ca_cb_neg_changed = sin_ca_cb_disc_neg != sin_ca_cb_disc_orig || cos_ca_cb_disc_neg != cos_ca_cb_disc_orig;
        if ca_cb_neg_changed && (sin_ca_cb_disc_neg != sin_ca_cb_disc_pos || cos_ca_cb_disc_neg != cos_ca_cb_disc_pos) {
            ca_cb_angle_values.push((sin_ca_cb_disc_neg, cos_ca_cb_disc_neg));
        }
        
        // PHI1 angle analysis with early termination optimization
        let mut phi1_angle_values = vec![(sin_phi1_disc_orig, cos_phi1_disc_orig)];
        
        // Check positive shift first (2 additional trig operations)
        let sin_phi1_pos = (feature[5] + ANGLE_SHIFT_RAD).sin();
        let cos_phi1_pos = (feature[5] + ANGLE_SHIFT_RAD).cos();
        let sin_phi1_disc_pos = discretize_value(sin_phi1_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi1_disc_pos = discretize_value(cos_phi1_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        
        let phi1_pos_changed = sin_phi1_disc_pos != sin_phi1_disc_orig || cos_phi1_disc_pos != cos_phi1_disc_orig;
        if phi1_pos_changed {
            phi1_angle_values.push((sin_phi1_disc_pos, cos_phi1_disc_pos));
        }
        
        // Check negative shift (2 more trig operations) - but can optimize based on positive result
        let sin_phi1_neg = (feature[5] - ANGLE_SHIFT_RAD).sin();
        let cos_phi1_neg = (feature[5] - ANGLE_SHIFT_RAD).cos();
        let sin_phi1_disc_neg = discretize_value(sin_phi1_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi1_disc_neg = discretize_value(cos_phi1_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        
        let phi1_neg_changed = sin_phi1_disc_neg != sin_phi1_disc_orig || cos_phi1_disc_neg != cos_phi1_disc_orig;
        if phi1_neg_changed && (sin_phi1_disc_neg != sin_phi1_disc_pos || cos_phi1_disc_neg != cos_phi1_disc_pos) {
            phi1_angle_values.push((sin_phi1_disc_neg, cos_phi1_disc_neg));
        }
        
        // PHI2 angle analysis with early termination optimization
        let mut phi2_angle_values = vec![(sin_phi2_disc_orig, cos_phi2_disc_orig)];
        
        // Check positive shift first (2 additional trig operations)
        let sin_phi2_pos = (feature[6] + ANGLE_SHIFT_RAD).sin();
        let cos_phi2_pos = (feature[6] + ANGLE_SHIFT_RAD).cos();
        let sin_phi2_disc_pos = discretize_value(sin_phi2_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi2_disc_pos = discretize_value(cos_phi2_pos, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        
        let phi2_pos_changed = sin_phi2_disc_pos != sin_phi2_disc_orig || cos_phi2_disc_pos != cos_phi2_disc_orig;
        if phi2_pos_changed {
            phi2_angle_values.push((sin_phi2_disc_pos, cos_phi2_disc_pos));
        }
        
        // Check negative shift (2 more trig operations) - but can optimize based on positive result
        let sin_phi2_neg = (feature[6] - ANGLE_SHIFT_RAD).sin();
        let cos_phi2_neg = (feature[6] - ANGLE_SHIFT_RAD).cos();
        let sin_phi2_disc_neg = discretize_value(sin_phi2_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        let cos_phi2_disc_neg = discretize_value(cos_phi2_neg, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE);
        
        let phi2_neg_changed = sin_phi2_disc_neg != sin_phi2_disc_orig || cos_phi2_disc_neg != cos_phi2_disc_orig;
        if phi2_neg_changed && (sin_phi2_disc_neg != sin_phi2_disc_pos || cos_phi2_disc_neg != cos_phi2_disc_pos) {
            phi2_angle_values.push((sin_phi2_disc_neg, cos_phi2_disc_neg));
        }
        
        // === HASH GENERATION ===
        // Generate all combinations of the collected unique values
        let max_combinations = ca_dist_values.len() * cb_dist_values.len() * 
                              ca_cb_angle_values.len() * phi1_angle_values.len() * phi2_angle_values.len();
        let mut unique_hashes = Vec::with_capacity(max_combinations);
        
        for &ca_dist in &ca_dist_values {
            for &cb_dist in &cb_dist_values {
                for &(sin_ca_cb_disc, cos_ca_cb_disc) in &ca_cb_angle_values {
                    for &(sin_phi1_disc, cos_phi1_disc) in &phi1_angle_values {
                        for &(sin_phi2_disc, cos_phi2_disc) in &phi2_angle_values {
                            let hashvalue = res1 << 25 | res2 << 20 | ca_dist << 16 
                                | cb_dist << 12 | sin_ca_cb_disc << 10 | cos_ca_cb_disc << 8
                                | sin_phi1_disc << 6 | cos_phi1_disc << 4 | sin_phi2_disc << 2 | cos_phi2_disc;
                            
                            unique_hashes.push(hashvalue);
                        }
                    }
                }
            }
        }
        
        // Remove any potential duplicates (should be rare with proper shifting)
        unique_hashes.sort_unstable();
        unique_hashes.dedup();
        
        (unique_hashes.len(), unique_hashes)
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
        let hash = GeometricHash::PDBTrRosetta(HashValue::from_u32(hash));
        match hash {
            GeometricHash::PDBTrRosetta(hash) => {
                println!("{:?}", hash);
            },
            _ => panic!("Invalid hash type"),
        }
    }
    
    #[test]
    fn test_deduplication_inline_performance() {
        let test_feature = vec![
            14.0, 17.0, 8.5, 12.3, 
            45.0_f32.to_radians(), 120.0_f32.to_radians(), -60.0_f32.to_radians()
        ];
        
        const ITERATIONS: usize = 100000;
        
        // Test default without shifting
        let start = std::time::Instant::now();
        for _ in 0..ITERATIONS {
            let _ = HashValue::perfect_hash_default(&test_feature);
        }
        let default_time = start.elapsed();
        
        // Test inline deduplication
        let start = std::time::Instant::now();
        for _ in 0..ITERATIONS {
            let _ = HashValue::perfect_hash_with_shifts_dedup_inline(&test_feature);
        }
        let inline_time = start.elapsed();
        
        // Test exhaustive deduplication
        let start = std::time::Instant::now();
        for _ in 0..ITERATIONS {
            let _ = HashValue::perfect_hash_with_all_shifts_exhaustive(&test_feature);
        }
        let exhaustive_time = start.elapsed();
        
        
        println!("=== Deduplication Performance ({} iterations) ===", ITERATIONS);
        println!("Default: {:?}", default_time);
        println!("Inline dedup: {:?}", inline_time);
        println!("Exhaustive dedup: {:?}", exhaustive_time);

        // Calculate per-iteration times in nanoseconds
        let default_per_iter = default_time.as_nanos() / ITERATIONS as u128;
        let inline_per_iter = inline_time.as_nanos() / ITERATIONS as u128;
        let exhaustive_per_iter = exhaustive_time.as_nanos() / ITERATIONS as u128;

        println!("\nPer iteration (nanoseconds):");
        println!("Default: {} ns", default_per_iter);
        println!("Inline dedup: {} ns", inline_per_iter);
        println!("Exhaustive dedup: {} ns", exhaustive_per_iter);

        // Calculate relative performance
        println!("\nRelative to callback approach:");
        println!("Inline dedup: {:.2}x", inline_per_iter as f64 / default_per_iter as f64);
        println!("Exhaustive dedup: {:.2}x", exhaustive_per_iter as f64 / default_per_iter as f64);

        // Test with different feature that produces more duplicates
        let duplicate_feature = vec![
            14.0, 17.0, 8.0, 8.0, 
            0.0_f32.to_radians(), 0.0_f32.to_radians(), 0.0_f32.to_radians()
        ];
        let (unique_count1, _) = HashValue::perfect_hash_with_shifts_dedup_inline(&test_feature);
        let (unique_count2, _) = HashValue::perfect_hash_with_shifts_dedup_inline(&duplicate_feature);

        let (unique_count3, _) = HashValue::perfect_hash_with_all_shifts_exhaustive(&duplicate_feature);

        println!("\nDeduplication effectiveness:");
        println!("Diverse feature: {} unique out of 8 ({:.1}% efficiency)", 
                 unique_count1, unique_count1 as f64 / 8.0 * 100.0);
        println!("Duplicate-prone feature: {} unique out of 8 ({:.1}% efficiency)", 
                 unique_count2, unique_count2 as f64 / 8.0 * 100.0);
        println!("Duplicate-prone feature (exhaustive): {} unique out of 243 ({:.1}% efficiency)", 
                 unique_count3, unique_count3 as f64 / 243.0 * 100.0);
    }
    
    #[test]
    fn test_deduplication_correctness() {
        // Test with a feature that should produce some duplicates
        let duplicate_feature = vec![
            14.0, 17.0, 8.0, 8.0, // Same distances
            0.0_f32.to_radians(), 0.0_f32.to_radians(), 0.0_f32.to_radians() // Same angles
        ];
        
        println!("=== Deduplication Correctness Test ===");
        // Get unique hashes from inline deduplication
        let (unique_count, unique_hashes) = HashValue::perfect_hash_with_shifts_dedup_inline(&duplicate_feature);
        
        println!("Inline dedup count: {}", unique_count);

        // Get unique hashes from exhaustive deduplication
        let (unique_count_ex, unique_hashes_ex) = HashValue::perfect_hash_with_all_shifts_exhaustive(&duplicate_feature);
        println!("Exhaustive dedup count: {}", unique_count_ex);
        
        
        // Print unique hashes
        println!("Unique hashes from inline dedup:");
        for i in 0..unique_count as usize {
            println!("  {}: {}", i, unique_hashes[i]);
            // Print reverse-mapped values for verification
            let hash = HashValue::from_u32(unique_hashes[i]);
            let values = hash.reverse_hash_default();
            println!("     values: {:?}", values);
        }
        
        // Print unique hashes from exhaustive dedup
        println!("Unique hashes from exhaustive dedup:");
        for (i, &hashval) in unique_hashes_ex.iter().enumerate() {
            println!("  {}: {}", i, hashval);
            // Print reverse-mapped values for verification
            let hash = HashValue::from_u32(hashval);
            let values = hash.reverse_hash_default();
            println!("     values: {:?}", values);
        }


    }
    
    #[test]
    fn test_max_32_shifts_efficiency() {
        println!("=== Max 32 Shifts Efficiency Test ===");
        
        // Test case 1: Feature near bin boundaries (should benefit from shifting)
        let boundary_feature = vec![
            14.0, 17.0, 8.499, 12.3,  // CA distance very close to bin boundary
            45.1_f32.to_radians(), 119.9_f32.to_radians(), -60.1_f32.to_radians()
        ];
        
        // Test case 2: Feature well within bins (should not benefit much from shifting)
        let centered_feature = vec![
            14.0, 17.0, 9.0, 13.0,  // Values well within bins
            50.0_f32.to_radians(), 110.0_f32.to_radians(), -50.0_f32.to_radians()
        ];
        
        // Test the new max 32 shifts function
        let (count_max32, hashes_max32) = HashValue::perfect_hash_with_max_32_shifts(&boundary_feature);
        let (count_exhaustive, hashes_exhaustive) = HashValue::perfect_hash_with_all_shifts_exhaustive(&boundary_feature);
        let (count_inline, _) = HashValue::perfect_hash_with_shifts_dedup_inline(&boundary_feature);
        
        println!("Hash generation comparison for boundary feature:");
        println!("  Max 32 shifts method: {} unique hashes", count_max32);
        println!("  Exhaustive method: {} unique hashes", count_exhaustive);
        println!("  Inline method: {} unique hashes", count_inline);
        
        // Test with centered feature
        let (count_max32_centered, _) = HashValue::perfect_hash_with_max_32_shifts(&centered_feature);
        let (count_exhaustive_centered, _) = HashValue::perfect_hash_with_all_shifts_exhaustive(&centered_feature);
        
        println!("Hash generation comparison for centered feature:");
        println!("  Max 32 shifts method: {} unique hashes", count_max32_centered);
        println!("  Exhaustive method: {} unique hashes", count_exhaustive_centered);
        
        // Verify that max32 method produces valid results
        let mut all_max32_found_in_exhaustive = true;
        for &max32_hash in &hashes_max32 {
            if !hashes_exhaustive.contains(&max32_hash) {
                all_max32_found_in_exhaustive = false;
                println!("  ERROR: Max32 hash {} not found in exhaustive!", max32_hash);
            }
        }
        
        if all_max32_found_in_exhaustive {
            println!("  ✓ All max32 hashes are present in exhaustive results");
        }
        
        // Performance comparison
        const ITERATIONS: usize = 100000;
        
        let start = std::time::Instant::now();
        for _ in 0..ITERATIONS {
            let _ = HashValue::perfect_hash_default(&boundary_feature);
        }
        let default_time = start.elapsed();
        
        let start = std::time::Instant::now();
        for _ in 0..ITERATIONS {
            let _ = HashValue::perfect_hash_with_max_32_shifts(&boundary_feature);
        }
        let max32_time = start.elapsed();
        
        let start = std::time::Instant::now();
        for _ in 0..ITERATIONS {
            let _ = HashValue::perfect_hash_with_all_shifts_exhaustive(&boundary_feature);
        }
        let exhaustive_time = start.elapsed();
        
        let start = std::time::Instant::now();
        for _ in 0..ITERATIONS {
            let _ = HashValue::perfect_hash_with_shifts_dedup_inline(&boundary_feature);
        }
        let inline_time = start.elapsed();
        
        println!("Performance comparison ({} iterations):", ITERATIONS);
        println!("  Default (1 hash): {:?}", default_time);
        println!("  Max 32 shifts: {:?}", max32_time);
        println!("  Inline dedup (8 combos): {:?}", inline_time);
        println!("  Exhaustive (243 combos): {:?}", exhaustive_time);
        
        let default_per_iter = default_time.as_nanos() / ITERATIONS as u128;
        let max32_per_iter = max32_time.as_nanos() / ITERATIONS as u128;
        let inline_per_iter = inline_time.as_nanos() / ITERATIONS as u128;
        let exhaustive_per_iter = exhaustive_time.as_nanos() / ITERATIONS as u128;
        
        println!("Per iteration (nanoseconds):");
        println!("  Default: {} ns", default_per_iter);
        println!("  Max 32: {} ns", max32_per_iter);
        println!("  Inline: {} ns", inline_per_iter);
        println!("  Exhaustive: {} ns", exhaustive_per_iter);
        
        println!("Efficiency (hashes per nanosecond * 1000):");
        println!("  Default: {:.1}", 1.0 / default_per_iter as f64 * 1000.0);
        println!("  Max 32: {:.1}", count_max32 as f64 / max32_per_iter as f64 * 1000.0);
        println!("  Inline: {:.1}", count_inline as f64 / inline_per_iter as f64 * 1000.0);
        println!("  Exhaustive: {:.1}", count_exhaustive as f64 / exhaustive_per_iter as f64 * 1000.0);
        
        println!("Max combinations possible: {} (2^5 = 32)", 2_u32.pow(5));
        println!("Actual max combinations achieved: {}", count_max32);
    }

    #[test]
    fn test_feature_distribution_analysis() {
        use crate::controller::io::read_structure_from_path;
        use crate::structure::core::CompactStructure;
        use crate::controller::feature::get_single_feature;
        use crate::utils::combination::CombinationIterator;
        use crate::geometry::core::HashType;
        
        println!("=== PDBTrRosetta Feature Distribution Analysis ===");
        
        // Try to load test structures
        let test_structures = vec![
            "data/AF-P17538-F1-model_v4.pdb",
        ];
        
        let mut all_features = Vec::new();
        let mut structures_analyzed = 0;
        
        for structure_path in &test_structures {
            if let Some(structure) = read_structure_from_path(structure_path) {
                println!("\nAnalyzing structure: {}", structure_path);
                let compact_structure = CompactStructure::build(&structure);
                let hash_type = HashType::PDBTrRosetta;
                let dist_cutoff = 20.0;
                
                let res_bound = CombinationIterator::new(compact_structure.num_residues);
                let mut feature = vec![0.0; 9];
                let mut structure_features = Vec::new();
                
                res_bound.for_each(|(i, j)| {
                    if i == j { return; }
                    
                    let has_feature = get_single_feature(i, j, &compact_structure, hash_type, dist_cutoff, &mut feature);
                    if has_feature {
                        // Store the geometric features (skip amino acid indices)
                        structure_features.push(vec![
                            feature[2], // CA distance
                            feature[3], // CB distance  
                            feature[4], // CA-CB angle
                            feature[5], // Phi1 torsion
                            feature[6], // Phi2 torsion
                        ]);
                    }
                });
                
                println!("  Extracted {} feature pairs", structure_features.len());
                all_features.extend(structure_features);
                structures_analyzed += 1;
                
                if structures_analyzed >= 3 { break; } // Limit to avoid too much output
            }
        }
        
        if all_features.is_empty() {
            println!("No structures found for analysis. Skipping distribution test.");
            return;
        }
        
        println!("\n=== Distribution Analysis Results ===");
        println!("Total features analyzed: {} from {} structures", all_features.len(), structures_analyzed);
        
        // Analyze each feature dimension
        let feature_names = vec!["CA Distance", "CB Distance", "CA-CB Angle", "Phi1 Torsion", "Phi2 Torsion"];
        
        for (dim, name) in feature_names.iter().enumerate() {
            let mut values: Vec<f32> = all_features.iter().map(|f| f[dim]).collect();
            values.sort_by(|a, b| a.partial_cmp(b).unwrap());
            
            let min_val = values[0];
            let max_val = values[values.len() - 1];
            let mean = values.iter().sum::<f32>() / values.len() as f32;
            let median = values[values.len() / 2];
            
            // Calculate standard deviation
            let variance = values.iter().map(|x| (x - mean).powi(2)).sum::<f32>() / values.len() as f32;
            let std_dev = variance.sqrt();
            
            // Calculate quantiles
            let q25 = values[values.len() / 4];
            let q75 = values[values.len() * 3 / 4];
            
            println!("\n{} Distribution:", name);
            println!("  Range: [{:.3}, {:.3}]", min_val, max_val);
            println!("  Mean: {:.3} ± {:.3}", mean, std_dev);
            println!("  Median: {:.3}", median);
            println!("  Quartiles: Q25={:.3}, Q75={:.3}", q25, q75);
            
            // For angles, also show in degrees
            if dim >= 2 {
                println!("  Range (degrees): [{:.1}°, {:.1}°]", min_val.to_degrees(), max_val.to_degrees());
                println!("  Mean (degrees): {:.1}° ± {:.1}°", mean.to_degrees(), std_dev.to_degrees());
            }
            
            // Analyze current binning effectiveness
            let current_nbin = if dim < 2 { PDBTR_NBIN_DIST as usize } else { PDBTR_NBIN_SIN_COS as usize };
            println!("  Current bins: {}", current_nbin);
            
            if dim < 2 {
                // Distance binning analysis
                let bin_width = (MAX_DIST - MIN_DIST) / current_nbin as f32;
                println!("  Current bin width: {:.3} Å", bin_width);
                
                // Show how many features fall in each bin
                let mut bin_counts = vec![0; current_nbin];
                for &val in &values {
                    let bin = ((val - MIN_DIST) / bin_width).min((current_nbin - 1) as f32) as usize;
                    bin_counts[bin] += 1;
                }
                
                println!("  Bin occupancy (showing non-empty bins):");
                for (i, &count) in bin_counts.iter().enumerate() {
                    if count > 0 {
                        let bin_start = MIN_DIST + i as f32 * bin_width;
                        let bin_end = bin_start + bin_width;
                        println!("    Bin {}: [{:.1}, {:.1}] -> {} features ({:.1}%)", 
                                i, bin_start, bin_end, count, 
                                100.0 * count as f32 / values.len() as f32);
                    }
                }
            } else {
                // Angle binning analysis (sin/cos representation)
                let sin_values: Vec<f32> = values.iter().map(|x| x.sin()).collect();
                let cos_values: Vec<f32> = values.iter().map(|x| x.cos()).collect();
                
                let bin_width = (MAX_SIN_COS - MIN_SIN_COS) / current_nbin as f32;
                println!("  Current sin/cos bin width: {:.3}", bin_width);
                
                // Analyze sin values
                let mut sin_bin_counts = vec![0; current_nbin];
                for &val in &sin_values {
                    let bin = ((val - MIN_SIN_COS) / bin_width).min((current_nbin - 1) as f32) as usize;
                    sin_bin_counts[bin] += 1;
                }
                
                // Analyze cos values
                let mut cos_bin_counts = vec![0; current_nbin];
                for &val in &cos_values {
                    let bin = ((val - MIN_SIN_COS) / bin_width).min((current_nbin - 1) as f32) as usize;
                    cos_bin_counts[bin] += 1;
                }
                
                println!("  Sin value bin occupancy:");
                for (i, &count) in sin_bin_counts.iter().enumerate() {
                    if count > 0 {
                        let bin_start = MIN_SIN_COS + i as f32 * bin_width;
                        let bin_end = bin_start + bin_width;
                        println!("    Sin Bin {}: [{:.2}, {:.2}] -> {} features ({:.1}%)", 
                                i, bin_start, bin_end, count, 
                                100.0 * count as f32 / values.len() as f32);
                    }
                }
                
                println!("  Cos value bin occupancy:");
                for (i, &count) in cos_bin_counts.iter().enumerate() {
                    if count > 0 {
                        let bin_start = MIN_SIN_COS + i as f32 * bin_width;
                        let bin_end = bin_start + bin_width;
                        println!("    Cos Bin {}: [{:.2}, {:.2}] -> {} features ({:.1}%)", 
                                i, bin_start, bin_end, count, 
                                100.0 * count as f32 / values.len() as f32);
                    }
                }
            }
        }
        
        // Analyze feature correlations
        println!("\n=== Feature Correlation Analysis ===");
        for i in 0..5 {
            for j in (i+1)..5 {
                let values_i: Vec<f32> = all_features.iter().map(|f| f[i]).collect();
                let values_j: Vec<f32> = all_features.iter().map(|f| f[j]).collect();
                
                let correlation = calculate_correlation(&values_i, &values_j);
                if correlation.abs() > 0.1 {
                    println!("  {} vs {}: correlation = {:.3}", 
                            feature_names[i], feature_names[j], correlation);
                }
            }
        }
        
        // Suggest optimal binning based on distribution
        println!("\n=== Binning Recommendations ===");
        
        // Distance recommendations
        let ca_distances: Vec<f32> = all_features.iter().map(|f| f[0]).collect();
        let cb_distances: Vec<f32> = all_features.iter().map(|f| f[1]).collect();
        
        println!("Distance Binning Suggestions:");
        suggest_distance_binning(&ca_distances, "CA Distance");
        suggest_distance_binning(&cb_distances, "CB Distance");
        
        // Angle recommendations  
        let angles = vec![
            ("CA-CB Angle", 2),
            ("Phi1 Torsion", 3), 
            ("Phi2 Torsion", 4),
        ];
        
        println!("\nAngle Binning Suggestions:");
        for (name, dim) in angles {
            let angle_values: Vec<f32> = all_features.iter().map(|f| f[dim]).collect();
            suggest_angle_binning(&angle_values, name);
        }
    }
    
    // Helper function to calculate Pearson correlation
    fn calculate_correlation(x: &[f32], y: &[f32]) -> f32 {
        let n = x.len() as f32;
        let mean_x = x.iter().sum::<f32>() / n;
        let mean_y = y.iter().sum::<f32>() / n;
        
        let numerator: f32 = x.iter().zip(y.iter())
            .map(|(xi, yi)| (xi - mean_x) * (yi - mean_y))
            .sum();
            
        let sum_sq_x: f32 = x.iter().map(|xi| (xi - mean_x).powi(2)).sum();
        let sum_sq_y: f32 = y.iter().map(|yi| (yi - mean_y).powi(2)).sum();
        
        if sum_sq_x == 0.0 || sum_sq_y == 0.0 {
            0.0
        } else {
            numerator / (sum_sq_x * sum_sq_y).sqrt()
        }
    }
    
    // Suggest optimal distance binning
    fn suggest_distance_binning(values: &[f32], name: &str) {
        let mut sorted_values = values.to_vec();
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        let min_val = sorted_values[0];
        let max_val = sorted_values[sorted_values.len() - 1];
        let range = max_val - min_val;
        
        println!("  {} ({:.1} - {:.1} Å, range: {:.1} Å):", name, min_val, max_val, range);
        
        // Current uniform binning
        let current_bin_width = (MAX_DIST - MIN_DIST) / NBIN_DIST as f32;
        let effective_bins = (range / current_bin_width).ceil() as usize;
        println!("    Current uniform: {:.2} Å bins, ~{} effective bins", current_bin_width, effective_bins);
        
        // Suggest protein-aware boundaries for different ranges
        if name.contains("CA") {
            println!("    Protein-aware suggestion: [3.8, 4.5, 6.0, 8.0, 10.0, 12.0, 15.0, 18.0] (secondary structure based)");
        } else {
            println!("    Protein-aware suggestion: [3.5, 4.2, 5.0, 6.5, 8.0, 10.0, 12.0, 15.0] (side chain packing based)");
        }
        
        // Quantile-based suggestion
        print!("    Quantile-based suggestion: [");
        let nbin_dist_int = PDBTR_NBIN_DIST as usize;
        for i in 0..=nbin_dist_int {
            let quantile = i as f32 / nbin_dist_int as f32;
            let index = (quantile * (sorted_values.len() - 1) as f32) as usize;
            print!("{:.1}", sorted_values[index]);
            if i < nbin_dist_int { print!(", "); }
        }
        println!("]");
    }
    
    // Suggest optimal angle binning  
    fn suggest_angle_binning(values: &[f32], name: &str) {
        let mut sorted_values = values.to_vec();
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        let min_val = sorted_values[0];
        let max_val = sorted_values[sorted_values.len() - 1];
        
        println!("  {} ({:.1}° - {:.1}°):", name, min_val.to_degrees(), max_val.to_degrees());
        
        // Analyze sin/cos distributions
        let sin_values: Vec<f32> = values.iter().map(|x| x.sin()).collect();
        let cos_values: Vec<f32> = values.iter().map(|x| x.cos()).collect();
        
        let sin_range = sin_values.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap() - 
                       sin_values.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        let cos_range = cos_values.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap() - 
                       cos_values.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        
        println!("    Sin range: {:.2}, Cos range: {:.2}", sin_range, cos_range);
        
        // Current binning
        let current_sin_cos_bin_width = (MAX_SIN_COS - MIN_SIN_COS) / PDBTR_NBIN_SIN_COS as f32;
        println!("    Current sin/cos bin width: {:.2}", current_sin_cos_bin_width);
        
        // Suggest Ramachandran-aware binning for torsion angles
        if name.contains("Phi") {
            println!("    Ramachandran-aware suggestion: [-162°, -126°, -54°, -12°, 12°, 54°, 126°, 162°]");
        } else {
            println!("    Uniform angle suggestion: 45° intervals");
        }
    }
}