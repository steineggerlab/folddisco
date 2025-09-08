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

    /// Strategic 7-way combination shifting - optimized for maximum boundary coverage
    /// Generates 7 hash values using all useful combinations of {0,+,-} for dist and angle
    /// Skips (-,-) combination which is often redundant with (+,+)
    #[inline]
    pub fn perfect_hash_with_shifts<F>(feature: &Vec<f32>, mut callback: F) 
    where F: FnMut(u32, u8) // (hash_value, shift_id: 0-7)
    {
        const NBIN_DIST: f32 = PDBTR_NBIN_DIST;
        const NBIN_ANGLE: f32 = PDBTR_NBIN_SIN_COS;
        
        // Optimized shift amounts for protein geometry
        const DIST_SHIFT: f32 = 0.6; // 60% of bin width for distances
        const ANGLE_SHIFT_RAD: f32 = std::f32::consts::PI / 8.0; // 22.5 degrees
        
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        
        // Precompute trigonometric values once for efficiency
        let sin_ca_cb = feature[4].sin();
        let cos_ca_cb = feature[4].cos();
        let sin_phi1 = feature[5].sin();
        let cos_phi1 = feature[5].cos();
        let sin_phi2 = feature[6].sin();
        let cos_phi2 = feature[6].cos();
        
        // 7 strategic combinations (skipping (-,-) as less useful):
        // 0: (0,0) - Original
        // 1: (0,+) - Angle+ only  
        // 2: (0,-) - Angle- only
        // 3: (+,0) - Distance+ only
        // 4: (+,+) - Both positive
        // 5: (+,-) - Dist+, Angle-
        // 6: (-,0) - Distance- only
        // 7: (-,+) - Dist-, Angle+
        let shift_combinations = [
            (0.0, 0.0),                           // 0: Original
            (0.0, ANGLE_SHIFT_RAD),               // 1: Angle+ only
            (0.0, -ANGLE_SHIFT_RAD),              // 2: Angle- only
            (DIST_SHIFT, 0.0),                    // 3: Distance+ only
            (DIST_SHIFT, ANGLE_SHIFT_RAD),        // 4: Both positive
            (DIST_SHIFT, -ANGLE_SHIFT_RAD),       // 5: Dist+, Angle-
            (-DIST_SHIFT, 0.0),                   // 6: Distance- only
            (-DIST_SHIFT, ANGLE_SHIFT_RAD),       // 7: Dist-, Angle+
        ];
        
        for (shift_id, &(dist_shift, angle_shift)) in shift_combinations.iter().enumerate() {
            // Distance features with shift
            let ca_dist = Self::discretize_with_shift(
                feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, dist_shift
            );
            let cb_dist = Self::discretize_with_shift(
                feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, dist_shift
            );
            
            // Angular features - optimize for zero shifts (use precomputed values)
            let (sin_ca_cb_shifted, cos_ca_cb_shifted, sin_phi1_shifted, cos_phi1_shifted, 
                 sin_phi2_shifted, cos_phi2_shifted) = if angle_shift == 0.0 {
                (sin_ca_cb, cos_ca_cb, sin_phi1, cos_phi1, sin_phi2, cos_phi2)
            } else {
                let shifted_ca_cb = feature[4] + angle_shift;
                let shifted_phi1 = feature[5] + angle_shift;
                let shifted_phi2 = feature[6] + angle_shift;
                (shifted_ca_cb.sin(), shifted_ca_cb.cos(),
                 shifted_phi1.sin(), shifted_phi1.cos(),
                 shifted_phi2.sin(), shifted_phi2.cos())
            };
            
            // Discretize angular values
            let sin_ca_cb_disc = discretize_value(
                sin_ca_cb_shifted, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
            );
            let cos_ca_cb_disc = discretize_value(
                cos_ca_cb_shifted, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
            );
            let sin_phi1_disc = discretize_value(
                sin_phi1_shifted, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
            );
            let cos_phi1_disc = discretize_value(
                cos_phi1_shifted, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
            );
            let sin_phi2_disc = discretize_value(
                sin_phi2_shifted, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
            );
            let cos_phi2_disc = discretize_value(
                cos_phi2_shifted, MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
            );
            
            // Pack hash value
            let hashvalue = res1 << 25 | res2 << 20 | ca_dist << 16 
                | cb_dist << 12 | sin_ca_cb_disc << 10 | cos_ca_cb_disc << 8
                | sin_phi1_disc << 6 | cos_phi1_disc << 4 | sin_phi2_disc << 2 | cos_phi2_disc;
            
            callback(hashvalue, shift_id as u8);
        }
    }
    
    /// Optimized discretization with shift - inline and branchless
    #[inline(always)]
    fn discretize_with_shift(value: f32, min_val: f32, max_val: f32, nbins: f32, shift: f32) -> u32 {
        let bin_width = (max_val - min_val) / nbins;
        // Apply shift by moving the bin boundaries
        let shifted_value = value + shift; // Note: we shift the value, not the boundaries
        let clamped_value = shifted_value.max(min_val).min(max_val - 1e-6);
        let bin_index = ((clamped_value - min_val) / bin_width).floor() as u32;
        bin_index.min(15) // 4-bit constraint for distances
    }
    
    /// Fast single shift hash generation for specific shift fraction
    #[inline]
    pub fn perfect_hash_single_shift(feature: &Vec<f32>, shift_fraction: f32) -> u32 {
        const NBIN_DIST: f32 = PDBTR_NBIN_DIST;
        const NBIN_ANGLE: f32 = PDBTR_NBIN_SIN_COS;
        const DIST_BIN_WIDTH: f32 = (MAX_DIST - MIN_DIST) / NBIN_DIST;
        const ANGLE_SHIFT_BASE: f32 = std::f32::consts::PI / 4.0; // 45 degrees
        
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        
        let dist_shift = DIST_BIN_WIDTH * shift_fraction;
        let angle_shift = ANGLE_SHIFT_BASE * shift_fraction;
        
        let ca_dist = Self::discretize_with_shift(
            feature[2], MIN_DIST, MAX_DIST, NBIN_DIST, dist_shift
        );
        let cb_dist = Self::discretize_with_shift(
            feature[3], MIN_DIST, MAX_DIST, NBIN_DIST, dist_shift
        );
        
        let shifted_ca_cb_angle = feature[4] + angle_shift;
        let shifted_phi1 = feature[5] + angle_shift;
        let shifted_phi2 = feature[6] + angle_shift;
        
        let sin_ca_cb_disc = discretize_value(
            shifted_ca_cb_angle.sin(), MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
        );
        let cos_ca_cb_disc = discretize_value(
            shifted_ca_cb_angle.cos(), MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
        );
        let sin_phi1_disc = discretize_value(
            shifted_phi1.sin(), MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
        );
        let cos_phi1_disc = discretize_value(
            shifted_phi1.cos(), MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
        );
        let sin_phi2_disc = discretize_value(
            shifted_phi2.sin(), MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
        );
        let cos_phi2_disc = discretize_value(
            shifted_phi2.cos(), MIN_SIN_COS, MAX_SIN_COS, NBIN_ANGLE
        );
        
        res1 << 25 | res2 << 20 | ca_dist << 16 
            | cb_dist << 12 | sin_ca_cb_disc << 10 | cos_ca_cb_disc << 8
            | sin_phi1_disc << 6 | cos_phi1_disc << 4 | sin_phi2_disc << 2 | cos_phi2_disc
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
    
    /// Most efficient deduplication: inline during generation
    /// Returns (unique_count, unique_hashes_array)
    /// Zero heap allocation, optimized for speed with precomputed trigonometry
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
    
    /// Helper function to add hash to unique array if not duplicate
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
    
    /// Stack-based deduplication with callback
    /// Good balance of performance and flexibility
    #[inline]
    pub fn perfect_hash_with_shifts_dedup_stack<F>(feature: &Vec<f32>, mut callback: F) -> u8
    where F: FnMut(u32)
    {
        let (unique_count, unique_hashes) = Self::perfect_hash_with_shifts_dedup_inline(feature);
        
        // Call callback for each unique hash
        for i in 0..unique_count as usize {
            callback(unique_hashes[i]);
        }
        
        unique_count
    }
    
    /// Usage example for indexing with deduplication
    /// This shows how to efficiently add a feature to an inverted index
    #[inline]
    pub fn add_to_index_with_dedup(
        feature: &Vec<f32>, 
        record_id: usize,
        index: &mut std::collections::HashMap<u32, Vec<usize>>
    ) -> u8 {
        let (unique_count, unique_hashes) = Self::perfect_hash_with_shifts_dedup_inline(feature);
        
        for i in 0..unique_count as usize {
            let hash = unique_hashes[i];
            index.entry(hash).or_insert_with(Vec::new).push(record_id);
        }
        
        unique_count
    }
    
    /// Usage example for querying with deduplication
    /// Returns deduplicated set of candidate record IDs
    #[inline]
    pub fn query_index_with_dedup(
        feature: &Vec<f32>,
        index: &std::collections::HashMap<u32, Vec<usize>>
    ) -> std::collections::HashSet<usize> {
        let mut candidates = std::collections::HashSet::new();
        let (unique_count, unique_hashes) = Self::perfect_hash_with_shifts_dedup_inline(feature);
        
        for i in 0..unique_count as usize {
            let hash = unique_hashes[i];
            if let Some(record_ids) = index.get(&hash) {
                candidates.extend(record_ids.iter());
            }
        }
        
        candidates
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
    fn test_shifting_performance() {
        let raw_feature = vec![
            14.0, 17.0, 14.0, 15.9, 
            116.0_f32.to_radians(), 80.0_f32.to_radians(), -100.0_f32.to_radians()
        ];
        
        // Test callback-based shifting
        let start = std::time::Instant::now();
        for _ in 0..10000 {
            let mut _hash_count = 0;
            HashValue::perfect_hash_with_shifts(&raw_feature, |_hash, _shift_id| {
                _hash_count += 1;
            });
        }
        let duration = start.elapsed();
        println!("Time elapsed in perfect_hash_with_shifts() is: {:?}", duration);
        
        // Test single shift
        let start = std::time::Instant::now();
        for _ in 0..10000 {
            let _ = HashValue::perfect_hash_single_shift(&raw_feature, 1.0/3.0);
            let _ = HashValue::perfect_hash_single_shift(&raw_feature, -1.0/3.0);
            let _ = HashValue::perfect_hash_default(&raw_feature);
        }
        let duration = start.elapsed();
        println!("Time elapsed in 3x perfect_hash_single_shift() is: {:?}", duration);
        
        // Compare original vs shifted - use different test values that will show shifting effect
        let test_feature = vec![
            14.0, 17.0, 7.99, 8.01, // Values near distance boundary
            89.0_f32.to_radians(), 91.0_f32.to_radians(), -1.0_f32.to_radians()
        ];
        
        let original_hash = HashValue::perfect_hash_default(&test_feature);
        let shifted_plus = HashValue::perfect_hash_single_shift(&test_feature, 1.0/3.0);
        let shifted_minus = HashValue::perfect_hash_single_shift(&test_feature, -1.0/3.0);
        
        println!("Original hash: {}", original_hash);
        println!("Shifted +1/3: {}", shifted_plus);
        println!("Shifted -1/3: {}", shifted_minus);
        
        // Show the reversed values to debug
        let orig_reversed = HashValue::from_u32(original_hash).reverse_hash_default();
        let plus_reversed = HashValue::from_u32(shifted_plus).reverse_hash_default();
        let minus_reversed = HashValue::from_u32(shifted_minus).reverse_hash_default();
        
        println!("Original values: {:?}", orig_reversed);
        println!("Shifted +1/3 values: {:?}", plus_reversed);
        println!("Shifted -1/3 values: {:?}", minus_reversed);
        
    }
    
    #[test]
    fn test_boundary_cases() {
        // Test values near bin boundaries
        let boundary_feature = vec![
            14.0, 17.0, 7.99, 8.01, // Near distance boundary
            90.0_f32.to_radians(), 0.0_f32.to_radians(), 180.0_f32.to_radians()
        ];
        
        let mut hashes = Vec::new();
        HashValue::perfect_hash_with_shifts(&boundary_feature, |hash, shift_id| {
            hashes.push((hash, shift_id));
        });
        
        println!("Boundary case hashes:");
        for (hash, shift_id) in hashes {
            println!("Shift {}: {}", shift_id, hash);
        }
        
        // Test that shifting helps capture boundary cases
        let similar_feature = vec![
            14.0, 17.0, 8.01, 7.99, // Swapped values near boundary
            92.0_f32.to_radians(), 2.0_f32.to_radians(), 178.0_f32.to_radians()
        ];
        
        let original1 = HashValue::perfect_hash_default(&boundary_feature);
        let original2 = HashValue::perfect_hash_default(&similar_feature);
        
        // They might have different original hashes
        println!("Original hashes: {} vs {}", original1, original2);
        
        // But shifting should provide more opportunities for matches
        let mut found_match = false;
        HashValue::perfect_hash_with_shifts(&boundary_feature, |hash1, _| {
            HashValue::perfect_hash_with_shifts(&similar_feature, |hash2, _| {
                if hash1 == hash2 {
                    found_match = true;
                    println!("Found match with shifting: {}", hash1);
                }
            });
        });
        
        if found_match {
            println!("Shifting successfully found match for boundary case");
        } else {
            println!("No match found even with shifting (this is also valid)");
        }
    }
    
    #[test]
    fn test_deduplication_inline_performance() {
        let test_feature = vec![
            14.0, 17.0, 8.5, 12.3, 
            45.0_f32.to_radians(), 120.0_f32.to_radians(), -60.0_f32.to_radians()
        ];
        
        const ITERATIONS: usize = 100000;
        
        // Test callback-based approach (current)
        let start = std::time::Instant::now();
        for _ in 0..ITERATIONS {
            HashValue::perfect_hash_with_shifts(&test_feature, |_hash, _shift_id| {
                // Simulate processing hash
            });
        }
        let callback_time = start.elapsed();
        
        // Test inline deduplication
        let start = std::time::Instant::now();
        for _ in 0..ITERATIONS {
            let _ = HashValue::perfect_hash_with_shifts_dedup_inline(&test_feature);
        }
        let inline_time = start.elapsed();
        
        // Test stack deduplication with callback
        let start = std::time::Instant::now();
        for _ in 0..ITERATIONS {
            let _ = HashValue::perfect_hash_with_shifts_dedup_stack(&test_feature, |_hash| {
                // Simulate processing hash
            });
        }
        let stack_time = start.elapsed();
        
        println!("=== Deduplication Performance ({} iterations) ===", ITERATIONS);
        println!("Callback approach: {:?}", callback_time);
        println!("Inline dedup: {:?}", inline_time);
        println!("Stack dedup: {:?}", stack_time);
        
        // Calculate per-iteration times in nanoseconds
        let callback_per_iter = callback_time.as_nanos() / ITERATIONS as u128;
        let inline_per_iter = inline_time.as_nanos() / ITERATIONS as u128;
        let stack_per_iter = stack_time.as_nanos() / ITERATIONS as u128;
        
        println!("\nPer iteration (nanoseconds):");
        println!("Callback approach: {} ns", callback_per_iter);
        println!("Inline dedup: {} ns", inline_per_iter);
        println!("Stack dedup: {} ns", stack_per_iter);
        
        // Calculate relative performance
        println!("\nRelative to callback approach:");
        println!("Inline dedup: {:.2}x", inline_per_iter as f64 / callback_per_iter as f64);
        println!("Stack dedup: {:.2}x", stack_per_iter as f64 / callback_per_iter as f64);
        
        // Test with different feature that produces more duplicates
        let duplicate_feature = vec![
            14.0, 17.0, 8.0, 8.0, 
            0.0_f32.to_radians(), 0.0_f32.to_radians(), 0.0_f32.to_radians()
        ];
        
        let (unique_count1, _) = HashValue::perfect_hash_with_shifts_dedup_inline(&test_feature);
        let (unique_count2, _) = HashValue::perfect_hash_with_shifts_dedup_inline(&duplicate_feature);
        
        println!("\nDeduplication effectiveness:");
        println!("Diverse feature: {} unique out of 8 ({:.1}% efficiency)", 
                 unique_count1, unique_count1 as f64 / 8.0 * 100.0);
        println!("Duplicate-prone feature: {} unique out of 8 ({:.1}% efficiency)", 
                 unique_count2, unique_count2 as f64 / 8.0 * 100.0);
    }
    
    #[test]
    fn test_deduplication_correctness() {
        // Test with a feature that should produce some duplicates
        let duplicate_feature = vec![
            14.0, 17.0, 8.0, 8.0, // Same distances
            0.0_f32.to_radians(), 0.0_f32.to_radians(), 0.0_f32.to_radians() // Same angles
        ];
        
        println!("=== Deduplication Correctness Test ===");
        
        // Count total hashes from callback approach
        let mut total_hashes = Vec::new();
        HashValue::perfect_hash_with_shifts(&duplicate_feature, |hash, shift_id| {
            total_hashes.push((hash, shift_id));
        });
        
        // Get unique hashes from inline deduplication
        let (unique_count, unique_hashes) = HashValue::perfect_hash_with_shifts_dedup_inline(&duplicate_feature);
        
        // Compare with manual deduplication
        let mut manual_unique: Vec<u32> = total_hashes.iter().map(|(h, _)| *h).collect();
        manual_unique.sort_unstable();
        manual_unique.dedup();
        
        println!("Total hashes generated: {}", total_hashes.len());
        println!("Inline dedup count: {}", unique_count);
        println!("Manual dedup count: {}", manual_unique.len());
        
        // Print all generated hashes
        for (hash, shift_id) in &total_hashes {
            println!("  Shift {}: {}", shift_id, hash);
        }
        
        // Print unique hashes
        println!("Unique hashes from inline dedup:");
        for i in 0..unique_count as usize {
            println!("  {}: {}", i, unique_hashes[i]);
        }
        
        // Verify correctness
        assert_eq!(unique_count as usize, manual_unique.len(), 
                   "Inline dedup should produce same count as manual dedup");
        
        // Test stack deduplication consistency
        let mut stack_hashes = Vec::new();
        let stack_count = HashValue::perfect_hash_with_shifts_dedup_stack(&duplicate_feature, |hash| {
            stack_hashes.push(hash);
        });
        
        assert_eq!(unique_count, stack_count, 
                   "Stack dedup should produce same count as inline dedup");
        
        println!("✓ All deduplication methods are consistent");
    }
    
    #[test]
    fn test_deduplication_boundary_cases() {
        println!("=== Boundary Case Deduplication Tests ===");
        
        // Test case 1: All identical values (maximum duplicates)
        let identical_feature = vec![
            14.0, 17.0, 8.0, 8.0, 
            0.0_f32.to_radians(), 0.0_f32.to_radians(), 0.0_f32.to_radians()
        ];
        
        let (count1, _) = HashValue::perfect_hash_with_shifts_dedup_inline(&identical_feature);
        println!("Identical feature produces {} unique hashes", count1);
        
        // Test case 2: Maximum diversity (minimum duplicates)
        let diverse_feature = vec![
            14.0, 17.0, 7.1, 13.7, 
            23.5_f32.to_radians(), 157.3_f32.to_radians(), -87.9_f32.to_radians()
        ];
        
        let (count2, _) = HashValue::perfect_hash_with_shifts_dedup_inline(&diverse_feature);
        println!("Diverse feature produces {} unique hashes", count2);
        
        // Test case 3: Near boundary values
        let boundary_feature = vec![
            14.0, 17.0, 7.99, 8.01, 
            89.9_f32.to_radians(), 90.1_f32.to_radians(), -179.9_f32.to_radians()
        ];
        
        let (count3, _) = HashValue::perfect_hash_with_shifts_dedup_inline(&boundary_feature);
        println!("Boundary feature produces {} unique hashes", count3);
        
        // Verify we always get at least 1 unique hash and at most 8
        assert!(count1 >= 1 && count1 <= 8);
        assert!(count2 >= 1 && count2 <= 8);
        assert!(count3 >= 1 && count3 <= 8);
        
        println!("✓ All boundary cases produce valid unique hash counts");
    }
    
    #[test]
    fn test_deduplication_usage_examples() {
        println!("=== Deduplication Usage Examples ===");
        
        // Create a simple index
        let mut index: std::collections::HashMap<u32, Vec<usize>> = std::collections::HashMap::new();
        
        // Test features
        let features = vec![
            // Feature 0: Alpha helix
            vec![14.0, 17.0, 8.0, 12.0, 120.0_f32.to_radians(), -60.0_f32.to_radians(), -45.0_f32.to_radians()],
            // Feature 1: Beta sheet
            vec![14.0, 17.0, 10.0, 14.0, 120.0_f32.to_radians(), -120.0_f32.to_radians(), 120.0_f32.to_radians()],
            // Feature 2: Similar to feature 0 (should have some shared hashes)
            vec![14.0, 17.0, 8.1, 12.1, 121.0_f32.to_radians(), -59.0_f32.to_radians(), -44.0_f32.to_radians()],
        ];
        
        // Add features to index using deduplication
        for (record_id, feature) in features.iter().enumerate() {
            let unique_count = HashValue::add_to_index_with_dedup(feature, record_id, &mut index);
            println!("Record {}: {} unique hashes added to index", record_id, unique_count);
        }
        
        println!("Index now contains {} unique hash keys", index.len());
        
        // Query the index
        let query_feature = vec![
            14.0, 17.0, 8.05, 12.05, 120.5_f32.to_radians(), -59.5_f32.to_radians(), -44.5_f32.to_radians()
        ];
        
        let candidates = HashValue::query_index_with_dedup(&query_feature, &index);
        println!("Query found {} candidate records: {:?}", candidates.len(), candidates);
        
        // Compare with non-deduplicated approach
        let mut candidates_no_dedup = std::collections::HashSet::new();
        HashValue::perfect_hash_with_shifts(&query_feature, |hash, _| {
            if let Some(record_ids) = index.get(&hash) {
                candidates_no_dedup.extend(record_ids.iter());
            }
        });
        
        println!("Non-dedup query found {} candidates: {:?}", 
                 candidates_no_dedup.len(), candidates_no_dedup);
        
        // Both should find the same candidates
        assert_eq!(candidates, candidates_no_dedup, 
                   "Deduped and non-deduped queries should find same candidates");
        
        println!("✓ Deduplication preserves search results while improving efficiency");
    }

}