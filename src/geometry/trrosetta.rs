// Geometric features from trRosetta paper

use std::fmt;
use crate::geometry::core::HashType;
use crate::geometry::util::discretize_f32_value_into_u64 as discretize_value;
use crate::geometry::util::continuize_u64_value_into_f32 as continuize_value;
use crate::geometry::util::*;

// Constants
const NBIN_OMEGA_SIN_COS: f32 = 3.0;
const NBIN_THETA_SIN_COS: f32 = 3.0;
const NBIN_PHI_SIN_COS: f32 = 3.0;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u64);

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
    
    pub fn hash_type(&self) -> HashType {
        HashType::TrRosetta
    }
    
    pub fn perfect_hash(feature: Vec<f32>, nbin_dist: usize, nbin_angle: usize) -> Self {
        let mut cb_dist = feature[0];
        if cb_dist > 20.0 {
            cb_dist = 20.0;
        }
        let sin_angles = feature.iter().skip(1).map(|&x| x.sin()).collect::<Vec<f32>>();
        let cos_angles = feature.iter().skip(1).map(|&x| x.cos()).collect::<Vec<f32>>();

        let nbin_dist = if nbin_dist > 16 { 16.0 } else { nbin_dist as f32 };
        let nbin_angle = if nbin_angle > 16 { 16.0 } else { nbin_angle as f32 };
        let h_cb_dist = discretize_value(cb_dist, MIN_DIST, MAX_DIST, nbin_dist);
        let h_sin_omega = discretize_value(sin_angles[0], MIN_SIN_COS, MAX_SIN_COS, nbin_angle);
        let h_sin_theta1 = discretize_value(sin_angles[1], MIN_SIN_COS, MAX_SIN_COS, nbin_angle);
        let h_sin_theta2 = discretize_value(sin_angles[2], MIN_SIN_COS, MAX_SIN_COS, nbin_angle);
        let h_sin_phi1 = discretize_value(sin_angles[3], MIN_SIN_COS, MAX_SIN_COS, nbin_angle);
        let h_sin_phi2 = discretize_value(sin_angles[4], MIN_SIN_COS, MAX_SIN_COS, nbin_angle);
        let h_cos_omega = discretize_value(cos_angles[0], MIN_SIN_COS, MAX_SIN_COS, nbin_angle);
        let h_cos_theta1 = discretize_value(cos_angles[1], MIN_SIN_COS, MAX_SIN_COS, nbin_angle);
        let h_cos_theta2 = discretize_value(cos_angles[2], MIN_SIN_COS, MAX_SIN_COS, nbin_angle);
        let h_cos_phi1 = discretize_value(cos_angles[3], MIN_SIN_COS, MAX_SIN_COS, nbin_angle);
        let h_cos_phi2 = discretize_value(cos_angles[4], MIN_SIN_COS, MAX_SIN_COS, nbin_angle);
        
        let hashvalue = h_cb_dist << 40
            | h_sin_omega << 36 | h_cos_omega << 32
            | h_sin_theta1 << 28 | h_cos_theta1 << 24
            | h_sin_theta2 << 20 | h_cos_theta2 << 16
            | h_sin_phi1 << 12 | h_cos_phi1 << 8
            | h_sin_phi2 << 4 | h_cos_phi2;
        HashValue(hashvalue)
    }
    
    pub fn perfect_hash_default(feature: Vec<f32>) -> Self {
        let mut cb_dist = feature[0];
        if cb_dist > 20.0 {
            cb_dist = 20.0;
        }
        let sin_angles = feature.iter().skip(1).map(|&x| x.sin()).collect::<Vec<f32>>();
        let cos_angles = feature.iter().skip(1).map(|&x| x.cos()).collect::<Vec<f32>>();

        let h_cb_dist = discretize_value(cb_dist, MIN_DIST, MAX_DIST, NBIN_DIST);
        let h_sin_omega = discretize_value(sin_angles[0], MIN_SIN_COS, MAX_SIN_COS, NBIN_OMEGA_SIN_COS);
        let h_sin_theta1 = discretize_value(sin_angles[1], MIN_SIN_COS, MAX_SIN_COS, NBIN_THETA_SIN_COS);
        let h_sin_theta2 = discretize_value(sin_angles[2], MIN_SIN_COS, MAX_SIN_COS, NBIN_THETA_SIN_COS);
        let h_sin_phi1 = discretize_value(sin_angles[3], MIN_SIN_COS, MAX_SIN_COS, NBIN_PHI_SIN_COS);
        let h_sin_phi2 = discretize_value(sin_angles[4], MIN_SIN_COS, MAX_SIN_COS, NBIN_PHI_SIN_COS);
        let h_cos_omega = discretize_value(cos_angles[0], MIN_SIN_COS, MAX_SIN_COS, NBIN_OMEGA_SIN_COS);
        let h_cos_theta1 = discretize_value(cos_angles[1], MIN_SIN_COS, MAX_SIN_COS, NBIN_THETA_SIN_COS);
        let h_cos_theta2 = discretize_value(cos_angles[2], MIN_SIN_COS, MAX_SIN_COS, NBIN_THETA_SIN_COS);
        let h_cos_phi1 = discretize_value(cos_angles[3], MIN_SIN_COS, MAX_SIN_COS, NBIN_PHI_SIN_COS);
        let h_cos_phi2 = discretize_value(cos_angles[4], MIN_SIN_COS, MAX_SIN_COS, NBIN_PHI_SIN_COS);
        
        let hashvalue = h_cb_dist << 40
            | h_sin_omega << 36 | h_cos_omega << 32
            | h_sin_theta1 << 28 | h_cos_theta1 << 24
            | h_sin_theta2 << 20 | h_cos_theta2 << 16
            | h_sin_phi1 << 12 | h_cos_phi1 << 8
            | h_sin_phi2 << 4 | h_cos_phi2;

        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self, nbin_dist: usize, nbin_angle: usize) -> Vec<f32> {
        let cb_dist = continuize_value(
            (self.0 >> 40) & BITMASK64_4BIT, MIN_DIST, MAX_DIST, nbin_dist as f32
        );
        let sin_omega = continuize_value(
            (self.0 >> 36) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_omega = continuize_value(
            (self.0 >> 32) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let sin_theta1 = continuize_value(
            (self.0 >> 28) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_theta1 = continuize_value(
            (self.0 >> 24) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let sin_theta2 = continuize_value(
            (self.0 >> 20) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_theta2 = continuize_value(
            (self.0 >> 16) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let sin_phi1 = continuize_value(
            (self.0 >> 12) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_phi1 = continuize_value(
            (self.0 >> 8) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let sin_phi2 = continuize_value(
            (self.0 >> 4) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let cos_phi2 = continuize_value(
            self.0 & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, nbin_angle as f32
        );
        let angle_omega = sin_omega.atan2(cos_omega).to_degrees();
        let angle_theta1 = sin_theta1.atan2(cos_theta1).to_degrees();
        let angle_theta2 = sin_theta2.atan2(cos_theta2).to_degrees();
        let angle_phi1 = sin_phi1.atan2(cos_phi1).to_degrees();
        let angle_phi2 = sin_phi2.atan2(cos_phi2).to_degrees();
        vec![cb_dist, angle_omega, angle_theta1, angle_theta2, angle_phi1, angle_phi2]
    }
    
    pub fn reverse_hash_default(&self) -> Vec<f32> {
        let cb_dist = ((self.0 >> 40) & BITMASK64_4BIT) as f32;
        let sin_omega = continuize_value(
            (self.0 >> 36) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, NBIN_OMEGA_SIN_COS
        );
        let cos_omega = continuize_value(
            (self.0 >> 32) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, NBIN_OMEGA_SIN_COS
        );
        let sin_theta1 = continuize_value(
            (self.0 >> 28) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, NBIN_THETA_SIN_COS
        );
        let cos_theta1 = continuize_value(
            (self.0 >> 24) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, NBIN_THETA_SIN_COS
        );
        let sin_theta2 = continuize_value(
            (self.0 >> 20) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, NBIN_THETA_SIN_COS
        );
        let cos_theta2 = continuize_value(
            (self.0 >> 16) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, NBIN_THETA_SIN_COS
        );
        let sin_phi1 = continuize_value(
            (self.0 >> 12) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, NBIN_PHI_SIN_COS
        );
        let cos_phi1 = continuize_value(
            (self.0 >> 8) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, NBIN_PHI_SIN_COS
        );
        let sin_phi2 = continuize_value(
            (self.0 >> 4) & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, NBIN_PHI_SIN_COS
        );
        let cos_phi2 = continuize_value(
            self.0 & BITMASK64_4BIT, MIN_SIN_COS, MAX_SIN_COS, NBIN_PHI_SIN_COS
        );
        let angle_omega = sin_omega.atan2(cos_omega).to_degrees();
        let angle_theta1 = sin_theta1.atan2(cos_theta1).to_degrees();
        let angle_theta2 = sin_theta2.atan2(cos_theta2).to_degrees();
        let angle_phi1 = sin_phi1.atan2(cos_phi1).to_degrees();
        let angle_phi2 = sin_phi2.atan2(cos_phi2).to_degrees();
        vec![cb_dist, angle_omega, angle_theta1, angle_theta2, angle_phi1, angle_phi2]
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

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let values = self.reverse_hash_default();
        write!(f, "HashValue({}), values={:?}", self.0, values)
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let val = self.reverse_hash_default();
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.0, val[0], val[1], val[2], val[3], val[4], val[5]
        )
    }   
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::util::map_aa_to_u8;
    use crate::geometry::core::GeometricHash;
    #[test]
    fn test_geometrichash_works() {
        // Test perfect hash
        let mut raw_feature = vec![14.0_f32, 15.9_f32, -116.0_f32, 30.0_f32, 45.0_f32, 87.0_f32];
        // Convert to radians
        raw_feature[1] = raw_feature[1].to_radians();
        raw_feature[2] = raw_feature[2].to_radians();
        raw_feature[3] = raw_feature[3].to_radians();
        raw_feature[4] = raw_feature[4].to_radians();
        raw_feature[5] = raw_feature[5].to_radians();
        let hash: GeometricHash = GeometricHash::TrRosetta(
            HashValue::perfect_hash(raw_feature, 8, 4)
        );
        let rev = hash.reverse_hash(8, 4);
        println!("{:?}", rev);
    }
}

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