// Default hasher using trrosetta features

// Implementation
// 5bit aa1, 5bit aa2, cbeta distance - 16 bins; 4 bits
// Angles - 6 bins each for sin and cos; total 25 bins; 4 bits for sin and cos each
 
use std::fmt;
use crate::geometry::core::HashType;
use crate::geometry::util::discretize_f32_value_into_u64 as discretize_value;
use crate::geometry::util::continuize_u64_value_into_f32 as continuize_value;

// Constants
// 1. for cb_dist
const MIN_DIST: f32 = 2.0;
const MAX_DIST: f32 = 20.0;
const NBIN_DIST: f32 = 16.0;
// 2. NEW IDEA for encoding angles; represent as sin and cos
const MIN_SIN_COS: f32 = -1.0;
const MAX_SIN_COS: f32 = 1.0;
const NBIN_SIN_COS: f32 = 6.0;
// Bitmasks
const BITMASK64_4BIT: u64 = 0x000000000000000F;
const BITMASK64_5BIT: u64 = 0x000000000000001F;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(pub u64);

impl HashValue {
    pub fn perfect_hash(feature: Vec<f32>) -> Self {
        let res1 = feature[0] as u64;
        let res2 = feature[1] as u64;
        HashValue::_perfect_hash(
            res1, res2, feature[2], feature[3],
            feature[4], feature[5], feature[6], feature[7]
        )
    }
    pub fn reverse_hash(&self) -> Vec<f32> {
        self._reverse_hash().to_vec()
    }
    pub fn hash_type(&self) -> super::core::HashType {
        HashType::FoldDiscoDefault
    }

    fn _perfect_hash(
        res1: u64, res2: u64, cb_dist: f32, omega: f32,
        theta1: f32, theta2: f32, phi1: f32, phi2: f32
    ) -> Self {
        let mut cbd = cb_dist;
        if cb_dist > 20.0 {
            cbd = 20.0;
        }
        let h_cb_dist = discretize_value(cbd, MIN_DIST, MAX_DIST, NBIN_DIST);
        
        // Convert angles to sin and cos
        let angles = [omega, theta1, theta2, phi1, phi2];
        let sin_cos_angles = angles.iter().map(
            |&x| (x.sin(), x.cos())
        ).collect::<Vec<(f32, f32)>>();
        // Discretize sin and cos
        let sin_cos_angles = sin_cos_angles.iter().map(
            |&(sin, cos)| {(
                discretize_value(sin, MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS),
                discretize_value(cos, MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS)
            )}
        ).collect::<Vec<(u64, u64)>>();
        
        // Combine all the hash values
        let hashvalue = res1 << 49 | res2 << 44 | h_cb_dist << 40
            | sin_cos_angles[0].0 << 36 | sin_cos_angles[0].1 << 32 // sin_omega, cos_omega
            | sin_cos_angles[1].0 << 28 | sin_cos_angles[1].1 << 24 // sin_theta1, cos_theta1
            | sin_cos_angles[2].0 << 20 | sin_cos_angles[2].1 << 16 // sin_theta2, cos_theta2
            | sin_cos_angles[3].0 << 12 | sin_cos_angles[3].1 << 8 // sin_phi1, cos_phi1
            | sin_cos_angles[4].0 << 4 | sin_cos_angles[4].1; // sin_phi2, cos_phi2

        HashValue(hashvalue)
    }

    fn _reverse_hash(&self) -> [f32; 8] {
        let res1 = (self.0 >> 49) as f32;
        // Mask bits
        let res2 = ((self.0 >> 44) & BITMASK64_5BIT) as f32;
        let cb_dist = ((self.0 >> 40) & BITMASK64_4BIT) as f32;
        let sin_cos_vec = [
            ((self.0 >> 36) & BITMASK64_4BIT), // sin_omega
            ((self.0 >> 32) & BITMASK64_4BIT), // cos_omega
            ((self.0 >> 28) & BITMASK64_4BIT), // sin_theta1
            ((self.0 >> 24) & BITMASK64_4BIT), // cos_theta1
            ((self.0 >> 20) & BITMASK64_4BIT), // sin_theta2
            ((self.0 >> 16) & BITMASK64_4BIT), // cos_theta2
            ((self.0 >> 12) & BITMASK64_4BIT), // sin_phi1
            ((self.0 >> 8) & BITMASK64_4BIT), // cos_phi1
            ((self.0 >> 4) & BITMASK64_4BIT), // sin_phi2
            (self.0 & BITMASK64_4BIT), // cos_phi2
        ];

        let sin_cos_vec = sin_cos_vec.iter().map(
            |&x| continuize_value(
                x, MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS
            )
        ).collect::<Vec<f32>>();
        
        // Restores original angles
        let omega = sin_cos_vec[0].atan2(sin_cos_vec[1]);
        let theta1 = sin_cos_vec[2].atan2(sin_cos_vec[3]);
        let theta2 = sin_cos_vec[4].atan2(sin_cos_vec[5]);
        let phi1 = sin_cos_vec[6].atan2(sin_cos_vec[7]);
        let phi2 = sin_cos_vec[8].atan2(sin_cos_vec[9]);
        [res1, res2, cb_dist, omega, theta1, theta2, phi1, phi2]
    }
    
    pub fn from_u64(hashvalue: u64) -> Self {
        HashValue(hashvalue)
    }
    pub fn as_u64(&self) -> u64 {
        self.0
    }
    pub fn as_usize(&self) -> usize {
        self.0 as usize
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
        let val = self.reverse_hash();
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.0, val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]
        )
        // write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::util::map_aa_to_u8;
    use crate::geometry::core::GeometricHash;
    #[test]
    fn test_default_hash_works() {
        // Test perfect hash
        let raw_feature = (
            b"ALA", b"ARG", 5.0_f32, -10.0_f32, 0.0_f32, 10.0_f32,
            345.0_f32, 15.0_f32,
        );
        let feature_input: Vec<f32> = vec![
            map_aa_to_u8(raw_feature.0) as f32, map_aa_to_u8(raw_feature.1) as f32,
            raw_feature.2, raw_feature.3.to_radians(), 
            raw_feature.4.to_radians(), raw_feature.5.to_radians(),
            raw_feature.6.to_radians(), raw_feature.7.to_radians(),
        ];
        let hash = GeometricHash::perfect_hash(feature_input, HashType::FoldDiscoDefault);
        assert_eq!(
            hash.downcast_fold_disco_default(),
            HashValue::from_u64(21050527393077_u64)
        ); //  NBIN_SIN_COS = 6

        let rev = hash.reverse_hash();
        assert_eq!(rev[0], 0.0);
        assert_eq!(rev[1], 1.0);
    }

}