// Default hasher using trrosetta features
// 32 bit implementation default hashing scheme

// Implementation
// 5bit aa1, 5bit aa2, cbeta distance - 16 bins; 4 bits
// Angles - 6 bins each for sin and cos; total 25 bins; 4 bits for sin and cos each
 
use std::fmt;
use crate::geometry::core::HashType;
use crate::geometry::util::discretize_f32_value_into_u32 as discretize_value;
use crate::geometry::util::continuize_u32_value_into_f32 as continuize_value;
use crate::geometry::util::*;


use super::util::map_u32_to_aa_u32_pair;

// TODO: IMPORTANT: Implement this to reduce execution time
// 9 bit for AA pairs, 3 bit for distance, 2 bits for sin & cos (4 bits for one angle)
// TOTAL: 9 + 3 + (2 * 2 * 5) = 32 bits

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(pub u32);

impl HashValue {
    pub fn perfect_hash(feature: Vec<f32>, nbin_dist: usize, nbin_angle: usize) -> Self {
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        HashValue::_perfect_hash(
            res1, res2, feature[2], feature[3],
            feature[4], feature[5], feature[6], feature[7],
            nbin_dist as f32, nbin_angle as f32
        )
    }
    pub fn perfect_hash_default(feature: Vec<f32>) -> Self {
        let res1 = feature[0] as u32;
        let res2 = feature[1] as u32;
        HashValue::_perfect_hash(
            res1, res2, feature[2], feature[3],
            feature[4], feature[5], feature[6], feature[7],
            NBIN_DIST as f32, NBIN_SIN_COS as f32
        )
    }
    pub fn reverse_hash(&self, nbin_dist:usize, nbin_angle:usize) -> Vec<f32> {
        self._reverse_hash(nbin_dist as f32, nbin_angle as f32).to_vec()
    }
    pub fn reverse_hash_default(&self) -> Vec<f32> {
        self._reverse_hash(NBIN_DIST, NBIN_SIN_COS).to_vec()
    }
    pub fn hash_type(&self) -> super::core::HashType {
        HashType::Default32bit
    }

    fn _perfect_hash(
        res1: u32, res2: u32, cb_dist: f32, omega: f32, theta1: f32, theta2: f32,
        phi1: f32, phi2: f32, nbin_dist: f32, nbin_angle: f32
    ) -> Self {
        // By default, bit for the distance is 3 and angle is 2
        let nbin_dist = if nbin_dist > 8.0 { 8.0 } else { nbin_dist };
        let nbin_angle = if nbin_angle > 4.0 { 4.0 } else { nbin_angle };
        
        let res_pair = map_aa_u32_pair_to_u32(res1, res2);
        let h_cb_dist = discretize_value(cb_dist, MIN_DIST, MAX_DIST, nbin_dist);
        
        // Convert angles to sin and cos
        let angles = [omega, theta1, theta2, phi1, phi2];
        let sin_cos_angles = angles.iter().map(
            |&x| (x.sin(), x.cos())
        ).collect::<Vec<(f32, f32)>>();
        // Discretize sin and cos
        let sin_cos_angles = sin_cos_angles.iter().enumerate().map(
            |(i, &(sin, cos))| {
                if i < 3 {
                    (
                        discretize_value(sin, MIN_SIN_COS, MAX_SIN_COS, nbin_angle),
                        discretize_value(cos, MIN_SIN_COS, MAX_SIN_COS, nbin_angle)
                    )
                } else {
                    (
                        discretize_value(sin, MIN_SIN_COS, MAX_SIN_COS, nbin_angle),
                        discretize_value(cos, MIN_SIN_COS, MAX_SIN_COS, nbin_angle)
                    )
                }
            }
        ).collect::<Vec<(u32, u32)>>();
        
        // Combine all the hash values
        let hashvalue = res_pair << 23 | h_cb_dist << 20
            | sin_cos_angles[0].0 << 18 | sin_cos_angles[0].1 << 16 // sin_omega, cos_omega
            | sin_cos_angles[1].0 << 14 | sin_cos_angles[1].1 << 12 // sin_theta1, cos_theta1
            | sin_cos_angles[2].0 << 10 | sin_cos_angles[2].1 << 8 // sin_theta2, cos_theta2
            | sin_cos_angles[3].0 << 6 | sin_cos_angles[3].1 << 4 // sin_phi1, cos_phi1
            | sin_cos_angles[4].0 << 2 | sin_cos_angles[4].1; // sin_phi2, cos_phi2

        HashValue(hashvalue)
    }

    fn _reverse_hash(&self, _nbin_dist: f32, nbin_angle: f32) -> [f32; 8] {
        let res_pair = ((self.0 >> 23) & BITMASK32_9BIT) as u32;
        let (res1, res2) = map_u32_to_aa_u32_pair(res_pair);
        // Mask bits
        let cb_dist = ((self.0 >> 20) & BITMASK32_3BIT) as f32;
        let sin_cos_vec = [
            ((self.0 >> 18) & BITMASK32_2BIT), // sin_omega
            ((self.0 >> 16) & BITMASK32_2BIT), // cos_omega
            ((self.0 >> 14) & BITMASK32_2BIT), // sin_theta1
            ((self.0 >> 12) & BITMASK32_2BIT), // cos_theta1
            ((self.0 >> 10) & BITMASK32_2BIT), // sin_theta2
            ((self.0 >> 8) & BITMASK32_2BIT), // cos_theta2
            ((self.0 >> 6) & BITMASK32_2BIT), // sin_phi1
            ((self.0 >> 4) & BITMASK32_2BIT), // cos_phi1
            ((self.0 >> 2) & BITMASK32_2BIT), // sin_phi2
            (self.0 & BITMASK32_2BIT), // cos_phi2
        ];

        let sin_cos_vec = sin_cos_vec.iter().enumerate().map(
            |(i, &x)| {
                if i < 6 {
                    continuize_value(
                        x, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
                    )
                } else {
                    continuize_value(
                        x, MIN_SIN_COS, MAX_SIN_COS, nbin_angle
                    )
                }
            }
        ).collect::<Vec<f32>>();
        
        // Restores original angles
        let omega = sin_cos_vec[0].atan2(sin_cos_vec[1]).to_degrees();
        let theta1 = sin_cos_vec[2].atan2(sin_cos_vec[3]).to_degrees();
        let theta2 = sin_cos_vec[4].atan2(sin_cos_vec[5]).to_degrees();
        let phi1 = sin_cos_vec[6].atan2(sin_cos_vec[7]).to_degrees();
        let phi2 = sin_cos_vec[8].atan2(sin_cos_vec[9]).to_degrees();
        [res1 as f32, res2 as f32, cb_dist, omega, theta1, theta2, phi1, phi2]
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
    pub fn as_usize(&self) -> usize {
        self.0 as usize
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
            b"ALA", b"ARG", 5.0_f32, -10.0_f32, 0.0_f32, 10.0_f32, 45.0_f32, 15.0_f32,
        );
        let feature_input: Vec<f32> = vec![
            map_aa_to_u8(raw_feature.0) as f32, map_aa_to_u8(raw_feature.1) as f32,
            raw_feature.2, raw_feature.3.to_radians(), 
            raw_feature.4.to_radians(), raw_feature.5.to_radians(),
            raw_feature.6.to_radians(), raw_feature.7.to_radians(),
        ];
        let hash = GeometricHash::perfect_hash_default(feature_input, HashType::Default32bit);
        //
        println!("{:?}", hash);

        let rev = hash.reverse_hash_default();
        assert_eq!(rev[0], 0.0);
        assert_eq!(rev[1], 1.0);
    }

    #[test]
    fn test_multiple_bins() {
        let bin_pairs = vec![(8, 4), (6, 3), (4, 2)];
        for (nbin_dist, nbin_angle) in bin_pairs {
            let raw_feature = (
                b"ALA", b"ARG", 10.0_f32, -10.0_f32, 0.0_f32, 10.0_f32, 45.0_f32, 15.0_f32,
            );
            let feature_input: Vec<f32> = vec![
                map_aa_to_u8(raw_feature.0) as f32, map_aa_to_u8(raw_feature.1) as f32,
                raw_feature.2, raw_feature.3.to_radians(), 
                raw_feature.4.to_radians(), raw_feature.5.to_radians(),
                raw_feature.6.to_radians(), raw_feature.7.to_radians(),
            ];
            let hash = GeometricHash::perfect_hash(
                feature_input, HashType::Default32bit, nbin_dist, nbin_angle
            );
            let rev = hash.reverse_hash(nbin_dist, nbin_angle);
            println!("{:?}", hash);
            println!("{:?}", rev);
        }
    }
    
}