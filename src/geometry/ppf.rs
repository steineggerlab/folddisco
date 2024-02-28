// Point Pair Features (PPF):

use std::fmt;
use crate::geometry::core::HashType;
use crate::geometry::util::discretize_f32_value_into_u32 as discretize_value;
use crate::geometry::util::continuize_u32_value_into_f32 as continuize_value;
use crate::geometry::util::map_u8_to_aa;

// Constants
// 1. for cb_dist
const MIN_DIST: f32 = 2.0;
const MAX_DIST: f32 = 20.0;
const NBIN_DIST: f32 = 16.0;
// 2. NEW IDEA for encoding angles; represent as sin and cos
const MIN_SIN_COS: f32 = -1.0;
const MAX_SIN_COS: f32 = 1.0;
const NBIN_SIN_COS: f32 = 8.0;
// Bitmasks
const BITMASK32_3BIT: u32 = 0x00000007;
const BITMASK32_4BIT: u32 = 0x0000000F;
const BITMASK32_5BIT: u32 = 0x0000001F;

// 5 bit for AA, 4 bit for distance, 3 bit for sin & cos

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u32);

impl HashValue {
    pub fn perfect_hash(ppf: Vec<f32>) -> Self {
        let aa1 = ppf[0] as u32;
        let aa2 = ppf[1] as u32;
        let h_dist = discretize_value(ppf[2], MIN_DIST, MAX_DIST, NBIN_DIST); 
        let sin_vec = vec![ppf[3].sin(), ppf[4].sin(), ppf[5].sin()];
        let cos_vec = vec![ppf[3].cos(), ppf[4].cos(), ppf[5].cos()];
        let h_sin_vec = sin_vec.iter().map(
            |&x| discretize_value(x, MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS)
        ).collect::<Vec<u32>>();
        let h_cos_vec = cos_vec.iter().map(
            |&x| discretize_value(x, MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS)
        ).collect::<Vec<u32>>();
        let hashvalue = (aa1 << 27) | (aa2 << 22) | (h_dist << 18)
            | (h_sin_vec[0] << 15) | (h_cos_vec[0] << 12) | (h_sin_vec[1] << 9)
            | (h_cos_vec[1] << 6) | (h_sin_vec[2] << 3) | (h_cos_vec[2]);
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> Vec<f32> {
        let aa1 = (self.0 >> 27) & BITMASK32_5BIT as u32;
        let aa2 = (self.0 >> 22) & BITMASK32_5BIT as u32;
        let dist = continuize_value(
            (self.0 >> 18) & BITMASK32_4BIT as u32, 
            MIN_DIST, MAX_DIST, NBIN_DIST
        );
        let sin_n1_d = continuize_value(
            (self.0 >> 15) & BITMASK32_3BIT as u32, 
            MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS
        );
        let cos_n1_d = continuize_value(
            (self.0 >> 12) & BITMASK32_3BIT as u32, 
            MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS
        );
        let sin_n2_d = continuize_value(
            (self.0 >> 9) & BITMASK32_3BIT as u32, 
            MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS
        );
        let cos_n2_d = continuize_value(
            (self.0 >> 6) & BITMASK32_3BIT as u32, 
            MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS
        );
        let sin_n3_d = continuize_value(
            (self.0 >> 3) & BITMASK32_3BIT as u32, 
            MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS
        );
        let cos_n3_d = continuize_value(
            self.0 & BITMASK32_3BIT as u32, 
            MIN_SIN_COS, MAX_SIN_COS, NBIN_SIN_COS
        );
        vec![
            aa1 as f32, aa2 as f32, dist, 
            sin_n1_d.atan2(cos_n1_d).to_degrees(),
            sin_n2_d.atan2(cos_n2_d).to_degrees(),
            sin_n3_d.atan2(cos_n3_d).to_degrees()
        ]
    }

    pub fn hash(&self) -> u32 {
        self.0
    }
    
    pub fn from_u32(hashvalue: u32) -> Self {
        HashValue(hashvalue)
    }

    pub fn from_u64(hashvalue: u64) -> Self {
        HashValue(hashvalue as u32)
    }
    
    pub fn as_u32(&self) -> u32 {
        self.0
    }
    
    pub fn as_u64(&self) -> u64 {
        self.0 as u64
    }
    
    pub fn hash_type(&self) -> HashType {
        HashType::PointPairFeature
    }

}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ppf = self.reverse_hash();
        write!(
            f,
            "HashValue({}), ({},{}) d={}, ({},{},{})",
            self.0, map_u8_to_aa(ppf[0] as u8),  map_u8_to_aa(ppf[1] as u8),
            ppf[2], ppf[3], ppf[4], ppf[5]
        )
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ppf = self.reverse_hash();
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.0, ppf[0], ppf[1], ppf[2], ppf[3], ppf[4], ppf[5]
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::util::map_aa_to_u8;
    use super::*;

    #[test]
    fn test_hashvalue() {
        let aa1 = b"ALA";
        let aa2 = b"GLY";
        let mut ppf = vec![
            7.5, 120.0_f32.to_radians(), 45.0_f32.to_radians(), -60.0_f32.to_radians()
        ];
        ppf.insert(0, map_aa_to_u8(aa1) as f32);
        ppf.insert(1, map_aa_to_u8(aa2) as f32);
        let hashvalue = HashValue::perfect_hash(ppf);
        println!("{:?}", hashvalue);
    }
}