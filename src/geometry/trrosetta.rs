// Geometric features from trRosetta paper
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

use std::fmt;
use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u64);

impl HashValue {
    pub fn from_u64(hashvalue: u64) -> Self {
        HashValue(hashvalue)
    }

    pub fn perfect_hash(
        cb_dist: f32,
        omega: f32,
        theta1: f32,
        theta2: f32,
        phi1: f32,
        phi2: f32,
    ) -> Self {
        let mut cbd = cb_dist;
        if cb_dist > 20.0 {
            cbd = 20.0;
        }
        let h_cb_dist = discretize_value(cbd, 2.0, 20.0, 10.0);
        // Torsion angles
        let h_omega = discretize_value(omega, -1.0, 1.0, 6.0);
        let h_theta1 = discretize_value(theta1, -1.0, 1.0, 6.0);
        let h_theta2 = discretize_value(theta2, -1.0, 1.0, 6.0);
        // Planar angles
        let h_phi1 = discretize_value(phi1, 0.0, 180.0, 6.0);
        let h_phi2 = discretize_value(phi2, 0.0, 180.0, 6.0);

        assert!(h_cb_dist < 256);
        assert!(h_omega < 256);
        assert!(h_theta1 < 256);
        assert!(h_theta2 < 256);
        assert!(h_phi1 < 256);
        assert!(h_phi2 < 256);
        let hashvalue = h_cb_dist << 40
            | h_omega << 32
            | h_theta1 << 24
            | h_theta2 << 16
            | h_phi1 << 8
            | h_phi2;
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> [f32; 6] {
        let h_cb_dist = (self.0 >> 40) as u8;
        let h_omega = (self.0 >> 32) as u8;
        let h_theta1 = (self.0 >> 24) as u8;
        let h_theta2 = (self.0 >> 16) as u8;
        let h_phi1 = (self.0 >> 8) as u8;
        let h_phi2 = (self.0 & 0x000000FF) as u8;
        [
            h_cb_dist as f32,
            h_omega as f32,
            h_theta1 as f32,
            h_theta2 as f32,
            h_phi1 as f32,
            h_phi2 as f32,
        ]
    }
}

// impl Hasher for HashValue {
//     fn finish(&self) -> u64 {
//         self.0
//     }

//     fn write(&mut self, _bytes: &[u8]) {
//         unimplemented!()
//     }
// }

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
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.0, val[0], val[1], val[2], val[3], val[4], val[5]
        )
        // write!(f, "{}", self.0)
    }
}

pub type HashCollection = Vec<HashValue>;

pub fn discretize_value(val: f32, min: f32, max: f32, num_bin: f32) -> u64 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    (val * (disc_f) + 0.5) as u64
}

pub fn continuize_value(val: u64, min: f32, max: f32, num_bin: f32) -> f32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    (val as f32) * (cont_f)
}

pub fn normalize_angle_degree(val: f32, min: f32, max: f32) -> f32 {
    (val - min) / (max - min)
}
