// Point Pair Features (PPF):

use std::fmt;
use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u32);

impl HashValue {
    pub fn from_u32(hashvalue: u32) -> Self {
        HashValue(hashvalue)
    }

    pub fn perfect_hash(aa1: u8, aa2: u8, i: usize, j: usize, dist: f32) -> Self {
        let h1 = aa1 as u32;
        let h2 = aa2 as u32;
        assert!(i != j);
        let h3 = ((i as i32 - j as i32).abs() + 1).ilog2();
        let h4 = discretize_value(dist, 2.0, 20.0, 4.0);

        assert!(h1 < 256);
        assert!(h2 < 256);
        assert!(h3 < 256);
        assert!(h4 < 256);
        let hashvalue = h1 << 24 | h2 << 16 | h3 << 8 | h4;
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> [f32; 4] {
        let h1 = (self.0 >> 24) as u8;
        let h2 = (self.0 >> 16) as u8;
        let h3 = (self.0 >> 8) as u8;
        let h4 = (self.0 & 0x000000FF) as u8;
        [h1 as f32, h2 as f32, h3 as f32, h4 as f32]
    }

    pub fn neighbors(&self) -> Vec<HashValue> {
        let mut neighbors = Vec::new();
        let mut hash = self.0;
        // Append hashes with 1 bit difference for h3 & h4
        // h3-1, h3, h3+1; h4-1, h4, h4+1 (if possible) total 9 hashes
        // Explain code with comments
        for i in 0..9 {
            let mut h3 = (hash >> 8) as u8;
            let mut h4 = (hash & 0x000000FF) as u8;
            // h3-1, h3, h3+1; h4-1, h4, h4+1 (if possible) total 9 hashes
            // Explain code with comments
            if i < 3 {
                h3 = h3.wrapping_sub(1);
            } else if i > 5 {
                h3 = h3.wrapping_add(1);
            }
            if i % 3 == 0 {
                h4 = h4.wrapping_sub(1);
            } else if i % 3 == 2 {
                h4 = h4.wrapping_add(1);
            }
            // Check if h3 and h4 are valid

            let hashvalue = (h3 as u32) << 8 | h4 as u32;
            neighbors.push(HashValue(hashvalue));

        }
        neighbors
    }
    
}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let feature = self.reverse_hash();
        write!(
            f,
            "HashValue({}), aa1={}, aa2={}, logdist={}, dist={}",
            self.0, feature[0], feature[1], feature[2], feature[3]
        )
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let feature  = self.reverse_hash();
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.0, feature[0], feature[1], feature[2], feature[3]
        )
    }
}

pub fn discretize_value(val: f32, min: f32, max: f32, num_bin: f32) -> u32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    (val * (disc_f) + 0.5) as u32
}


pub type HashCollection = Vec<HashValue>;

pub fn map_aa_to_u8(aa: &[u8; 3]) -> u8 {
    match aa {
        b"ALA" => 0,
        b"ARG" => 1,
        b"ASN" => 2,
        b"ASP" => 3,
        b"CYS" => 4,
        b"GLN" => 5,
        b"GLU" => 6,
        b"GLY" => 7,
        b"HIS" => 8,
        b"ILE" => 9,
        b"LEU" => 10,
        b"LYS" => 11,
        b"MET" => 12,
        b"PHE" => 13,
        b"PRO" => 14,
        b"SER" => 15,
        b"THR" => 16,
        b"TRP" => 17,
        b"TYR" => 18,
        b"VAL" => 19,
        _ => panic!("Invalid AA"),
    }
}
pub fn map_u8_to_aa(aa: u8) -> &'static str {
    match aa {
        0 => "ALA",
        1 => "ARG",
        2 => "ASN",
        3 => "ASP",
        4 => "CYS",
        5 => "GLN",
        6 => "GLU",
        7 => "GLY",
        8 => "HIS",
        9 => "ILE",
        10 => "LEU",
        11 => "LYS",
        12 => "MET",
        13 => "PHE",
        14 => "PRO",
        15 => "SER",
        16 => "THR",
        17 => "TRP",
        18 => "TYR",
        19 => "VAL",
        _ => panic!("Invalid AA"),
    }
}

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

        // let h_cb_dist = discretize_value(cbd, 2.0, 20.0, 8.0); 3
        // // Torsion angles
        // let h_omega = discretize_value(omega, -1.0, 1.0, 4.0);
        // let h_theta1 = discretize_value(theta1, -1.0, 1.0, 4.0);
        // let h_theta2 = discretize_value(theta2, -1.0, 1.0, 4.0);
        // // Planar angles
        // let h_phi1 = discretize_value(phi1, 0.0, 180.0, 8.0);
        // let h_phi2 = discretize_value(phi2, 0.0, 180.0, 8.0);

use std::fmt;
use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u64);

impl HashValue {
    pub fn from_u64(hashvalue: u64) -> Self {
        HashValue(hashvalue)
    }
    pub fn as_u64(&self) -> u64 {
        self.0
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
        let h_cb_dist = discretize_value(cbd, 2.0, 20.0, 8.0);
        // Torsion angles
        let h_omega = discretize_value(omega, -1.0, 1.0, 4.0);
        let h_theta1 = discretize_value(theta1, -1.0, 1.0, 4.0);
        let h_theta2 = discretize_value(theta2, -1.0, 1.0, 4.0);
        // Planar angles
        let h_phi1 = discretize_value(phi1, 0.0, 180.0, 8.0);
        let h_phi2 = discretize_value(phi2, 0.0, 180.0, 8.0);

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

use std::fmt;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u16);

impl HashValue {
    pub fn perfect_hash(edge1: f32, edge2: f32, edge3: f32) -> Self {
        let mut hashvalue = 0;
        let edges = vec![edge1, edge2, edge3];
        // If edges are not sorted, panic
        if edges[0] > edges[1] || edges[0] > edges[2] || edges[1] > edges[2] {
            panic!("Edges are not sorted"); // Not panic, but return None
        }
        // If edges are not in range, panic
        for edge in edges.iter() {
            if *edge < 3.5 || *edge > 19.5 {
                panic!("Edges are not in range");
            }
        }

        let edge1 = (edge1.round() - 4.0_f32) as u16;
        // Assert that higher bits are empty
        let edge2 = (edge2.round() - 4.0_f32) as u16;
        let edge3 = (edge3.round() - 4.0_f32) as u16;
        hashvalue += edge1 << 8;
        hashvalue += edge2 << 4;
        hashvalue += edge3;
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> (u16, u16, u16) {
        let edge1 = ((self.0 >> 8) & 0b1111) + 4;
        let edge2 = ((self.0 >> 4) & 0b1111) + 4;
        let edge3 = (self.0 & 0b1111) + 4;
        if edge1 < 4 || edge1 > 19 {
            println!("{}", edge1);
            panic!("06 Edges are not in range");
        } else if edge2 < 4 || edge2 > 19 {
            panic!("07 Edges are not in range");
        } else if edge3 < 4 || edge3 > 19 {
            println!("{}", edge3);
            panic!("08 Edges are not in range");
        }
        (edge1, edge2, edge3)
    }
}

#[cfg(test)]
mod temp_test {
    #[test]
    fn test_one_hash() {
        let hash = super::HashValue::perfect_hash(5.0, 8.0, 11.0);
        assert_eq!(hash.0, 0b0000000101000111);
        let reverse = hash.reverse_hash();
        assert_eq!(reverse, (5, 8, 11));
    }
}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (edge1, edge2, edge3) = self.reverse_hash();
        write!(
            f,
            "HashValue({}), edge1={}, edge2={}, edge3={}",
            self.0, edge1, edge2, edge3
        )
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (edge1, edge2, edge3) = self.reverse_hash();
        write!(f, "{}\t{}\t{}\t{}", self.0, edge1, edge2, edge3)
    }
}

pub type HashCollection = Vec<HashValue>;

pub struct PossibleTriangles {
    pub dist_vec: Vec<u8>,
    pub triangles: Vec<(u8, u8, u8)>,
}

impl PossibleTriangles {
    pub fn new(dist_vec: Vec<u8>) -> Self {
        let mut triangles = Vec::new();
        for i in 0..dist_vec.len() {
            for j in i + 1..dist_vec.len() {
                for k in j + 1..dist_vec.len() {
                    if dist_vec[i] + dist_vec[j] > dist_vec[k]
                        && dist_vec[i] + dist_vec[k] > dist_vec[j]
                        && dist_vec[j] + dist_vec[k] > dist_vec[i]
                    {
                        triangles.push((dist_vec[i], dist_vec[j], dist_vec[k]));
                    }
                }
            }
        }
        PossibleTriangles {
            dist_vec,
            triangles,
        }
    }
    pub fn get_triangles(&self) -> &Vec<(u8, u8, u8)> {
        &self.triangles
    }
    pub fn get_dist_vec(&self) -> &Vec<u8> {
        &self.dist_vec
    }

    pub fn all_triangles(&self) -> Vec<(u8, u8, u8)> {
        let mut triangles = Vec::new();
        for i in 0..self.dist_vec.len() {
            for j in i + 1..self.dist_vec.len() {
                for k in j + 1..self.dist_vec.len() {
                    if self.dist_vec[i] + self.dist_vec[j] > self.dist_vec[k]
                        && self.dist_vec[i] + self.dist_vec[k] > self.dist_vec[j]
                        && self.dist_vec[j] + self.dist_vec[k] > self.dist_vec[i]
                    {
                        triangles.push((self.dist_vec[i], self.dist_vec[j], self.dist_vec[k]));
                    }
                }
            }
        }
        triangles
    }
}


// Point Pair Features (PPF):

use std::fmt;
use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u32);

impl HashValue {
    pub fn from_u32(hashvalue: u32) -> Self {
        HashValue(hashvalue)
    }

    pub fn perfect_hash(aa1: u8, aa2: u8, i: usize, j: usize, dist: f32) -> Self {
        let h1 = aa1 as u32;
        let h2 = aa2 as u32;
        assert!(i != j);
        let h3 = ((i as i32 - j as i32).abs() + 1).ilog2();
        let h4 = discretize_value(dist, 2.0, 20.0, 4.0);

        assert!(h1 < 256);
        assert!(h2 < 256);
        assert!(h3 < 256);
        assert!(h4 < 256);
        let hashvalue = h1 << 24 | h2 << 16 | h3 << 8 | h4;
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> [f32; 4] {
        let h1 = (self.0 >> 24) as u8;
        let h2 = (self.0 >> 16) as u8;
        let h3 = (self.0 >> 8) as u8;
        let h4 = (self.0 & 0x000000FF) as u8;
        [h1 as f32, h2 as f32, h3 as f32, h4 as f32]
    }

    pub fn neighbors(&self) -> Vec<HashValue> {
        let mut neighbors = Vec::new();
        let mut hash = self.0;
        // Append hashes with 1 bit difference for h3 & h4
        // h3-1, h3, h3+1; h4-1, h4, h4+1 (if possible) total 9 hashes
        // Explain code with comments
        for i in 0..9 {
            let mut h3 = (hash >> 8) as u8;
            let mut h4 = (hash & 0x000000FF) as u8;
            // h3-1, h3, h3+1; h4-1, h4, h4+1 (if possible) total 9 hashes
            // Explain code with comments
            if i < 3 {
                h3 = h3.wrapping_sub(1);
            } else if i > 5 {
                h3 = h3.wrapping_add(1);
            }
            if i % 3 == 0 {
                h4 = h4.wrapping_sub(1);
            } else if i % 3 == 2 {
                h4 = h4.wrapping_add(1);
            }
            // Check if h3 and h4 are valid

            let hashvalue = (h3 as u32) << 8 | h4 as u32;
            neighbors.push(HashValue(hashvalue));

        }
        neighbors
    }
    
}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let feature = self.reverse_hash();
        write!(
            f,
            "HashValue({}), aa1={}, aa2={}, logdist={}, dist={}",
            self.0, feature[0], feature[1], feature[2], feature[3]
        )
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let feature  = self.reverse_hash();
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.0, feature[0], feature[1], feature[2], feature[3]
        )
    }
}

pub fn discretize_value(val: f32, min: f32, max: f32, num_bin: f32) -> u32 {
    let cont_f = (max - min) / (num_bin - 1.0_f32);
    let disc_f = 1.0_f32 / cont_f;
    (val * (disc_f) + 0.5) as u32
}


pub type HashCollection = Vec<HashValue>;

pub fn map_aa_to_u8(aa: &[u8; 3]) -> u8 {
    match aa {
        b"ALA" => 0,
        b"ARG" => 1,
        b"ASN" => 2,
        b"ASP" => 3,
        b"CYS" => 4,
        b"GLN" => 5,
        b"GLU" => 6,
        b"GLY" => 7,
        b"HIS" => 8,
        b"ILE" => 9,
        b"LEU" => 10,
        b"LYS" => 11,
        b"MET" => 12,
        b"PHE" => 13,
        b"PRO" => 14,
        b"SER" => 15,
        b"THR" => 16,
        b"TRP" => 17,
        b"TYR" => 18,
        b"VAL" => 19,
        _ => panic!("Invalid AA"),
    }
}
pub fn map_u8_to_aa(aa: u8) -> &'static str {
    match aa {
        0 => "ALA",
        1 => "ARG",
        2 => "ASN",
        3 => "ASP",
        4 => "CYS",
        5 => "GLN",
        6 => "GLU",
        7 => "GLY",
        8 => "HIS",
        9 => "ILE",
        10 => "LEU",
        11 => "LYS",
        12 => "MET",
        13 => "PHE",
        14 => "PRO",
        15 => "SER",
        16 => "THR",
        17 => "TRP",
        18 => "TYR",
        19 => "VAL",
        _ => panic!("Invalid AA"),
    }
}







// -------------EXAMPLES---------------- //
// example/vae_test/01.rs
/*
 * File: 01.rs
 * Project: vae_test
 * Created: 2023-06-29 14:37:33
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code is written as part of project "vae_test".
 * ---
 * Last Modified: 2023-07-11 20:10:13
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2023 Hyunbin Kim, All rights reserved
 */

use std::fmt;
use tch::{
    nn::Module,
    nn::OptimizerConfig,
    nn::{self, VarStore},
    IndexOp, Kind, Reduction, Tensor,
};

use motifsearch::{geometry::two_float::{HashValue, HashCollection}, index::IndexTablePrinter};
use motifsearch::controller::{self, Controller, GeometryHashCollector};
use motifsearch::index::builder::IndexBuilder;
use motifsearch::PDBReader;

// Use hashmap
use std::collections::HashMap;

struct Encoder {
    fc1: nn::Linear,
    fc2: nn::Linear,
    fc3: nn::Linear,
}

impl Encoder {
    fn forward(&self, x: &Tensor) -> Tensor {
        x.apply(&self.fc1).relu().apply(&self.fc2).relu().apply(&self.fc3)
    }
}

fn get_tensor_from_vec(vec: &Vec<(String, Tensor)>, name: &str) -> Tensor {
    for (n, t) in vec {
        if n == name {
            let mut out = Tensor::empty(&t.size(), (Kind::Float, t.device()));
            out.copy_(&t);
            return out;
        }
    }
    panic!("No tensor named {}", name);
}

fn read_and_build_model() -> Encoder {
    let path = "data/encoder.safetensors";
    let loaded = tch::Tensor::read_safetensors(path).expect("Failed to read safetensors");

    let mut vs = VarStore::new(tch::Device::Cpu);
    let mut fc1 = nn::linear(&vs.root() / "fc1", 7, 64, Default::default());
    let mut fc2 = nn::linear(&vs.root() / "fc2", 64, 64, Default::default());
    let mut fc3 = nn::linear(&vs.root() / "fc3", 64, 2, Default::default());
    
    let fc1_weight = get_tensor_from_vec(&loaded, "0.weight");
    let fc1_bias = get_tensor_from_vec(&loaded, "0.bias");
    let fc2_weight = get_tensor_from_vec(&loaded, "2.weight");
    let fc2_bias = get_tensor_from_vec(&loaded, "2.bias");
    let fc3_weight = get_tensor_from_vec(&loaded, "4.weight");
    let fc3_bias = get_tensor_from_vec(&loaded, "4.bias");

    fc1.ws = fc1_weight;
    fc1.bs = Some(fc1_bias);
    fc2.ws = fc2_weight;
    fc2.bs = Some(fc2_bias);
    fc3.ws = fc3_weight;
    fc3.bs = Some(fc3_bias);
    
    Encoder {
        fc1: fc1,
        fc2: fc2,
        fc3: fc3,
    }
}

fn load_homeobox_toy() -> Vec<String> {
    vec![
        "data/homeobox/1akha-.pdb".to_string(),
        "data/homeobox/1b72a-.pdb".to_string(),
        "data/homeobox/1b72b-.pdb".to_string(),
        "data/homeobox/1ba5--.pdb".to_string(),
    ]
}
fn load_path(dir: &str) -> Vec<String> {
    // Load all pdbs in given path
    let mut pdb_paths = Vec::new();
    let paths = std::fs::read_dir(dir).expect("Unable to read pdb directory");
    for path in paths {
        let path = path.expect("Unable to read path");
        let path = path.path();
        let path = path.to_str().expect("Unable to convert path to string");
        // If the path is a pdb file, add it to the list
        if path.ends_with(".pdb") {
            pdb_paths.push(path.to_string());
        }
    }
    pdb_paths
}

fn main() {
    // IMPORTANT: Model should be saved as safetensors
    // Load model
    let path = "data/encoder.safetensors";
    let loaded = tch::Tensor::read_safetensors(path);
    let enc = read_and_build_model();
    // Test if encoder works
    // let temp = enc.forward(&Tensor::ones(&[1, 7], (Kind::Float, tch::Device::Cpu)));
    // println!("temp: {:?}", temp);

    // Load dataset
    // let dataset = load_homeobox_toy();
    let dataset = load_path("data/serine_peptidases_filtered");
    let mut controller = Controller::new(dataset.clone());
    controller.fill_numeric_id_vec();
    controller.path_vec = dataset.clone();
    for pdb_path in dataset {
        // Start measure time
        let start = std::time::Instant::now();
        let pdb_reader = PDBReader::from_file(&pdb_path).expect("Failed to read PDB file");
        let structure = pdb_reader.read_structure().expect("Failed to read structure");
        let compact = structure.to_compact();

        let mut hash_collector = GeometryHashCollector::new();

        let N: usize = compact.num_residues.try_into().unwrap();
        let mut trr_input_vec = Vec::with_capacity(N * N - N);
        let mut counter: usize = 0;
        let mut inner_res_pair_vec =  Vec::with_capacity(N * N - N);
        for i in 0..compact.num_residues {
            for j in 0..compact.num_residues {
                counter += 1;
                if i == j {
                    continue;
                }
                let trr = compact
                    .get_trrosetta_feature2(i, j);
                if trr.is_none() {
                    let empty: [f32; 7] = [0.0; 7];
                    let trr = empty;
                    trr_input_vec.push(trr);
                    // Fill inner_res_pair_vec with (i, j)
                    let inner_res_pair = (
                        compact.residue_serial[i], compact.residue_serial[j],
                        compact.residue_name[i], compact.residue_name[j],
                    );
                    inner_res_pair_vec.push(inner_res_pair);
                    counter += 1;
                    continue;
                }
                // Fill trr_input_vec with trr
                trr_input_vec.push(trr.unwrap());
                
                // Fill inner_res_pair_vec with (i, j)
                let inner_res_pair = (
                    compact.residue_serial[i], compact.residue_serial[j],
                    compact.residue_name[i], compact.residue_name[j],
                );
                inner_res_pair_vec.push(inner_res_pair);
                // dbg!("{},{}, {:?}, {}", i, j, trr, trr_input_vec.len());
                // print!("{}/{} ", trr, encoded);
                    // let hash_value =
                //     HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);
                // hash_collector.collect_hash(hash_value);
                counter += 1;
            }
        }
        // Convert trr_input_vec to tensor
        // Print trr_input_vec is filled
        // println!("trr_input_vec: {:?}", trr_input_vec);
        // Flatten
        let trr_input_vec = trr_input_vec.into_iter().flatten().collect::<Vec<_>>();
        let trr_input = Tensor::of_slice(&trr_input_vec).reshape(&[(N*N - N) as i64, 7]);
        let encoded = enc.forward(&trr_input);
        // encoded.print();
        // Convert torch tensor to vector
        let encoded_vec: Vec<Vec<f32>> = TryFrom::try_from(encoded).unwrap();
        // 2023-07-11 17:36:44
        let mut hash_vec = Vec::with_capacity(encoded_vec.len());
        for i in 0..encoded_vec.len() {
            let hash_value = HashValue::perfect_hash(encoded_vec[i][0], encoded_vec[i][1]);
            hash_vec.push(hash_value);
        }
        controller.hash_new_collection_vec.push(hash_vec);
        controller.res_pair_vec.push(inner_res_pair_vec);
    }
    println!("Numeric id length {:?}", controller.numeric_id_vec.len());
    println!("Hashvec length {:?}", controller.hash_new_collection_vec.len());
    let index_builder = IndexBuilder::new();
    let index_table = index_builder.concat(
        &controller.numeric_id_vec,
        &controller.hash_new_collection_vec,
    );
    let table_printer = IndexTablePrinter::Debug;
    table_printer.print(&index_table, "data/serine_vae_hash_table_rounded.tsv");
    println!("Checking");
    let mut serine_filter: HashMap<String, Vec<u64>> = HashMap::new();
    serine_filter.insert("1aq2.pdb".to_string(), vec![250, 232, 269]);
    serine_filter.insert("1wab.pdb".to_string(), vec![47, 195, 192]);
    serine_filter.insert("1sc9.pdb".to_string(), vec![80, 235, 207]);
    serine_filter.insert("2o7r.pdb".to_string(), vec![169, 306, 276]);
    serine_filter.insert("1bs9.pdb".to_string(), vec![90, 187, 175]);
    serine_filter.insert("1ju3.pdb".to_string(), vec![117, 287, 259]);
    serine_filter.insert("1uk7.pdb".to_string(), vec![34, 252, 224]);
    serine_filter.insert("1okg.pdb".to_string(), vec![255, 75, 61]);
    serine_filter.insert("1qfm.pdb".to_string(), vec![554, 680, 641]);
    controller.save_hash_per_pair("data/serine_vae_per_pair_hash.tsv");
    controller.save_filtered_hash_pair("data/serine_vae_per_pair_hash_filtered.tsv", &serine_filter);
    println!("DONE");
}


// example/index_io_old.rs

use std::{sync::{atomic::{AtomicUsize, Ordering}, Arc}, time::Instant, io::{BufWriter, Write}, fs::File};

use byteorder::LittleEndian;
// This example is to check if the IndexAllocator works as intended
use motifsearch::{prelude::*, structure::atom::Atom, index::alloc::resize_allocation};
use rayon::prelude::*;

fn main() {
    const NUM_THREADS: usize = 6;
    
    // Load directory
    let yeast_pdb_paths = motifsearch::utils::loader::load_path("data/serine_peptidases_filtered");
    let mut controller = motifsearch::controller::Controller::new(yeast_pdb_paths);
    controller.num_threads_file = NUM_THREADS;
    controller.num_threads_hash = NUM_THREADS;
    controller.fill_numeric_id_vec();
    
    let alloc_size = 
    
    let mut allocator = IndexAllocator::new(controller.num_threads_hash, alloc_size);
    let mut dedup_size: AtomicUsize = AtomicUsize::new(0);
    
    let start = Instant::now();
    let mut result = (&controller.path_vec).into_par_iter().map(
        |pdb_path| {
            // Read structure from PDB file
            let pdb_reader = PDBReader::from_file(&pdb_path).expect(
                &log_msg(FAIL, "Failed to read PDB file")
            );
            let compact = pdb_reader.read_structure().expect(
                &log_msg(FAIL, "Failed to read PDB file")
            );
            let compact = compact.to_compact();
            // Iterate over all residue pairs
            let res_bound = get_all_combination(compact.num_residues, false);
            // Single
            let mut hash_collection: Vec<usize> = res_bound.into_par_iter().enumerate().map(
                |(i, (res1, res2))| {
                        if res1 == res2 {
                            return 0usize;
                        }
                        let feature = compact.get_trrosetta_feature(res1, res2).unwrap_or(
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        );
                        if feature[0] < 2.0 || feature[0] > 20.0 {
                            return 0usize;
                        }
                        let hash = HashValue::perfect_hash(
                            feature[0], feature[1], feature[2], feature[3], feature[4], feature[5]
                        );
                        hash.as_usize()
                    }
            ).collect::<Vec<_>>();
            // Filter out invalid residue pairs & deduplicate
            // Remove zero
            hash_collection.sort_unstable();
            hash_collection.dedup();
            hash_collection.retain(|&x| x != 0);
            hash_collection.shrink_to_fit();
            
            // Check if reallocation is needed
            let len = hash_collection.len();
            let cursor = allocator.cursor.load(Ordering::Relaxed);
            let size = allocator.data_size.load(Ordering::Relaxed);
            if cursor + len >= size {
                // Use unsafe code to resize the allocation vector
                // Allow mutable method to be called in here
                unsafe {
                    allocator.data_size.store(size * 2, Ordering::Relaxed);
                    let size_as_mb = size as f32 * 8.0 / 1024.0 / 1024.0;
                    println!("Resized allocation vector to {}MB", size * 2);
                    // Get mutable reference to allocator with raw pointer
                    let mut ptr = allocator.allocation.as_ptr() as &mut [AtomicUsize];
                    // ptr.set_len(size * 2);
                    // Resize the allocation vector
                    resize_allocation(ptr, size * 2);
                    // Update the allocation vector
                }
            }
            
            // Allocate
            // allocator.fill_and_update(hash_collection);
            allocator.fill_usize_vec(hash_collection);
            
            // allocator.fill_usize_vec(hash_collection);
            
        }
    ).collect::<Vec<_>>();

    /** 
     * IMPLEMENTATION WITH VECTOR START
     *     let mut hash_collection: Vec<usize> = res_bound.into_par_iter().map(
     *           |(res1, res2)| {
     *                if res1 == res2 {
     *                    return 0;
     *                }
     *                let feature = compact.get_trrosetta_feature(res1, res2).unwrap_or(
     *                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
     *                );
     *                if feature[0] < 2.0 || feature[0] > 20.0 {
     *                    return 0;
     *                }
     *                let hash = HashValue::perfect_hash(
     *                    feature[0], feature[1], feature[2], feature[3], feature[4], feature[5]
     *                );
     *                hash.as_usize()
     *            }
     *        ).collect();
     *        
     *        // Filter out invalid residue pairs & deduplicate
     *        // Remove zero
     *        hash_collection.sort_unstable();
     *        hash_collection.dedup();
     *        hash_collection.retain(|&x| x != 0);
     *        // allocator.fill_and_update(hash_collection);
     *        // allocator.fill_usize_vec(hash_collection);
     *        // Drop everything
     *        drop(compact);
     *        drop(pdb_reader);
     *        hash_collection.shrink_to_fit();
     *        hash_collection
     *    }
     *).flatten().collect::<Vec<_>>();
     *result.shrink_to_fit();
     *let result = unsafe {
     *    std::slice::from_raw_parts(
     *        result.as_ptr() as *const u8,
     *        result.len() * std::mem::size_of::<usize>()
     *    )
     *};
     *   // Print size of result
     *   println!("Result size: {}", result.len());
     */

    // Resize the allocation
    allocator.save_to_file("analysis/raw_ecoli.bin").expect(
        &log_msg(FAIL, "Failed to save index file")
    );
    // // Save to file
    // let mut file = File::create("analysis/raw_ecoli.bin").expect(
    //     &log_msg(FAIL, "Failed to create index file")
    // );
    // let mut buf = BufWriter::new(file);
    // buf.write_all(result).expect(
    //     &log_msg(FAIL, "Failed to write index file")
    // );
    // buf.flush().expect(
    //     &log_msg(FAIL, "Failed to flush index file")
    // );
    
    let end = Instant::now();
    println!("Time elapsed with {} threads: {:?}", NUM_THREADS, end - start);
    // Clear everything
}

// example/kmeans.rs
// Generated
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

fn main() {
    // Read data
    let file = File::open("yeast_raw_feature_100000.tsv").unwrap();
    let reader = BufReader::new(file);
    let mut data = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        let mut row = Vec::new();
        for val in line.split("\t") {
            row.push(val.parse::<f32>().unwrap());
        }
        data.push(row);
    }

    // K means clustering
    let k = 1024;
    let mut centroids = Vec::new();
    for i in 0..k {
        centroids.push(data[i].clone());
    }
    loop {
        let mut clusters = vec![Vec::new(); k];
        for row in &data {
            let mut min_distance = std::f32::INFINITY;
            let mut closest_centroid = 0;
            for (i, centroid) in centroids.iter().enumerate() {
                let distance = euclidean_distance(row, centroid);
                if distance < min_distance {
                    min_distance = distance;
                    closest_centroid = i;
                }
            }
            clusters[closest_centroid].push(row.clone());
        }
        let mut new_centroids = Vec::new();
        for cluster in &clusters {
            let mut centroid = vec![0.0; 6];
            for row in cluster {
                for (i, val) in row.iter().enumerate() {
                    centroid[i] += val;
                }
            }
            for val in &mut centroid {
                *val /= cluster.len() as f32;
            }
            new_centroids.push(centroid);
        }
        if centroids == new_centroids {
            // Save clusters into file
            let mut file = File::create("yeast_clusters.tsv").unwrap();
            for (i, cluster) in clusters.iter().enumerate() {
                file.write(format!("{}\t", i).as_bytes()).unwrap();
                for row in cluster {
                    for val in row {
                        file.write(format!("{}\t", val).as_bytes()).unwrap();
                    }
                    file.write("\n".as_bytes()).unwrap();
                }
                file.write("\n".as_bytes()).unwrap();
            }
            break;
        }
        centroids = new_centroids;
    }

    // Save centroids into file
    let mut file = File::create("yeast_centroids.tsv").unwrap();
    for (i, centroid) in centroids.iter().enumerate() {
        file.write(format!("{}\t", i).as_bytes()).unwrap();
        for val in centroid {
            file.write(format!("{}\t", val).as_bytes()).unwrap();
        }
        file.write("\n".as_bytes()).unwrap();
    }
}

fn euclidean_distance(a: &[f32], b: &[f32]) -> f32 {
    let mut sum = 0.0;
    for i in 0..a.len() {
        sum += (a[i] - b[i]).powi(2);
    }
    sum.sqrt()
}

fn fast_assign_raw_val_to_cluster(raw_val: f32, centroids: &[Vec<f32>]) -> usize {
    let mut min_distance = std::f32::INFINITY;
    let mut closest_centroid = 0;
    for (i, centroid) in centroids.iter().enumerate() {
        let distance = (raw_val - centroid[0]).powi(2);
        if distance < min_distance {
            min_distance = distance;
            closest_centroid = i;
        }
    }
    closest_centroid
}

// example/train.rs
// use motifsearch::*;

// type EncoderLayer = (
//     (Linear<6, 32>, ReLU),
//     (Linear<32, 64>, ReLU),
//     (Linear<64, 6>, ReLU),
// );

// type DecoderLayer = (
//     (Linear<6, 64>, ReLU),
//     (Linear<64, 32>, ReLU),
//     (Linear<32, 6>, Sigmoid),
// );

use dfdx::optim::{Momentum, Sgd, SgdConfig, WeightDecay};
use dfdx::prelude::*;

fn main() {
    // Log start time
    let start = std::time::Instant::now();
    type EncoderLayer = (
        (Linear<6, 32>, BatchNorm1D<32>, ReLU),
        (Linear<32, 32>, BatchNorm1D<32>, ReLU),
        Linear<32, 2>,
    );
    type DecoderLayer = (
        (Linear<2, 32>, BatchNorm1D<32>, ReLU),
        (Linear<32, 32>, BatchNorm1D<32>, ReLU),
        Linear<32, 6>,
    );

    type Mlp = (EncoderLayer, DecoderLayer);

    let dev = AutoDevice::default();
    let mut mlp = dev.build_module::<Mlp, f32>();
    let mut grads = mlp.alloc_grads();
    let mut sgd = Sgd::new(
        &mlp,
        SgdConfig {
            lr: 1e-2,
            momentum: Some(Momentum::Nesterov(0.9)),
            weight_decay: Some(WeightDecay::L2(1e-4)),
        },
    );

    let x: Tensor<Rank2<4, 6>, f32, _> = dev.sample_normal();
    // let y: Tensor<Rank2<4, 2>, f32, _> = dev.sample_normal();

    for i in 0..50 {
        let prediction = mlp.forward_mut(x.trace(grads));
        let loss = mse_loss(prediction, x.clone());
        println!("Epoch {:03}: {:?}", i, loss.array());
        grads = loss.backward();
        sgd.update(&mut mlp, &grads).expect("SGD failed");
        mlp.zero_grads(&mut grads);
    }
    let end = std::time::Instant::now();
    println!("{:?}", mlp.forward(x).array());
    println!("Time elapsed: {:?}", end - start);
    mlp.save("test.npz").expect("Failed to save");
}

// example/vae.rs
// // VQ-VQE written with dfdx crate

// // Import

use dfdx::nn::modules::{BatchNorm1D, Linear, Module, ModuleVisitor, ReLU, TensorCollection};
use dfdx::optim::{Momentum, Sgd, SgdConfig, WeightDecay};
use dfdx::shapes::{Dtype, Rank1, Rank2};
use dfdx::tensor::{AutoDevice, SampleTensor, Tape, Tensor, Trace};
use dfdx::tensor_ops::Device;

//
struct VAE<
    const IN: usize,
    const HIDDEN1: usize,
    const HIDDEN2: usize,
    const OUT: usize,
    E: Dtype,
    D: Device<E>,
> {
    encoder: (
        (Linear<IN, HIDDEN1, E, D>, BatchNorm1D<HIDDEN1, E, D>, ReLU),
        (Linear<HIDDEN1, HIDDEN2, E, D>, BatchNorm1D<HIDDEN2, E, D>, ReLU),
    ),
    mean: Linear<HIDDEN2, OUT, E, D>,
    logvar: Linear<HIDDEN2, OUT, E, D>,
    decoder: (
        (Linear<OUT, HIDDEN2, E, D>, BatchNorm1D<HIDDEN2, E, D>, ReLU),
        ( Linear<HIDDEN2, HIDDEN1, E, D>, BatchNorm1D<HIDDEN1, E, D>, ReLU),
        Linear<HIDDEN1, IN, E, D>,
    ),
}

impl<
        const IN: usize,
        const HIDDEN1: usize,
        const HIDDEN2: usize,
        const OUT: usize,
        E: Dtype,
        D: Device<E>,
    > TensorCollection<E, D> for VAE<IN, HIDDEN1, HIDDEN2, OUT, E, D>
{
    type To<E2: Dtype, D2: Device<E2>> = VAE<IN, HIDDEN1, HIDDEN2, OUT, E2, D2>;

    fn iter_tensors<V: ModuleVisitor<Self, E, D>>(
            visitor: &mut V,
        ) -> Result<Option<Self::To<V::E2, V::D2>>, V::Err> {
        visitor.visit_fields(
            // Define name of each field and how to access it,
            // using ModuleField for Modules, and TensorField for Tensors
            (
                Self::module("encoder", |s| &s.encoder, |s| &mut s.encoder),
                Self::module("mean", |s| &s.mean, |s| &mut s.mean),
                Self::module("logvar", |s| &s.logvar, |s| &mut s.logvar),
                Self::module("decoder", |s| &s.decoder, |s| &mut s.decoder)
            ),
            |(encoder, mean, logvar, decoder)| VAE {
                encoder,
                mean,
                logvar,
                decoder,
            },
        )
    }

}
