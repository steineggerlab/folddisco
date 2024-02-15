// File: new_mod.rs
// Created: 2024-01-18 15:47:16
// Author: Hyunbin Kim (khb7840@gmail.com)
// Description:
//    a new controller implementation that supports multiple hash types
// Copyright © 2024 Hyunbin Kim, All rights reserved

use std::io::Write;

// External imports
use dashmap::DashMap;
use rayon::prelude::*;

// Internal imports
use crate::PDBReader;
use crate::geometry::core::{GeometricHash, HashType};
use crate::index::new_alloc::{ IndexBuilder, HashableSync };
use crate::structure::core::CompactStructure;
use crate::utils::log::{ log_msg, INFO, FAIL, WARN, DONE };

pub struct FoldDisco<K: HashableSync> {
    pub path_vec: Vec<String>,
    pub numeric_id_vec: Vec<usize>,
    pub hash_collection: Vec<Vec<K>>,
    pub index_builder: IndexBuilder<K, usize>,
    pub remove_redundancy: bool,
    pub num_threads: usize,
}

impl<K: HashableSync> FoldDisco<K> {
    pub fn new(path_vec: Vec<String>) -> FoldDisco<K> {
        FoldDisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::new(),
            hash_collection: Vec::new(),
            index_builder: IndexBuilder::empty(),
            remove_redundancy: false,
            num_threads: 4,
        }
    }

    pub fn collect_hash(&mut self) {
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating files");
        // For iterating files, apply multi-threading with num_threads_for_file
        let output = pool.install(|| {
            self.path_vec
                .par_iter()
                .map(|pdb_path| {
                    let pdb_reader = PDBReader::from_file(pdb_path).expect(
                        log_msg(FAIL, "PDB file not found").as_str()
                    );
                    let compact = pdb_reader.read_structure().expect(
                        log_msg(FAIL, "Failed to read structure").as_str()
                    );
                    let hash_vec = get_geometric_hash_from_structure(
                        &compact, HashType::FoldDiscoDefault
                    );
                    // Drop intermediate variables
                    drop(compact);
                    drop(pdb_reader);
                    hash_vec
                })
            })
            .collect::<Vec<Vec<K>>>();
        self.hash_collection = output;
    }

    pub fn get_allocation_size(&self) -> usize {
        let mut allocation_size = 0usize;
        self.path_vec.iter().for_each(|pdb_path| {
            let pdb_reader = PDBReader::from_file(pdb_path).expect(
            log_msg(FAIL, "PDB file not found").as_str()
            );
            let compact = pdb_reader.read_structure().expect(
                log_msg(FAIL, "Failed to read structure").as_str()
            );
            allocation_size += compact.num_residues * (compact.num_residues - 1);
        });
        allocation_size
    }

    pub fn save_raw_feature(&mut self, path: &str, discretize: bool) {}

    pub fn fill_numeric_id_vec(&mut self) {
        string_vec_to_numeric_id_vec(&self.path_vec, &mut self.numeric_id_vec);
    }

    pub fn save_id_vec(&self, path: &str) {
        // Save numeric_id_vec & path_vec as headerless tsv
        let mut file = std::fs::File::create(path).expect("Unable to create file");
        for i in 0..self.numeric_id_vec.len() {
            let numeric_id = self.numeric_id_vec[i];
            let path = &self.path_vec[i];
            file.write_all(format!("{}\t{}\n", numeric_id, path).as_bytes())
                .expect("Unable to write data");
        }
    }
}

fn string_vec_to_numeric_id_vec(string_vec: &Vec<String>, numeric_id_vec: &mut Vec<usize>) {
    for i in 0..string_vec.len() {
        numeric_id_vec.push(i);
    }
}

fn get_all_combination(n: usize, include_same: bool) -> Vec<(usize, usize)> {
    let mut res = Vec::new();
    for i in 0..n {
        for j in 0..n {
            if i == j && !include_same {
                continue;
            }
            res.push((i, j));
        }
    }
    res
}

pub fn get_geometric_hash_from_structure(structure: &CompactStructure, hash_type: HashType) -> Vec<GeometricHash> {
    let mut hash_vec = Vec::new();

    let res_bound = get_all_combination(
        structure.num_residues, false
    );

    res_bound.iter().for_each(|(i, j)| {
        let feature =
        match hash_type {
            HashType::PDBMotif => {
                
                let dist = structure.get_distance(*i, *j);
                let angle = structure.get_angle(*i, *j);
                if dist.is_some() && angle.is_some() {
                    let hash = HashType::perfect_hash(
                        vec![dist.unwrap(), angle.unwrap()]
                    );
                    hash_vec.push(hash);
                }
                // hash_vec.push(hash);
                vec![0.0_f32, 0.0_f32]
            },
            HashType::FoldDiscoDefault => {
                let feature = structure.get_trrosetta_feature2(*i, *j);
                let hash = H::perfect_hash(feature);
                hash_vec.push(hash);
            },
            _ => {
                unimplemented!();
            }
        }
        
        
        let hash = GeometricHash::perfect_hash(feature, hash_type);
        hash_vec.push(hash);
    });
    

    hash_vec
}




//


#[cfg(test)]
mod tests {
    use super::*;

}