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
use crate::{PDBReader, HashableSync};
use crate::geometry::core::{GeometricHash, HashType};
use crate::geometry::util::map_aa_to_u8;
use crate::index::new_alloc::IndexBuilder;
use crate::structure::core::CompactStructure;
use crate::utils::log::{ log_msg, INFO, FAIL, WARN, DONE };

use super::io::convert_hashmap_to_offset_and_values;

pub struct FoldDisco {
    pub path_vec: Vec<String>,
    pub numeric_id_vec: Vec<usize>,
    pub hash_collection: Vec<Vec<GeometricHash>>,
    pub index_builder: IndexBuilder<usize, GeometricHash>,
    pub remove_redundancy: bool,
    pub num_threads: usize,
    pub output_path: String,
}

impl FoldDisco {
    pub fn new(path_vec: Vec<String>) -> FoldDisco {
        FoldDisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::new(),
            hash_collection: Vec::new(),
            index_builder: IndexBuilder::empty(),
            remove_redundancy: false,
            num_threads: 4,
            output_path: String::new(),
        }
    }
    pub fn set_output_path(&mut self, output_path: String) {
        self.output_path = output_path;
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
                        &compact.to_compact(), HashType::FoldDiscoDefault
                    );
                    // Drop intermediate variables
                    drop(compact);
                    drop(pdb_reader);
                    hash_vec
                })
            })
            .collect::<Vec<Vec<GeometricHash>>>();
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

    pub fn save_raw_feature(&mut self, path: &str, discretize: bool) {
        todo!("Implement save_raw_feature method!!!");
    }

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
    
    pub fn set_index_table(&mut self) {
        let mut index_builder = IndexBuilder::new(
            &self.numeric_id_vec, &self.hash_collection,
            self.num_threads, self.get_allocation_size(),
            format!("{}.offset", self.output_path), // offset file path
            format!("{}.index", self.output_path), // data file path
        );
        self.index_builder = index_builder;
    }
    
    // TODO: IMPORTANT:
    // pub fn fill_index_table(&mut self) {
    //     let index_map = self.index_builder.fill_and_return_dashmap()
    //     let (offset_map, value_vec) = convert_hashmap_to_offset_and_values(index_map);
    // }
    
    // pub fn save_index_table(&self) {
    //     self.index_builder.save();
    // }
    
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
    let res_bound = get_all_combination(
        structure.num_residues, false
    );
    let mut hash_vec = Vec::new();

    res_bound.iter().for_each(|(i, j)| {
        // Get residue name and convert to f32
        let res1 = structure.get_res_name(*i);
        let res2 = structure.get_res_name(*j);
        let res1 = map_aa_to_u8(res1) as f32;
        let res2 = map_aa_to_u8(res2) as f32;
        
        let feature = match &hash_type {
            HashType::PDBMotif => {
                let ca_dist = structure.get_ca_distance(*i, *j);
                let cb_dist = structure.get_cb_distance(*i, *j);
                let ca_cb_angle = structure.get_ca_cb_angle(*i, *j); // degree
                if ca_dist.is_some() && cb_dist.is_some() && ca_cb_angle.is_some() {
                    let feature = vec![
                        res1, res2, ca_dist.unwrap(), cb_dist.unwrap(), ca_cb_angle.unwrap()
                    ];
                    feature
                } else {
                    vec![0.0; 5]
                }
            },
            HashType::FoldDiscoDefault => {
                let feature = structure.get_default_feature(*i, *j);
                if feature.is_some() {
                    // Concatenate res1 and res2 to the feature
                    let mut feature = feature.unwrap();
                    feature.insert(0, res1);
                    feature.insert(1, res2);
                    feature
                } else {
                    vec![0.0; 8]
                }
            },
            _ => {
                todo!("Implement feature-collection methods for other hash types here");
            }
        };
        let hash = GeometricHash::perfect_hash(feature, hash_type);
        hash_vec.push(hash);
    });
    
    hash_vec
}

//

                        // // Filter out invalid residue pairs & deduplicate
                        // hash_collection.sort_unstable();
                        // hash_collection.dedup();
                        
                        // // drop unnecessary variables
                        // drop(compact);
                        // drop(pdb_reader);
                        // hash_collection.shrink_to_fit();
// let (offset_table, value_vec) = convert_hashmap_to_offset_and_values(index_table);

