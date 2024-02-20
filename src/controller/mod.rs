// File: mod.rs
// Created: 2024-01-18 15:47:16
// Description:
//    a new controller implementation that supports multiple hash types
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

pub mod feature;
pub mod io;
pub mod query;

use std::io::Write;
// External imports
use rayon::prelude::*;

// Internal imports
use crate::prelude::print_log_msg;
use crate::{measure_time, PDBReader};
use crate::geometry::core::{GeometricHash, HashType};
use crate::index::new_alloc::IndexBuilder;
use crate::utils::log::{ log_msg, INFO, FAIL, WARN, DONE };
use crate::controller::feature::get_geometric_hash_from_structure;

const DEFAULT_REMOVE_REDUNDANCY: bool = true;
const DEFAULT_NUM_THREADS: usize = 4;
const DEFAULT_HASH_TYPE: HashType = HashType::FoldDiscoDefault;

pub struct FoldDisco {
    pub path_vec: Vec<String>,
    pub numeric_id_vec: Vec<usize>,
    pub hash_collection: Vec<Vec<GeometricHash>>,
    pub index_builder: IndexBuilder<usize, GeometricHash>,
    pub hash_type: HashType,
    pub remove_redundancy: bool,
    pub num_threads: usize,
    pub output_path: String,
    pub flags: SubProcessFlags,
}

pub struct SubProcessFlags {
    pub fill_numeric_id_vec: bool,
    pub collect_hash: bool,
    pub set_index_table: bool,
    pub fill_index_table: bool,
    pub convert_index: bool,
    pub save_index_table: bool,
}
impl SubProcessFlags {
    pub fn new() -> SubProcessFlags {
        SubProcessFlags {
            fill_numeric_id_vec: false,
            collect_hash: false,
            set_index_table: false,
            fill_index_table: false,
            convert_index: false,
            save_index_table: false,
        }
    }
}

impl FoldDisco {
    pub fn new(path_vec: Vec<String>) -> FoldDisco {
        FoldDisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::new(),
            hash_collection: Vec::new(),
            index_builder: IndexBuilder::empty(),
            hash_type: DEFAULT_HASH_TYPE,
            remove_redundancy: DEFAULT_REMOVE_REDUNDANCY,
            num_threads: DEFAULT_NUM_THREADS,
            output_path: String::new(),
            flags: SubProcessFlags::new(),
        }
    }
    pub fn new_with_hash_type(path_vec: Vec<String>, hash_type: HashType) -> FoldDisco {
        FoldDisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::new(),
            hash_collection: Vec::new(),
            index_builder: IndexBuilder::empty(),
            hash_type: hash_type,
            remove_redundancy: DEFAULT_REMOVE_REDUNDANCY,
            num_threads: DEFAULT_NUM_THREADS,
            output_path: String::new(),
            flags: SubProcessFlags::new(),
        }
    }

    pub fn new_with_params(
        path_vec: Vec<String>, hash_type: HashType, remove_redundancy: bool,
        num_threads: usize, output_path: String
    ) -> FoldDisco {
        FoldDisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::new(),
            hash_collection: Vec::new(),
            index_builder: IndexBuilder::empty(),
            hash_type: hash_type,
            remove_redundancy: remove_redundancy,
            num_threads: num_threads,
            output_path: output_path,
            flags: SubProcessFlags::new(),
        }
    }
    // Setters
    pub fn set_path_vec(&mut self, path_vec: Vec<String>) {
        self.path_vec = path_vec;
    }
    pub fn set_hash_type(&mut self, hash_type: HashType) {
        self.hash_type = hash_type;
    }
    pub fn set_remove_redundancy(&mut self, remove_redundancy: bool) {
        self.remove_redundancy = remove_redundancy;
    }
    pub fn set_num_threads(&mut self, num_threads: usize) {
        self.num_threads = num_threads;
    }
    pub fn set_output_path(&mut self, output_path: String) {
        self.output_path = output_path;
    }

    // Main methods
    pub fn fill_numeric_id_vec(&mut self) {
        string_vec_to_numeric_id_vec(&self.path_vec, &mut self.numeric_id_vec);
        self.flags.fill_numeric_id_vec = true;
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
                        &compact.to_compact(), self.hash_type
                    );
                    // Drop intermediate variables
                    drop(compact);
                    drop(pdb_reader);
                    // If remove_redundancy is true, remove duplicates
                    if self.remove_redundancy {
                        let mut hash_vec = hash_vec;
                        hash_vec.sort_unstable();
                        hash_vec.dedup();
                        hash_vec
                    } else {
                        hash_vec
                    }
                })
            })
            .collect::<Vec<Vec<GeometricHash>>>();
        
        self.hash_collection = output;
        self.flags.collect_hash = true;
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
        // Check if hash_collection is filled
        if !self.flags.collect_hash {
            print_log_msg(FAIL, "Hash collection is not filled yet");
            return;
        }
        if !self.flags.fill_numeric_id_vec {
            print_log_msg(WARN, "Numeric ID vector is not filled yet. Filling numeric ID vector...");
            self.fill_numeric_id_vec();
        }
        let index_builder = IndexBuilder::new(
            &self.numeric_id_vec, &self.hash_collection,
            self.num_threads, self.get_allocation_size(),
            format!("{}.offset", self.output_path), // offset file path
            format!("{}.index", self.output_path), // data file path
        );
        self.index_builder = index_builder;
        self.flags.set_index_table = true;
    }
    
    pub fn fill_index_table(&mut self) {
        // Check if index_table is set
        if !self.flags.set_index_table {
            print_log_msg(WARN, "Index table is not set yet. Setting index table...");
            self.set_index_table();
        }
        // self.index_builder.fill_with_dashmap();
        let index_map = measure_time!(self.index_builder.fill_and_return_dashmap());
        self.flags.fill_index_table = true;

        // THIS IS NOT FINISHED YET
        // TODO: IMPORTANT: Move this.
        let (offset_map, value_vec) = measure_time!(self.index_builder.convert_hashmap_to_offset_and_values(index_map));
        self.flags.convert_index = true;
        
        // self.index_builder.offset_table = offset_map;
        // println!("offset_map: {:?}", offset_map);
        println!("value_vec: {:?}", value_vec.len());
    }
    
    
    pub fn save_offset_map(&self) {
        // Check if index_table is filled
        if !self.flags.fill_index_table {
            print_log_msg(FAIL, "Index table is not filled yet");
            return;
        }
        todo!("Implement save_offset_map method");
        // self.index_builder.save_offset_map();
    }
    
    pub fn save_index_table(&self) {
        // self.index_builder.save();
        todo!("Implement save_index_table method");
    }
    
}

fn string_vec_to_numeric_id_vec(string_vec: &Vec<String>, numeric_id_vec: &mut Vec<usize>) {
    for i in 0..string_vec.len() {
        numeric_id_vec.push(i);
    }
}

