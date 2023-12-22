// File: build_index.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for building index table
// For building index table, we need a directory containing PDB files, 
// and the path to save the index table.

use crate::cli::*;
use crate::index::index_table::{Value, self};
use crate::index::lookup::{save_lookup_to_file, load_lookup_from_file};
use crate::index::io::*;
use rayon::prelude::*;
// use std::collections::HashMap;
use dashmap::DashMap;
// Import unique
use std::collections::HashSet;
use std::hash;
use std::sync::Arc;

use crate::prelude::*;

// TODO: Generated. need to be

const HELP_INDEX: &str = "\
USAGE: motifsearch index [OPTIONS]
Options:
    -d, --pdb-dir <PDB_DIR>     Directory containing PDB files
    -i, --index-path <INDEX_PATH>   Path to save the index table
    -t, --threads <THREADS>     Number of threads to use
    -v, --verbose               Print verbose messages
    -h, --help                  Print this help menu
";

pub fn build_index(env: AppArgs) {
    match env {
        AppArgs::Index {
            pdb_dir,
            pdb_path_vec,
            index_path,
            num_threads,
            verbose,
            help,
        } => {
            // Check if arguments are valid
            if pdb_dir.is_none() && pdb_path_vec.is_empty() {
                eprintln!("{}", HELP_INDEX);
                std::process::exit(1);
            }
            if help {
                eprintln!("{}", HELP_INDEX);
                std::process::exit(0);
            } else {
                if verbose { print_log_msg(INFO, "Building index table..."); }
                // Setup multithreading
                let start = std::time::Instant::now();
                rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global().unwrap();
                // Load PDB files
                let pdb_path_vec = if pdb_dir.is_some() {
                    load_path(&pdb_dir.unwrap())
                } else {
                    pdb_path_vec
                };
                let mut controller = Controller::new(pdb_path_vec);
                controller.fill_numeric_id_vec();
                // Iterate per file parallelly
                let hashes = (&controller.path_vec).into_par_iter().map(
                    |pdb_path| {
                        // Read structure from PDB file
                        let pdb_reader = PDBReader::from_file(&pdb_path).expect(
                            &log_msg(FAIL, "Failed to read PDB file")
                        );
                        let compact = pdb_reader.read_structure().expect(
                            &log_msg(FAIL, "Failed to read PDB file")
                        );
                        let compact = compact.to_compact();
                        // let mut hash_vec: Vec<u64> = Vec::with_capacity(compact.num_residues.pow(2));
                        let mut res_pair: Vec<(u16, u16)> = Vec::with_capacity(compact.num_residues.pow(2));
                        // Iterate over all residue pairs
                        let res_bound = get_all_combination(compact.num_residues, false);
                        // Single
                        let mut hash_collection: Vec<u64> = res_bound.into_par_iter().map(
                            |(n, m)| {
                                if n == m {
                                    return 0u64
                                }
                                let trr = compact.get_trrosetta_feature(n, m).unwrap_or([0.0; 6]);
                                if trr[0] < 2.0 || trr[0] > 20.0 {
                                    return 0u64
                                }
                                let hash_value = HashValue::perfect_hash(
                                    compact.residue_serial[n] as u64, compact.residue_serial[m] as u64,
                                    trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]
                                );
                                hash_value.as_u64()
                            }
                        ).collect();
                        
                        // Filter out invalid residue pairs & deduplicate
                        hash_collection.sort_unstable();
                        hash_collection.dedup();
                        
                        // drop unnecessary variables
                        drop(compact);
                        drop(pdb_reader);
                        hash_collection.shrink_to_fit();
                        
                        // Return
                        hash_collection
                    }
                ).collect::<Vec<_>>();
                // Remove invalid hash
                
                
                let lap1 = std::time::Instant::now();
                if verbose { 
                    print_log_msg(INFO, &format!("Time elapsed for getting features {:?}", lap1 - start));
                }

                let mut index_table = DashMap::<u64, Vec<u64>>::new();
                let mut index_table = Arc::new(index_table);
                // for i in (0..hashes.len()).into_par_iter().rev() {
                //     let hash_vec = &hashes[i];
                //     for j in 0..hash_vec.len() {
                //         let hash = hash_vec[j];
                //         if hash == 0 {
                //             continue
                //         }
                //         let numeric_id = controller.numeric_id_vec[i] as u64;
                //         let mut value_vec = index_table.entry(hash).or_insert(Vec::new());
                //         value_vec.push(numeric_id);
                //     }
                // }
                // FILL IN INDEX TABLE WITH PARALLEL
                (&hashes).into_par_iter().enumerate().for_each(
                    |(i, hash_vec)| {
                        for j in 0..hash_vec.len() {
                            let hash = hash_vec[j];
                            if hash == 0 {
                                continue
                            }
                            let numeric_id = controller.numeric_id_vec[i] as u64;
                            let mut value_vec = index_table.entry(hash).or_insert(Vec::new());
                            value_vec.push(numeric_id);
                        }
                    }
                );

                let lap2 = std::time::Instant::now();

                if verbose { 
                    print_log_msg(INFO, &format!("Time elapsed for filling index table {:?}", lap2 - lap1));
                }
                
                // index_table.remove(&0u64); // Remove invalid hash
                // Convert to offset table
                let index_table = Arc::try_unwrap(index_table).unwrap();
                let (offset_table, value_vec) = convert_hashmap_to_offset_and_values(index_table);
                let lap3 = std::time::Instant::now();
                if verbose { 
                    print_log_msg(INFO, &format!("Time elapsed for converting index table {:?}", lap3 - lap2));
                }
                // Save offset table
                let offset_path = format!("{}.offset", index_path);
                save_offset_map(&offset_path, &offset_table).expect(
                    &log_msg(FAIL, "Failed to save offset table")
                );
                let lap4 = std::time::Instant::now();
                if verbose { 
                    print_log_msg(INFO, &format!("Time elapsed for saving offset table {:?}", lap4 - lap3));
                }
                
                // Save value vector
                let value_path = format!("{}.value", index_path);
                write_u64_vector(&value_path, &value_vec).expect(
                    &log_msg(FAIL, "Failed to save values")
                );
                let lap5 = std::time::Instant::now();
                if verbose { 
                    print_log_msg(INFO, &format!("Time elapsed for saving value vector {:?}", lap5 - lap4));
                }
                
                // Save lookup. The path to lookup table is the same as the index table with .lookup extension
                let lookup_path = format!("{}.lookup", index_path);
                save_lookup_to_file(&lookup_path, &controller.path_vec, &controller.numeric_id_vec, None);
                let lap6 = std::time::Instant::now();
                if verbose { 
                    print_log_msg(INFO, &format!("Time elapsed for saving lookup table {:?}", lap6 - lap5));
                }
                // TEMP: load index table saved
                let lap7 = std::time::Instant::now();
                let offset_table = read_offset_map(&offset_path).expect("[ERROR] Failed to load offset table");
                let lap8 = std::time::Instant::now();

                if verbose {
                    print_log_msg(INFO, &format!("Time elapsed for loading offset table {:?}", lap8 - lap7));
                }
                let value_vec = read_u64_vector(&value_path).expect("[ERROR] Failed to load value vector");
                let lap9 = std::time::Instant::now();
                if verbose {
                    print_log_msg(INFO, &format!("Time elapsed for loading value vector {:?}", lap9 - lap8));
                }
                let (path_vec, numeric_id_vec, optional_vec) = load_lookup_from_file(&lookup_path);
                let lap10 = std::time::Instant::now();
                if verbose {
                    print_log_msg(INFO, &format!("Time elapsed for loading lookup table {:?}", lap10 - lap9));
                }
                let lap11 = std::time::Instant::now();
                let offset = offset_table.get(&hashes[1][3]).unwrap();
                let values = get_values_with_offset(&value_vec, offset.0, offset.1);
                if verbose {
                    print_log_msg(INFO, &format!("Time elapsed for getting values {:?}", lap11 - lap10));
                }
                if verbose { 
                    print_log_msg(DONE, "Done.");
                }
            }
        }
        _ => {
            eprintln!("{}", HELP_INDEX);
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_index() {
        let pdb_dir = "data/serine_peptidases_filtered";
        let pdb_path_vec = load_path(pdb_dir);
        let index_path = "data/serine_peptidases_index";
        let num_threads = 6;
        let verbose = true;
        let help = false;
        let env = AppArgs::Index {
            pdb_dir: Some(pdb_dir.to_string()),
            pdb_path_vec,
            index_path: index_path.to_string(),
            num_threads,
            verbose,
            help,
        };
        build_index(env);
    }
}