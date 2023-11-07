// File: build_index.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for building index table
// For building index table, we need a directory containing PDB files, 
// and the path to save the index table.

use crate::cli::*;
use rayon::prelude::*;
use std::collections::HashMap;
// Import unique
use std::collections::HashSet;
use std::hash;

use crate::prelude::*;

// TODO: Generated. need to be

const HELP_INDEX: &str = "\
USAGE: motifsearch index [OPTIONS] <PDBS...>
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
            if help {
                println!("{}", HELP_INDEX);
            } else {
                if verbose { println!("[INFO ] Building index table..."); }
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
                let hashes = (controller.path_vec).into_par_iter().map(
                    |pdb_path| {
                        // Read structure from PDB file
                        let pdb_reader = PDBReader::from_file(pdb_path).expect("[ERROR] PDB file not found");
                        let compact = pdb_reader.read_structure().expect("[ERROR] Failed to read PDB file");
                        let compact = compact.to_compact();
                        let mut hash_vec: Vec<u64> = Vec::with_capacity(compact.num_residues.pow(2));
                        let mut res_pair: Vec<(u16, u16)> = Vec::with_capacity(compact.num_residues.pow(2));
                        let res_bound = get_all_combination(compact.num_residues, false);
                        let hash_collection: Vec<u64> = res_bound.iter().map(
                            |(n, m)| {
                                if n == m {
                                    return 0u64
                                }
                                let trr = compact.get_trrosetta_feature(*n, *m).unwrap_or([0.0; 6]);
                                if trr[0] < 2.0 || trr[0] > 20.0 {
                                    return 0u64
                                }
                                let hash_value = HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);
                                hash_value.as_u64()
                            }
                        ).collect();
                        hash_collection
                    }
                ).collect::<Vec<_>>();
                let lap1 = std::time::Instant::now();
                if verbose { println!("[INFO ] Time elapsed for building index table {:?}", lap1 - start); }
                let mut index_table = IndexTable::new();
                index_table.fill(&controller.numeric_id_vec, &hashes);
                index_table.save_to_bin(&index_path).expect("[ERROR] Failed to save index table");
                let lap2 = std::time::Instant::now();
                if verbose { println!("[INFO ] Time elapsed for saving index table {:?}", lap2 - lap1); }
                if verbose { println!("[INFO ] Done!"); }
            }
        }
        _ => {
            println!("Invalid subcommand. Try `motifsearch --help` for more information.");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_index() {
        let pdb_dir = "analysis/raw_ecoli";
        let pdb_path_vec = load_path(pdb_dir);
        let index_path = "data/index.bin";
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