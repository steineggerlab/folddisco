// File: build_index.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for building index table
// For building index table, we need a directory containing PDB files, 
// and the path to save the index table.


use crate::cli::*;
use crate::controller::io::save_offset_vec;
use crate::index::alloc::convert_sorted_pairs_to_offset_and_values_vec;
use crate::prelude::*;
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

pub const HELP_INDEX: &str = "\
USAGE: motifsearch index [OPTIONS]
Options:
    -d, --pdbs <PDB_DIR>         Directory containing PDB files
    -H, --hash <HASH_TYPE>     Hash type to use (pdb, trrosetta, default)
    -i, --index <INDEX_PATH>   Path to save the index table
    -t, --threads <THREADS>         Number of threads to use
    -v, --verbose                   Print verbose messages
    -h, --help                      Print this help menu
";
// TODO: ADD MONITOR_MEMORY AS A PARAMETER
// TODO: ADD NBIN_ANGLE, NBIN_DIST AS PARAMETERS

pub fn build_index(env: AppArgs) {
    match env {
        AppArgs::Index {
            pdb_dir,
            pdb_path_vec,
            hash_type,
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
                let pdb_dir_clone = pdb_dir.clone();
                if verbose { 
                    print_log_msg(
                        INFO,&format!(
                            "Indexing {} with {} threads", 
                            pdb_dir_clone.unwrap_or("None".to_string()),
                            num_threads
                        )
                    );
                }
                // Setup multithreading
                rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global().unwrap();
                
                // Load PDB files
                let pdb_path_vec = if pdb_dir.is_some() {
                    load_path(&pdb_dir.unwrap())
                } else {
                    pdb_path_vec
                };
                let hash_type = HashType::get_with_str(hash_type.as_str());

                let mut fold_disco = FoldDisco::new_with_params(
                    pdb_path_vec, hash_type, true, num_threads, index_path.clone()
                );
                
                // Main workflow
                // 2. Collect hash
                measure_time!(fold_disco.collect_hash_pairs());
                if verbose { print_log_msg(INFO, &format!("Hash collected (Allocated {}MB)", PEAK_ALLOC.current_usage_as_mb())); }
                measure_time!(fold_disco.sort_hash_pairs());
                if verbose { print_log_msg(INFO, &format!("Hash sorted (Allocated {}MB)", PEAK_ALLOC.current_usage_as_mb())); }
                measure_time!(fold_disco.fill_numeric_id_vec());

                let (offset_table, value_vec) =
                    // measure_time!(fold_disco.index_builder.convert_sorted_pairs_to_offset_and_values(fold_disco.hash_id_pairs));
                    measure_time!(convert_sorted_pairs_to_offset_and_values_vec(fold_disco.hash_id_pairs));
                if verbose { print_log_msg(INFO, &format!("Converted to offsets (Allocated {}MB)", PEAK_ALLOC.current_usage_as_mb())); }
                // Save offset table
                let offset_path = format!("{}.offset", index_path);
                // measure_time!(save_offset_map(&offset_path, &offset_table).expect(
                measure_time!(save_offset_vec(&offset_path, &offset_table).expect(
                    &log_msg(FAIL, "Failed to save offset table")
                ));
                drop(offset_table);

                // Save value vector  
                let value_path = format!("{}.value", index_path);
                measure_time!(write_usize_vector(&value_path, &value_vec).expect(
                    &log_msg(FAIL, "Failed to save values")
                ));
                drop(value_vec);

                // Save lookup. The path to lookup table is the same as the index table with .lookup extension
                let lookup_path = format!("{}.lookup", index_path);
                measure_time!(save_lookup_to_file(
                    &lookup_path, &fold_disco.path_vec,
                     &fold_disco.numeric_id_vec, None
                ));

                // Save hash type
                let hash_type_path = format!("{}.type", index_path);
                hash_type.save_to_file(&hash_type_path);
                
                // eprintln!("FINAL - {}MB", PEAK_ALLOC.current_usage_as_mb());
                if verbose { print_log_msg(DONE, "Done."); }
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
        let hash_type = "default32";
        let index_path = "data/serine_peptidases_default32";
        let num_threads = 4;
        let verbose = true;
        let help = false;
        let env = AppArgs::Index {
            pdb_dir: Some(pdb_dir.to_string()),
            pdb_path_vec,
            hash_type: hash_type.to_string(),
            index_path: index_path.to_string(),
            num_threads,
            verbose,
            help,
        };
        build_index(env);
    }
}