// File: build_index.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

//! This file contains the workflow for building index table
//! For building index table, we need a directory containing PDB files, 
//! and the path to save the index table.


use std::path::PathBuf;

use crate::cli::config::{write_index_config_to_file, IndexConfig};
use crate::controller::map::{convert_sorted_hash_pairs_to_simplemap, convert_sorted_hash_vec_to_simplemap};
use crate::controller::mode::{parse_path_vec_by_id_type, IdType, IndexMode};
use crate::cli::*;
use crate::controller::io::write_usize_vector_in_bits;
use crate::prelude::*;
use memmap2::MmapMut;
use peak_alloc::PeakAlloc;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

pub const HELP_INDEX: &str = "\
USAGE: folddisco index [OPTIONS]
Options:
    -p, --pdbs <PDB_DIR>         Directory containing PDB files
    -y, --type <HASH_TYPE>       Hash type to use (pdb, trrosetta, default)
    -i, --index <INDEX_PATH>     Path to save the index table
    -m, --mode <MODE>            Mode to index (default=id, id: index only id, grid: index id and grid, pos: index id and position)
    -t, --threads <THREADS>      Number of threads to use
    -d, --distance <NBIN_DIST>   Number of distance bins (default 0, zero means default)
    -a, --angle <NBIN_ANGLE>     Number of angle bins (default 0, zero means default)
    -g, --grid <GRID_WIDTH>      Grid width (default 30.0)
    -c, --chunk <CHUNK_SIZE>     Number of PDB files to index at once (default, max=65535)
    -r, --recursive              Index PDB files in subdirectories
    -n, --max-residue <MAX_RES>  Maximum number of residues in a PDB file (default=3000)
    --id <ID_TYPE>               ID type to use (pdb, uniprot, afdb, relpath, abspath, default=relpath)
    -v, --verbose                Print verbose messages
    -h, --help                   Print this help menu
";

pub fn build_index(env: AppArgs) {
    match env {
        AppArgs::Index {
            pdb_dir,
            hash_type,
            index_path,
            num_threads,
            mode,
            max_residue,
            num_bin_dist,
            num_bin_angle,
            grid_width,
            chunk_size,
            recursive,
            id_type,
            verbose,
            help,
        } => {
            // Check if arguments are valid
            if pdb_dir.is_none() {
                eprintln!("{}", HELP_INDEX);
                std::process::exit(1);
            }
            if help {
                eprintln!("{}", HELP_INDEX);
                std::process::exit(0);
            }
            let pdb_dir_clone = pdb_dir.clone();
            // Load PDB files
            let pdb_path_vec = if pdb_dir.is_some() {
                load_path(&pdb_dir.unwrap(), recursive)
            } else {
                print_log_msg(FAIL, "Directory containing PDB files is not provided");
                eprintln!("{}", HELP_INDEX);
                std::process::exit(1);
            };
            let chunk_size = if chunk_size > u16::max_value() as usize { u16::max_value() as usize } else { chunk_size };
            let num_chunks = if pdb_path_vec.len() < chunk_size { 1 } else { (pdb_path_vec.len() as f64 / chunk_size as f64).ceil() as usize };
            if verbose { 
                print_log_msg(
                    INFO,&format!(
                        "Indexing {} with {} threads and {} chunks",
                        pdb_dir_clone.unwrap_or("None".to_string()),
                        num_threads, num_chunks
                    )
                );
            }
            
            let hash_type = HashType::get_with_str(hash_type.as_str());
            if verbose { print_log_msg(INFO, &format!("Hash type: {:?}", hash_type)); }
            let pdb_path_chunks = pdb_path_vec.chunks(chunk_size);
            let index_mode = IndexMode::get_with_str(mode.as_str());
            let id_type = IdType::get_with_str(id_type.as_str());
            
            pdb_path_chunks.into_iter().enumerate().for_each(|(i, pdb_path_vec)| {
                let index_path = if num_chunks == 1 {
                    if verbose { print_log_msg(INFO, "Indexing all PDB files in one chunk"); }
                    index_path.clone()
                } else {
                    if verbose { print_log_msg(INFO, &format!("Indexing chunk {}", i)); }
                    format!("{}_{}", index_path, i)
                };
                print_log_msg(INFO, &format!("Before initializing (Allocated {}MB)", PEAK_ALLOC.current_usage_as_mb()));
                let mut fold_disco = FoldDisco::new_with_params(
                    pdb_path_vec.to_vec(), hash_type, true, num_threads, 
                    num_bin_dist, num_bin_angle, index_path.clone(), grid_width,
                );
                
                match index_mode {
                    IndexMode::Id => {
                        if verbose { print_log_msg(INFO, "Collecting ids of the structures"); }
                        measure_time!(fold_disco.collect_hash_vec());
                        if verbose {
                            print_log_msg(INFO, 
                                &format!("Total {} hashes collected (Allocated {}MB)", fold_disco.hash_id_vec.len(), PEAK_ALLOC.current_usage_as_mb())
                            );
                        }
                        measure_time!(fold_disco.sort_hash_vec());
                    }
                    IndexMode::Grid => {
                        if verbose { print_log_msg(INFO, "Collecting ids and grids with hashes observed"); }
                        measure_time!(fold_disco.collect_hash_grids());
                        if verbose {
                            print_log_msg(INFO, 
                                &format!("Total {} hashes collected (Allocated {}MB)", fold_disco.hash_id_grids.len(), PEAK_ALLOC.current_usage_as_mb())
                                // &format!("Total {} hashes collected", fold_disco.hash_id_grids.len())
                            );
                        }
                        measure_time!(fold_disco.sort_hash_grids());
                    }
                    IndexMode::Pos => {
                        if verbose { print_log_msg(INFO, "Collecting ids and positions with hashes"); }
                        //measure_time!(fold_disco.collect_hash_positions());
                        // if verbose {
                        //     print_log_msg(INFO, 
                        //         &format!("Total {} hashes collected (Allocated {}MB)", fold_disco.hash_id_pos.len(), PEAK_ALLOC.current_usage_as_mb())
                        //     );
                        // }
                        todo!("Implement this part");
                    }
                    IndexMode::Big => {
                        if verbose { print_log_msg(INFO, "Collecting ids of the structures"); }
                        measure_time!(fold_disco.collect_and_count());
                        if verbose {
                            print_log_msg(INFO, 
                                &format!("Hashes collected (Allocated {}MB)", PEAK_ALLOC.current_usage_as_mb())
                            );
                        }
                        measure_time!(fold_disco.fold_disco_index.allocate_entries());
                        measure_time!(fold_disco.add_entries());
                        measure_time!(fold_disco.fold_disco_index.finish_index());
                        measure_time!(fold_disco.fold_disco_index.save_offset_to_file());
                    }
                }
                if verbose { print_log_msg(INFO,
                    &format!("Hash sorted (Allocated {}MB)", PEAK_ALLOC.current_usage_as_mb())
                    // "Hash sorted"
                ); }
                fold_disco.fill_numeric_id_vec();

                let (offset_map, value_vec) = match index_mode {
                    IndexMode::Id => {
                        // measure_time!(convert_sorted_hash_pairs_to_simplemap(fold_disco.hash_id_pairs))
                        measure_time!(convert_sorted_hash_vec_to_simplemap(fold_disco.hash_id_vec))
                    }
                    IndexMode::Grid => {
                        measure_time!(convert_sorted_hash_pairs_to_simplemap(fold_disco.hash_id_grids))
                    }
                    IndexMode::Pos => {
                        todo!("Implement this part");
                    }
                    _ => {
                        unimplemented!();
                    }
                };

                if verbose { print_log_msg(INFO, 
                    &format!("Offset & values acquired (Allocated {}MB)", PEAK_ALLOC.current_usage_as_mb())
                    // &format!("Offset & values acquired")
                ); }
                let offset_path = format!("{}.offset", index_path);
                // measure_time!(save_offset_vec(&offset_path, &offset_table).expect(
                //     &log_msg(FAIL, "Failed to save offset table")
                // ));
                measure_time!(offset_map.dump_to_disk(&PathBuf::from(&offset_path)).expect(
                    &log_msg(FAIL, "Failed to save offset table")
                ));
                let value_path = format!("{}.value", index_path);
                match index_mode {
                    IndexMode::Id => {
                        measure_time!(write_usize_vector_in_bits(&value_path, &value_vec, 16).expect(
                            &log_msg(FAIL, "Failed to save values")
                        ));
                    }
                    IndexMode::Grid => {
                        measure_time!(write_usize_vector_in_bits(&value_path, &value_vec, 24).expect(
                            &log_msg(FAIL, "Failed to save values")
                        ));
                    }
                    IndexMode::Pos => {
                        todo!("Implement this part");
                        // measure_time!(write_usize_vector_in_bits(&value_path, &value_vec, 48).expect(
                        //     &log_msg(FAIL, "Failed to save values")
                        // ));
                    }
                    _ => {
                        unimplemented!();
                    }
                }
                let lookup_path = format!("{}.lookup", index_path);
                let id_vec = parse_path_vec_by_id_type(&fold_disco.path_vec, id_type.clone());
                measure_time!(save_lookup_to_file(
                    &lookup_path, &id_vec, &fold_disco.numeric_id_vec,
                    Some(&fold_disco.nres_vec), Some(&fold_disco.plddt_vec)
                ));
                let hash_type_path = format!("{}.type", index_path);
                let chunk_size = pdb_path_vec.len();
                let index_config = IndexConfig::new(
                    hash_type, num_bin_dist, num_bin_angle, index_mode.clone(),
                    grid_width, chunk_size, max_residue
                );
                write_index_config_to_file(&hash_type_path, index_config);
                if verbose { print_log_msg(DONE, &format!("Indexing done for chunk {} - {}", i, index_path)); }
            });
            if verbose { print_log_msg(DONE, "Done."); }
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
        // let pdb_dir = "analysis/e_coli/sample";
        let hash_type = "pdbtr";
        // let index_path = "analysis/e_coli/test_index";
        let index_path = "data/serine_peptidases_pdbtr_test";
        let index_mode = "big";
        let num_threads = 8;
        let num_bin_dist = 16;
        let num_bin_angle = 4;
        let chunk_size = 65535;
        let max_residue = 3000;
        let grid_width = 40.0;
        let recursive = true;
        let verbose = true;
        let help = false;
        let id_type = "relpath";
        let env = AppArgs::Index {
            pdb_dir: Some(pdb_dir.to_string()),
            hash_type: hash_type.to_string(),
            index_path: index_path.to_string(),
            mode: index_mode.to_string(),
            num_threads,
            num_bin_dist,
            num_bin_angle,
            grid_width,
            chunk_size,
            max_residue,
            recursive,
            id_type: id_type.to_string(),
            verbose,
            help,
        };
        build_index(env);
    }
}