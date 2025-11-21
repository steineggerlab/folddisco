// File: build_index.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

//! This file contains the workflow for building index table
//! For building index table, we need a directory containing PDB files, 
//! and the path to save the index table.


use std::path::PathBuf;

use crate::cli::config::{write_index_config_to_file, IndexConfig};
use crate::controller::map::convert_sorted_hash_vec_to_simplemap;
use crate::controller::mode::{parse_path_vec_by_id_type, IdType, IndexMode};
use crate::cli::*;
use crate::controller::io::{default_index_path, write_usize_vector_in_bits};
use crate::prelude::*;
use crate::structure::io::StructureFileFormat;

#[cfg(feature= "foldcomp")]
use crate::structure::io::fcz::*;

#[cfg(feature = "foldcomp")]
use rayon::prelude::ParallelSliceMut;

pub const HELP_INDEX: &str = "\
usage: folddisco index -p <i:PDB_DIR>|<i:FOLDCOMP_DB> -i <o:INDEX_PATH> [OPTIONS]

input/output:
 -p, --pdbs <PATH>                Directory or Foldcomp DB containing PDB files
 -i, --index <PATH>               Path to save the index table
 -r, --recursive                  Index PDB files in subdirectories recursively

indexing parameters:
 -t, --threads <INT>              Number of threads to use [1]
 -n, --max-residue <INT>          Maximum number of residues in a PDB file [50000]
 --id <STR>                       ID type to use (pdb, uniprot, afdb, relpath, abspath) [relpath]
 -m, --mode <MODE>                Mode to index [id]
    id: suitable for smaller dataset (N < 65536) with hashmap offset;
    big: 8GB fixed-size offset table, suitable for large dataset

hashing parameters:
 -y, --type STR                   Hash type to use (default, pdb, trrosetta, ppf, 3di) [default]
 -d, --distance INT               Number of distance bins [default, 16]
 -a, --angle INT                  Number of angle bins [default, 4]
 --multiple-bins STR              Multiple bins for distance and angle (dist1-ang1,dist2-ang2 e.g. 16-4,8-3)
                                  While increasing sensitivity, this option increases the size of the index.

general options:
 -v, --verbose                    Print verbose messages
 -h, --help                       Print this help menu
 
examples:
# Default. h_sapiens directory or foldcomp database is indexed with default parameters
folddisco index -p h_sapiens -i index/h_sapiens -t 12

# Indexing big protein dataset
folddisco index -p swissprot -i index/swissprot -t 64 -m big -v

# Indexing with custom hash type and parameters
folddisco index -p h_sapiens -i index/h_sapiens -t 12 -y default -d 16 -a 4 # Default
folddisco index -p h_sapiens -i index/h_sapiens -t 12 -y pdb -d 8 -a 3 # PDB
";

pub fn build_index(env: AppArgs) {
    match env {
        AppArgs::Index {
            pdb_container,
            hash_type,
            index_path,
            num_threads,
            mode,
            max_residue,
            num_bin_dist,
            num_bin_angle,
            multiple_bins,
            grid_width,
            chunk_size,
            recursive,
            mmap_on_disk,
            id_type,
            verbose,
            help: _,
        } => {
            if verbose { print_logo(); }
            // Check if arguments are valid
            if pdb_container.is_none() {
                eprintln!("{}", HELP_INDEX);
                std::process::exit(1);
            }
            // help is handled in the main function
            let pdb_container_clone = pdb_container.clone();
            let index_path = if index_path.is_empty() {
                // Get default index path
                default_index_path(pdb_container.as_ref().unwrap())
            } else { 
                index_path
            };
            #[cfg(feature = "foldcomp")]
            let pdb_container_name: &'static str = Box::leak(pdb_container.clone().unwrap().into_boxed_str());
            #[allow(unused_mut)]
            let mut input_format: StructureFileFormat = StructureFileFormat::PDB;
            // Load PDB files
            let pdb_path_vec = if pdb_container.is_some() {
                // Check if pdb_dir is a directory or db file
                // If not foldcomp, just load
                #[cfg(not(feature = "foldcomp"))]
                { load_path(&pdb_container.unwrap(), recursive) }
                #[cfg(feature = "foldcomp")]
                {
                    let pdb_container = pdb_container.unwrap();
                    // Check if pdb_container is a directory or a file
                    let is_dir = PathBuf::from(&pdb_container).is_dir();
                    if is_dir {
                        load_path(&pdb_container, recursive)
                    } else {
                        let mut lookup_vec = read_foldcomp_db_lookup(&pdb_container).expect(
                            &log_msg(FAIL, "Failed to read Foldcomp DB lookup")
                        );
                        lookup_vec.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
                        let mut index_vec = read_foldcomp_db_index(&pdb_container).expect(
                            &log_msg(FAIL, "Failed to read Foldcomp DB index")
                        );
                        index_vec.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
                        input_format = StructureFileFormat::FCZDB;
                        get_path_vector_out_of_lookup_and_index(&lookup_vec, &index_vec)
                    }
                }
            } else {
                print_log_msg(FAIL, "Directory containing PDB files is not provided");
                eprintln!("{}", HELP_INDEX);
                std::process::exit(1);
            };
            
            let index_mode = IndexMode::get_with_str(mode.as_str());
            // Set index_mode if pdb_path_vec.len() > 65536
            let index_mode = if pdb_path_vec.len() > u16::max_value() as usize {
                if verbose {
                    print_log_msg(
                        INFO, &format!(
                            "Number of PDB files ({}) exceeds 65536, switching to Big mode",
                            pdb_path_vec.len()
                        )
                    );
                }
                IndexMode::Big
            } else {
                index_mode
            };
            if index_mode == IndexMode::Big && verbose {
                print_log_msg(INFO, "Indexing in Big mode.");
            }
            let chunk_size = if chunk_size > u16::max_value() as usize { u16::max_value() as usize } else { chunk_size };
            // Overwrite chunk_size if index_mode is Big
            let chunk_size = if index_mode == IndexMode::Big { pdb_path_vec.len() } else { chunk_size };
            let num_chunks = if pdb_path_vec.len() <= chunk_size { 1 } else { (pdb_path_vec.len() as f64 / chunk_size as f64).ceil() as usize };
            if verbose { 
                print_log_msg(
                    INFO,&format!(
                        "Indexing {} with {} threads and {} chunks",
                        pdb_container_clone.unwrap_or("None".to_string()),
                        num_threads, num_chunks
                    )
                );
            }
            let hash_type = HashType::get_with_str(hash_type.as_str());
            if verbose { print_log_msg(INFO, &format!("Hash type: {:?}", hash_type)); }
            let pdb_path_chunks = pdb_path_vec.chunks(chunk_size);
            let id_type = IdType::get_with_str(id_type.as_str());
            
            let multiple_bins = if let Some(multiple_bins) = multiple_bins {
                Some(parse_pairs(&multiple_bins))
            } else {
                None
            };
            
            pdb_path_chunks.into_iter().enumerate().for_each(|(i, pdb_path_vec)| {
                // let pdb_container_name_inner: &'static str = pdb_container_name.clone();
                let index_path = if num_chunks == 1 {
                    if verbose { print_log_msg(INFO, "Indexing all PDB files in one chunk"); }
                    index_path.clone()
                } else {
                    if verbose { print_log_msg(INFO, &format!("Indexing chunk {}", i)); }
                    format!("{}_{}", index_path, i)
                };

                #[cfg(not(feature = "foldcomp"))]
                let mut folddisco = FoldDisco::new(
                    pdb_path_vec.to_vec(), hash_type, num_threads, 
                    num_bin_dist, num_bin_angle, index_path.clone(), 
                    grid_width, index_mode, multiple_bins.clone(),
                    mmap_on_disk,
                );
                #[cfg(feature = "foldcomp")]
                let mut folddisco = if PathBuf::from(&pdb_container_name).is_dir() {
                    FoldDisco::new(
                        pdb_path_vec.to_vec(), hash_type, num_threads, 
                        num_bin_dist, num_bin_angle, index_path.clone(), 
                        grid_width, index_mode, multiple_bins.clone(),
                        mmap_on_disk,
                    )
                } else {
                    FoldDisco::new_with_foldcomp_db(
                        pdb_path_vec.to_vec(), hash_type, num_threads, 
                        num_bin_dist, num_bin_angle, index_path.clone(), 
                        grid_width, index_mode, pdb_container_name,
                        multiple_bins.clone(), mmap_on_disk,
                    )
                };
                
                match index_mode {
                    IndexMode::Id => {
                        if verbose {
                            print_log_msg(INFO, "Collecting ids of the structures"); 
                            measure_time!(folddisco.collect_hash_vec());
                            print_log_msg(INFO, 
                                &format!("Total {} hashes collected", folddisco.hash_id_vec.len())
                            );
                            measure_time!(folddisco.sort_hash_vec());
                        } else {
                            folddisco.collect_hash_vec();
                            folddisco.sort_hash_vec();
                        }
                    }
                    IndexMode::Big => {
                        if verbose {
                            print_log_msg(INFO, "Collecting ids of the structures");
                            measure_time!(folddisco.collect_and_count());
                            print_log_msg(INFO, &format!("Hashes collected"));
                            measure_time!(folddisco.fold_disco_index.allocate_entries());
                            measure_time!(folddisco.add_entries());
                            measure_time!(folddisco.fold_disco_index.finish_index());
                            // TESTING here
                            let pruned_index = measure_time!(folddisco.fold_disco_index.prune_to_sparse());
                            // measure_time!(folddisco.fold_disco_index.save_offset_to_file());
                            measure_time!(pruned_index.save_offset_to_file());
                            
                        } else {
                            folddisco.collect_and_count();
                            folddisco.fold_disco_index.allocate_entries();
                            folddisco.add_entries();
                            folddisco.fold_disco_index.finish_index();
                            // TESTING here
                            let pruned_index = folddisco.fold_disco_index.prune_to_sparse();
                            // folddisco.fold_disco_index.save_offset_to_file();
                            pruned_index.save_offset_to_file();
                        }
                    }
                }
                if verbose { print_log_msg(INFO,"Hash sorted"); }
                folddisco.fill_numeric_id_vec();

                match index_mode {
                    IndexMode::Id => {
                        let offset_path = format!("{}.offset", index_path);
                        let value_path = index_path.clone(); // Changed to new format without .value extension 2025-08-08 21:33:08
                        if verbose {
                            let (offset_map, value_vec) = measure_time!(convert_sorted_hash_vec_to_simplemap(folddisco.hash_id_vec));
                            print_log_msg(INFO, &format!("Offset & values acquired"));
                            measure_time!(offset_map.dump_to_disk(&PathBuf::from(&offset_path)).expect(
                               &log_msg(FAIL, "Failed to save offset table")
                            ));
                            measure_time!(write_usize_vector_in_bits(&value_path, &value_vec, 16).expect(
                                &log_msg(FAIL, "Failed to save values")
                            ));
                        } else {
                            let (offset_map, value_vec) = convert_sorted_hash_vec_to_simplemap(folddisco.hash_id_vec);
                            offset_map.dump_to_disk(&PathBuf::from(&offset_path)).expect(
                                &log_msg(FAIL, "Failed to save offset table")
                            );
                            write_usize_vector_in_bits(&value_path, &value_vec, 16).expect(
                                &log_msg(FAIL, "Failed to save values")
                            );
                        }
                    }
                    IndexMode::Big => {}
                }
                let lookup_path = format!("{}.lookup", index_path);
                let id_vec = parse_path_vec_by_id_type(&folddisco.path_vec, &id_type);
                // Check if db_key_vec is empty or not
                let key_vec_input = if folddisco.numeric_db_key_vec.is_empty() {
                    None
                } else {
                    Some(&folddisco.numeric_db_key_vec)
                };
                if verbose {
                    measure_time!(save_lookup_to_file(
                        &lookup_path, &id_vec, &folddisco.numeric_id_vec,
                        Some(&folddisco.nres_vec), Some(&folddisco.plddt_vec), key_vec_input
                    ));
                } else {
                    save_lookup_to_file(
                        &lookup_path, &id_vec, &folddisco.numeric_id_vec,
                        Some(&folddisco.nres_vec), Some(&folddisco.plddt_vec),
                        key_vec_input
                    );
                }

                let hash_type_path = format!("{}.type", index_path);
                let chunk_size = pdb_path_vec.len();
                #[cfg(not(feature = "foldcomp"))]
                let index_config = IndexConfig::new(
                    hash_type, num_bin_dist, num_bin_angle, index_mode.clone(),
                    grid_width, chunk_size, max_residue, input_format.clone(),
                    None, multiple_bins.clone(),
                );
                #[cfg(feature = "foldcomp")]
                let index_config = IndexConfig::new(
                    hash_type, num_bin_dist, num_bin_angle, index_mode.clone(),
                    grid_width, chunk_size, max_residue, input_format.clone(), 
                    Some(pdb_container_name.to_string()), multiple_bins.clone(),
                );
                write_index_config_to_file(&hash_type_path, index_config);
                if verbose { print_log_msg(DONE, &format!("Indexing done for chunk {} - {}", i+1, index_path)); }
            });
            if verbose { print_log_msg(DONE, "Done."); }
        }
        _ => {
            eprintln!("{}", HELP_INDEX);
            std::process::exit(1);
        }
    }
}

fn parse_pairs(input: &str) -> Vec<(usize, usize)> {
    input
        .split(',')
        .filter_map(|pair_str| {
            let parts: Vec<_> = pair_str.trim().split('-').collect();
            if parts.len() == 2 {
                match (parts[0].parse(), parts[1].parse()) {
                    (Ok(a), Ok(b)) => Some((a, b)),
                    _ => None,
                }
            } else {
                None
            }
        })
        .collect()
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_build_index() {
        let env = AppArgs::Index {
            pdb_container: Some("data/serine_peptidases".to_string()),
            hash_type: "pdbtr".to_string(),
            index_path: "data/serine_peptidases_pdbtr_small".to_string(),
            mode: "id".to_string(),
            num_threads: 8,
            num_bin_dist: 16,
            num_bin_angle: 4,
            multiple_bins: None,
            grid_width: 20.0,
            chunk_size: 65536,
            max_residue: 3000,
            recursive: true,
            mmap_on_disk: false,
            id_type: "relpath".to_string(),
            verbose: true,
            help: false,
        };
        build_index(env);
    }
    #[test]
    fn test_build_index_of_foldcomp_db() {
        #[cfg(feature = "foldcomp")]
        {
            let env = AppArgs::Index {
                pdb_container: Some("data/foldcomp/example_db".to_string()),
                hash_type: "pdbtr".to_string(),
                index_path: "data/example_db_folddisco_db".to_string(),
                mode: "id".to_string(),
                num_threads: 8,
                num_bin_dist: 16,
                num_bin_angle: 4,
                multiple_bins: Some("16-4,8-3".to_string()),
                grid_width: 20.0,
                chunk_size: 65536,
                max_residue: 50000,
                recursive: true,
                mmap_on_disk: false,
                id_type: "relpath".to_string(),
                verbose: true,
                help: false,
            };
            build_index(env);
        }
    }
}