// File: query_pdb.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for querying PDB files
// When querying PDB files, we need index table and query file.

use core::{hash, num};
use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use crate::cli::workflows::config::read_config_from_file;
use crate::cli::*;
use crate::controller::io::{get_values_with_offset, get_values_with_offset_u16, read_offset_map, read_u16_vector, read_usize_vector};
use crate::controller::retrieve::{connected, hash_vec_to_aa_pairs, res_vec_as_string, retrieve_residue_with_hash};
use crate::index::lookup::{load_lookup_from_file};
use crate::prelude::*;


pub const HELP_QUERY: &str = "\
USAGE: motifsearch query [OPTIONS] <QUERY_PDB> <CHAIN1><RES1>,<CHAIN2><RES2>,<CHAIN3><RES3>...
Example: motifsearch query -i index_table.index -t 6 1aq2.pdb A250,A232,A269
Options:
    -d, --pdb <PDB_PATH>         Path of PDB file to query
    -q, --query <QUERY_STRING>   Query string
    -i, --index <INDEX_PATH>     Path of index table to load
    -t, --threads <THREADS>      Number of threads to use
    -v, --verbose                Print verbose messages
    -c, --check-nearby           Check nearby hashes (In development)
    -R, --retrieve               Retrieve matched residues (Need PDB files)
    -m, --match <MATCH_CUTOFF>   Match cutoff (default 0.0)
    -s, --score <SCORE_CUTOFF>   Score cutoff (default 0.0)
    -n, --num-res <NUM_RES_CUTOFF> Number of residues cutoff (default 3000)
    -l, --plddt <PLDDT_CUTOFF>   PLDDT cutoff (default 0.0)
    -h, --help                   Print this help menu
";


pub fn query_pdb(env: AppArgs) {
    match env {
        AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            check_nearby,
            retrieve,
            match_cutoff,
            score_cutoff,
            num_res_cutoff,
            plddt_cutoff,
            help,
        } => {
            if help {
                eprintln!("{}", HELP_QUERY);
                std::process::exit(0);
            }
            // Check if arguments are valid
            if index_path.is_none() {
                eprintln!("{}", HELP_QUERY);
                std::process::exit(1);
            }

            // Get path. formatting without quotation marks
            let index_path = index_path.unwrap();
            // Check if index_path_0 is a file.
            let index_chunk_prefix = format!("{}_0", index_path.clone());
            let index_chunk_path = format!("{}_0.offset", index_path.clone());
            let mut index_paths = Vec::new();
            if std::path::Path::new(&index_chunk_path).is_file() {
                print_log_msg(INFO, &format!("Index table is chunked"));
                let mut i = 0;
                loop {
                    let index_chunk_prefix = format!("{}_{}", index_path.clone(), i);
                    let index_chunk_path = format!("{}.offset", index_chunk_prefix);
                    if std::path::Path::new(&index_chunk_path).is_file() {
                        index_paths.push(index_chunk_prefix);
                        i += 1;
                    } else {
                        break;
                    }
                }
            } else {
                index_paths.push(index_path.clone());
            }
            // Limit threads to the number of index paths
            let threads = std::cmp::min(threads, index_paths.len());
            print_log_msg(INFO, &format!("Found {} index file(s). Querying with {} threads", index_paths.len(), threads));
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            pool.install(|| {
                index_paths.into_par_iter().for_each(|index_path| {
                    let offset_path = format!("{}.offset", index_path.clone()); 
                    let value_path = format!("{}.value", index_path.clone());
                    let lookup_path = format!("{}.lookup", index_path.clone());
                    let hash_type_path = format!("{}.type", index_path.clone());
                    assert!(std::path::Path::new(&offset_path).is_file());
                    assert!(std::path::Path::new(&value_path).is_file());
                    assert!(std::path::Path::new(&lookup_path).is_file());
                    assert!(std::path::Path::new(&hash_type_path).is_file());
                    let (hash_type, num_bin_dist, num_bin_angle) = read_config_from_file(&hash_type_path);
                    
                    // Make query with pdb
                    let query_residues = parse_query_string(&query_string);
                    let pdb_query =  make_query(
                        &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, check_nearby
                    );
                    
                    // Load index table
                    // let (mmap, value_vec) = measure_time!(read_usize_vector(&value_path).expect("[ERROR] Failed to load value vector"));
                    let (mmap, value_vec) = measure_time!(read_u16_vector(&value_path).expect("[ERROR] Failed to load value vector"));
                    
                    let offset_table = measure_time!(
                        read_offset_map(&offset_path, hash_type).expect("[ERROR] Failed to load offset table")
                    );
                    // let (path_vec, numeric_id_vec, optional_vec) = measure_time!(load_lookup_from_file(&lookup_path));
                    let lookup = measure_time!(load_lookup_from_file(&lookup_path));

                    // Get offset from offset_table with query
                    let mut offset_to_query = Vec::new();
                    for i in 0..pdb_query.len() {
                        let offset = offset_table.get(&pdb_query[i]);
                        // Get offset map values
                        if offset.is_none() {
                            continue;
                        } else {
                            offset_to_query.push(*offset.unwrap());
                        }
                    }
                    
                    // Make union of values queried
                    let mut query_count_map: HashMap<usize, (usize, f32)> = HashMap::new();
                    for i in 0..offset_to_query.len() {
                        // let single_queried_values = get_values_with_offset(&value_vec, offset_to_query[i].0, offset_to_query[i].1);
                        let single_queried_values = get_values_with_offset_u16(&value_vec, offset_to_query[i].0, offset_to_query[i].1);
                        let hash_count = offset_to_query[i].1;
                        for j in 0..single_queried_values.len() {
                            let nid = lookup.1[single_queried_values[j] as usize];
                            let count = query_count_map.get(&nid);
                            // Added calculation of idf
                            if count.is_none() {
                                query_count_map.insert(nid, (1, (lookup.0.len() as f32 / hash_count as f32).log2()));
                            } else {
                                query_count_map.insert(nid, (count.unwrap().0 + 1, (lookup.0.len() as f32 / hash_count as f32).log2() + count.unwrap().1));
                            }
                        }
                    }
                    // Sort by count and print
                    let mut query_count_vec: Vec<_> = query_count_map.iter().collect();
                    query_count_vec.sort_by(|a, b| b.1.1.partial_cmp(&a.1.1).unwrap());
                    let count_cut = (query_residues.len() * (query_residues.len() - 1)) / 2;
                    // let count_cut = 0;
                    let hash_set: HashSet<GeometricHash> = pdb_query.iter().cloned().collect();
                    let aa_filter = hash_vec_to_aa_pairs(&pdb_query);
                    for (nid, count) in query_count_vec {
                        if count.0 >= count_cut {
                            if retrieve {
                                let retrieved = retrieve_residue_with_hash(
                                    &hash_set, &aa_filter, &lookup.0[*nid], hash_type, num_bin_dist, num_bin_angle
                                );
                                if retrieved.is_some() {
                                    let retrieved = retrieved.unwrap();
                                    let connected = connected(&retrieved, query_residues.len());
                                    let total_matches = retrieved.len();
                                    if connected > 0 && total_matches < 4 * count.0 {
                                        println!(
                                            "{};uniq_matches={};idf={};total_matches={};connected={};{}",
                                            lookup.0[*nid], count.0, count.1, total_matches,
                                            connected, res_vec_as_string(&retrieved)
                                        );
                                    }
                                }
                            } else {
                                println!("{}\t{}\t{}", lookup.0[*nid], count.0, count.1);
                            }
                        }
                    }
                });
            });
        },
        _ => {
            eprintln!("{}", HELP_QUERY);
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    #[ignore]
    fn test_query_pdb_workflow() {
        // let pdb_path = String::from("data/serine_peptidases_filtered/1azw.pdb");
        // let query_string = String::from("A110,A294,A266");
        let pdb_path = String::from("data/serine_peptidases_filtered/4cha.pdb");
        let query_string = String::from("B57,B102,C195");
        let threads = 1;
        let index_path = Some(String::from("data/serine_peptidases_pdb"));
        let check_nearby = false;
        let retrieve = false;
        let help = false;
        let match_cutoff = 0.0;
        let score_cutoff = 0.0;
        let num_res_cutoff = 3000;
        let plddt_cutoff = 0.0;
        let env = AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            check_nearby,
            retrieve,
            match_cutoff,
            score_cutoff,
            num_res_cutoff,
            plddt_cutoff,
            help,
        };
        query_pdb(env);
    }
}