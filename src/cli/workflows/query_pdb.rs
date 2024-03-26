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
use crate::controller::query::{make_query_map, parse_threshold_string};
use crate::controller::retrieve::{connected, hash_vec_to_aa_pairs, res_vec_as_string, retrieve_residue_with_hash};
use crate::index::lookup::{load_lookup_from_file};
use crate::prelude::*;


pub const HELP_QUERY: &str = "\
USAGE: motifsearch query [OPTIONS] <QUERY_PDB> <CHAIN1><RES1>,<CHAIN2><RES2>,<CHAIN3><RES3>...
Example: motifsearch query -i index_table.index -t 6 1aq2.pdb A250,A232,A269
Options:
    -d, --pdb <PDB_PATH>             Path of PDB file to query
    -q, --query <QUERY_STRING>       Query string
    -i, --index <INDEX_PATH>         Path of index table to load
    -t, --threads <THREADS>          Number of threads to use
    -v, --verbose                    Print verbose messages
    -e, --exact                      Return only exact matches
    -v, --retrieve                   Retrieve matched residues (Need PDB files)
    -d, --distance <DIST_CUTOFF>     Distance cutoff (default 0.0)
    -a, --angle <ANGLE_CUTOFF>       Angle cutoff (default 0.0)
    -m, --match <MATCH_CUTOFF>       Match cutoff (default 0.0)
    -s, --score <SCORE_CUTOFF>       Score cutoff (default 0.0)
    -n, --num-res <NUM_RES_CUTOFF>   Number of residues cutoff (default 3000)
    -l, --plddt <PLDDT_CUTOFF>       PLDDT cutoff (default 0.0)
    -h, --help                       Print this help menu
";


pub fn query_pdb(env: AppArgs) {
    match env {
        AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            exact_match,
            retrieve,
            dist_threshold,
            angle_threshold,
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
                    
                    let pdb_loaded = PDBReader::from_file(&pdb_path).expect(
                        &log_msg(FAIL, &format!("Failed to read PDB file: {}", &pdb_path))
                    );
                    let structure = pdb_loaded.read_structure().expect(
                        &log_msg(FAIL, &format!("Failed to read the structure from {}", &pdb_path))
                    );
                    let dist_thresholds = parse_threshold_string(dist_threshold.clone());
                    let angle_thresholds = parse_threshold_string(angle_threshold.clone());
                    // Make query with pdb
                    let query_residues = parse_query_string(&query_string, structure.chains[0]);
                    // let pdb_query =  measure_time!(make_query(
                    //     &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, exact_match, dist_thresholds, angle_thresholds
                    // ));
                    // NOTE: testing with new querying scheme that checks the node coverage
                    let pdb_query_map =  measure_time!(make_query_map(
                        &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, dist_thresholds, angle_thresholds
                    ));
                    let pdb_query: Vec<GeometricHash> = pdb_query_map.keys().cloned().collect();

                    // Load index table
                    // let (mmap, value_vec) = measure_time!(read_usize_vector(&value_path).expect("[ERROR] Failed to load value vector"));
                    let (mmap, value_vec) = measure_time!(read_u16_vector(&value_path).expect("[ERROR] Failed to load value vector"));
                    
                    let offset_table = measure_time!(
                        read_offset_map(&offset_path, hash_type).expect("[ERROR] Failed to load offset table")
                    );
                    // let (path_vec, numeric_id_vec, optional_vec) = measure_time!(load_lookup_from_file(&lookup_path));
                    let lookup = measure_time!(load_lookup_from_file(&lookup_path));
                    
                    // Make union of values queried
                    // Query count map: nid -> (count, idf, nres, plddt, node_set, edge_set, exact_match_count, overflow_count)
                    let mut query_count_map: HashMap<usize, (usize, f32, usize, f32, HashSet<usize>, HashSet<(usize, usize)>, usize, usize)> = HashMap::new();
                    for i in 0..pdb_query.len() {
                        let offset = offset_table.get(&pdb_query[i]);
                        if offset.is_none() {
                            continue;
                        }
                        let offset = *offset.unwrap();
                        let edge_info = pdb_query_map.get(&pdb_query[i]).unwrap();
                        let is_exact = edge_info.1;
                        let edge = edge_info.0;
                        // let single_queried_values = get_values_with_offset(&value_vec, offset_to_query[i].0, offset_to_query[i].1);
                        let single_queried_values = get_values_with_offset_u16(&value_vec, offset.0, offset.1);
                        let hash_count = offset.1;
                        for j in 0..single_queried_values.len() {
                            let nid = lookup.1[single_queried_values[j] as usize];
                            let nres = lookup.2[single_queried_values[j] as usize];
                            let plddt = lookup.3[single_queried_values[j] as usize];
                            
                            if nres > num_res_cutoff || plddt < plddt_cutoff {
                                continue;
                            }
                            let count = query_count_map.get(&nid);
                            // Added calculation of idf
                            let idf = (lookup.0.len() as f32 / hash_count as f32).log2();
                            // TODO: Check if normalization is done per match or per query
                            let nres_norm = (nres as f32).log2() * -0.5 + 6.0;
                            if count.is_none() {
                                let mut node_set = HashSet::new();
                                node_set.insert(edge.0);
                                node_set.insert(edge.1);
                                let mut edge_set = HashSet::new();
                                edge_set.insert(edge);
                                let exact_match_count = if is_exact { 1usize } else { 0usize };
                                let overflow_count = 0usize;
                                query_count_map.insert(nid, (1, idf + nres_norm, nres, plddt, node_set, edge_set, exact_match_count, overflow_count));
                            } else {
                                // Update node set and edge set
                                let mut node_set = count.unwrap().4.clone();
                                node_set.insert(edge.0);
                                node_set.insert(edge.1);
                                let mut edge_set = count.unwrap().5.clone();
                                // If edge is already in the set, it is an overflow. Check edge set contains the edge
                                let is_overflow: bool = edge_set.contains(&edge);
                                edge_set.insert(edge);
                                let exact_match_count = if is_exact { count.unwrap().6 + 1 } else { count.unwrap().6 };
                                let overflow_count = if is_overflow { count.unwrap().7 + 1 } else { count.unwrap().7 };
                                let idf = if is_exact { idf + nres_norm } else { (idf + nres_norm) / (1.5 + (overflow_count * overflow_count) as f32) };
                                // query_count_map.insert(nid, (count.unwrap().0 + 1, idf + count.unwrap().1 + nres_norm, nres, plddt));
                                query_count_map.insert(nid, (count.unwrap().0 + 1, idf + count.unwrap().1, nres, plddt, node_set, edge_set, exact_match_count, overflow_count));
                            }
                        }
                    }
                    // Sort by count and print
                    let mut query_count_vec: Vec<_> = query_count_map.iter().collect();
                    // Nested sorting by idf & exact match count
                    query_count_vec.sort_by(|a, b| b.1.1.partial_cmp(&a.1.1).unwrap());

                    let count_cut = match_cutoff.round() as usize;
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
                                    println!(
                                        "{};uniq_matches={};idf={};total_matches={};connected={};{}",
                                        lookup.0[*nid], count.0, count.1, total_matches,
                                        connected, res_vec_as_string(&retrieved)
                                    );
                                }
                            } else {
                                if count.1 < score_cutoff {
                                    continue;
                                }
                                println!(
                                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                                    lookup.0[*nid], count.0, count.1, count.2, count.3,
                                    count.4.len(), count.5.len(), count.6, count.7
                                );
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
        let exact_match = false;
        let retrieve = false;
        let dist_threshold = Some(String::from("0.5"));
        let angle_threshold = Some(String::from("5.0,10.0"));
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
            exact_match: exact_match,
            retrieve,
            dist_threshold,
            angle_threshold,
            match_cutoff,
            score_cutoff,
            num_res_cutoff,
            plddt_cutoff,
            help,
        };
        query_pdb(env);
    }
}