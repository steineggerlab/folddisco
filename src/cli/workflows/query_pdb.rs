// File: query_pdb.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for querying PDB files
// When querying PDB files, we need index table and query file.

use core::{hash, num};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

use dashmap::DashMap;
use rayon::prelude::*;

use crate::cli::config::{read_index_config_from_file, IndexConfig, QueryConfig};
use crate::controller::mode::IndexMode;
use crate::cli::*;
use crate::controller::io::{
    get_values_with_offset, get_values_with_offset_u16, get_values_with_offset_u24, 
    get_values_with_offset_u8, read_offset_map, read_u16_vector, read_u32_vector, 
    read_u8_vector, read_usize_vector
};
use crate::controller::query::{self, check_and_get_indices, get_offset_value_lookup_type, make_query_map, parse_threshold_string};
use crate::controller::rank::{count_query_gridmode, count_query_idmode, QueryResult};
use crate::controller::retrieve::{connected, hash_vec_to_aa_pairs, res_vec_as_string, retrieve_residue_with_hash};
use crate::index::lookup::{load_lookup_from_file};
use crate::prelude::*;
use crate::structure::grid::{convert_to_id_grid_vector, grid_index_to_tuple, nearby, tuple_to_grid_index};

pub const HELP_QUERY: &str = "\
USAGE: motifsearch query [OPTIONS] <QUERY_PDB> <CHAIN1><RES1>,<CHAIN2><RES2>,<CHAIN3><RES3>...
Example: motifsearch query -i index_table.index -t 6 1aq2.pdb A250,A232,A269
Options:
    -p, --pdb <PDB_PATH>             Path of PDB file to query
    -q, --query <QUERY_STRING>       Query string
    -i, --index <INDEX_PATH>         Path of index table to load
    -t, --threads <THREADS>          Number of threads to use
    -v, --verbose                    Print verbose messages
    -r, --retrieve                   Retrieve matched residues (Need PDB files)
    -d, --distance <DIST_THRESHOLD>  Distance threshold (default 0.0)
    -a, --angle <ANGLE_THRESHOLD>    Angle threshold (default 0.0)
    -m, --match <MATCH_CUTOFF>       Match cutoff (default 0.0)
    -s, --score <SCORE_CUTOFF>       Score cutoff (default 0.0)
    -n, --num-res <NUM_RES_CUTOFF>   Number of residues cutoff (default 3000)
    -l, --plddt <PLDDT_CUTOFF>       PLDDT cutoff (default 0.0)
    -h, --help                       Print this help menu
    --amino-acid <MODE>              Amino acid mode (default 0=matching, 1=all, 2=similar)
";

pub fn query_pdb(env: AppArgs) {
    match env {
        AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            retrieve,
            amino_acid, // TODO:" Implement amino acid mode"
            dist_threshold,
            angle_threshold,
            match_cutoff,
            score_cutoff,
            num_res_cutoff,
            plddt_cutoff,
            verbose,
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
            // Print query information
            if verbose {
                print_log_msg(INFO, &format!("Querying {}:{} to {}", &pdb_path, &query_string, &index_path.clone().unwrap()));
                // NOTE: If needed, print filter information
            }
            // Get index paths
            let index_paths = check_and_get_indices(index_path.clone(), verbose);
            if verbose {
                print_log_msg(INFO, &format!("Found {} index file(s). Querying with {} threads", index_paths.len(), threads));
            }
            // Set thread pool
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();

            let queries = if query_string.ends_with(".txt") || query_string.ends_with(".tsv") {
                // Read file and get path, query, output by line
                let mut queries: Vec<(String, String, String)> = Vec::new();
                let file = std::fs::File::open(&query_string).expect(
                    &log_msg(FAIL, &format!("Failed to open query file: {}", &query_string))
                );
                let reader = std::io::BufReader::new(file);
                for line in reader.lines() {
                    let line = line.expect("Failed to read line");
                    let mut split = line.split('\t');
                    let pdb_path = split.next().expect("Failed to get pdb path").to_string();
                    let query_string = split.next().unwrap_or("").to_string();
                    let output_path = split.next().unwrap_or("").to_string();
                    queries.push((pdb_path, query_string, output_path));
                }
                queries
            } else {
                vec![(pdb_path.clone(), query_string.clone(), String::new())]
            };
            
            let query_map_vec: Vec<(Vec<(usize, QueryResult)>, String, String)> = pool.install(|| {
                index_paths.into_par_iter().map(|index_path| {
                    let (offset_path, value_path, lookup_path, hash_type_path) = get_offset_value_lookup_type(index_path);
                    let config = read_index_config_from_file(&hash_type_path);
                    let hash_type = config.hash_type;
                    let num_bin_dist = config.num_bin_dist;
                    let num_bin_angle = config.num_bin_angle;
                    let mode = config.mode;
                    let dist_thresholds = parse_threshold_string(dist_threshold.clone());
                    let angle_thresholds = parse_threshold_string(angle_threshold.clone());
                    let output: Vec<(Vec<(usize, QueryResult)>, String, String)> = match mode {
                        IndexMode::Id => {
                            let (mmap, value_vec) = measure_time!(
                                read_u16_vector(&value_path).expect(
                                    &log_msg(FAIL, &format!("Failed to load value vector: {}", &value_path))
                                )
                            );
                            let offset_table = measure_time!(
                                read_offset_map(&offset_path, hash_type).expect(
                                    &log_msg(FAIL, &format!("Failed to load offset table: {}", &offset_path))
                                )
                            );
                            let lookup = measure_time!(load_lookup_from_file(&lookup_path));
                            
                            let multi_query_map = measure_time!(
                                parse_multiple_queries(
                                    &queries, hash_type, num_bin_dist, num_bin_angle, &dist_thresholds, &angle_thresholds
                                )
                            );

                            multi_query_map.into_par_iter().map(|(pdb_query_map, pdb_query, query_residues, pdb_path, output_path)| {
                                let query_count_map = measure_time!(count_query_idmode(
                                    &pdb_query, &pdb_query_map ,&offset_table, value_vec, &lookup
                                ));
                                let match_count_filter = get_match_count_filter(match_cutoff.clone(), pdb_query.len(), query_residues.len());
                                let mut query_count_vec: Vec<(usize, QueryResult)> = query_count_map.into_iter().filter(|(k, v)| {
                                    v.total_match_count >= match_count_filter[0] && v.node_count >= match_count_filter[1] && 
                                    v.edge_count >= match_count_filter[2] && v.exact_match_count >= match_count_filter[3] &&
                                    v.overflow_count <= match_count_filter[4] && v.grid_count <= match_count_filter[5] &&
                                    v.idf >= score_cutoff && v.nres <= num_res_cutoff && v.plddt >= plddt_cutoff 
                                }).collect();
                                let query_info_string = format!("{}\t{}", pdb_path, res_chain_to_string(&query_residues));
                                (query_count_vec, query_info_string, output_path)
                            }).collect::<Vec<(Vec<(usize, QueryResult)>, String, String)>>()
                        },
                        IndexMode::Grid => {
                            let (mmap, value_vec) = measure_time!(
                                read_u8_vector(&value_path).expect(
                                    &log_msg(FAIL, &format!("Failed to load value vector: {}", &value_path))
                                )
                            );
                            let offset_table = measure_time!(
                                read_offset_map(&offset_path, hash_type).expect(
                                    &log_msg(FAIL, &format!("Failed to load offset table: {}", &offset_path))
                                )
                            );
                            let lookup = measure_time!(load_lookup_from_file(&lookup_path));
                            // count_query_gridmode(
                            //     pdb_query, pdb_query_map, offset_table, value_vec, &lookup
                            // )
                            let multi_query_map = measure_time!(parse_multiple_queries(
                                &queries, hash_type, num_bin_dist, num_bin_angle, &dist_thresholds, &angle_thresholds
                            ));
                            multi_query_map.into_par_iter().map(|(pdb_query_map, pdb_query, query_residues, pdb_path, output_path)| {
                                let mut query_count_map = measure_time!(count_query_gridmode(
                                    &pdb_query, &pdb_query_map, &offset_table, value_vec, &lookup
                                ));
                                let match_count_filter = get_match_count_filter(
                                    match_cutoff.clone(), pdb_query.len(), query_residues.len()
                                );
                                let mut query_count_vec: Vec<(usize, QueryResult)> = query_count_map.into_iter().filter(|(k, v)| {
                                    v.total_match_count >= match_count_filter[0] && v.node_count >= match_count_filter[1] && 
                                    v.edge_count >= match_count_filter[2] && v.exact_match_count >= match_count_filter[3] &&
                                    v.overflow_count <= match_count_filter[4] && v.grid_count <= match_count_filter[5] &&
                                    v.idf >= score_cutoff && v.nres <= num_res_cutoff && v.plddt >= plddt_cutoff 
                                }).collect();
                                // Convert map to Vec<(usize, QueryResult)>
                                let query_info_string = format!("{}\t{}", pdb_path, res_chain_to_string(&query_residues));
                                (query_count_vec, query_info_string, output_path)
                            }).collect::<Vec<(Vec<(usize, QueryResult)>, String, String)>>()
                        }
                        IndexMode::Pos => {
                            todo!()
                        },
                    };
                    output
                }).reduce(|| Vec::new(), |mut acc, x| {
                    acc.extend(x);
                    acc
                })
            });
            // BEFORE

            // Merge query_map_vec based on same query_info_string
            let mut query_count_map: DashMap<(String, String), Vec<(usize, QueryResult)>> = DashMap::new();
            // for (query_count_vec, query_info_string, output_path) in query_map_vec.iter() {
            query_map_vec.into_par_iter().for_each(|(query_count_vec, query_info_string, output_path)| {
                let key = (query_info_string.clone(), output_path.clone());
                if query_count_map.contains_key(&key) {
                    query_count_map.get_mut(&key).unwrap().extend(query_count_vec.clone());
                } else {
                    query_count_map.insert(key, query_count_vec.clone());
                }
            });
            let mut query_map_vec = query_count_map.into_par_iter().collect::<Vec<((String, String), Vec<(usize, QueryResult)>)>>();
            
            // for ((query_residues, output_path),query_count_vec) in query_map_vec.iter_mut() {
            query_map_vec.par_iter_mut().for_each(|((query_residues, output_path), query_count_vec)| {
                measure_time!(query_count_vec.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap()));
                // If output path is not empty, write to file
                if !output_path.is_empty() {
                    // Create file
                    let file = std::fs::File::create(&output_path).expect(
                        &log_msg(FAIL, &format!("Failed to create file: {}", &output_path))
                    );
                    let mut writer = std::io::BufWriter::new(file);
                    for (k, v) in query_count_vec.iter() {
                        writer.write_all(format!("{:?}\t{}\t{}\n", v, query_residues, index_path.clone().unwrap()
                    ).as_bytes()).expect(
                            &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
                        );
                    }
                } else {
                    for (k, v) in query_count_vec.iter() {
                        println!("{:?}\t{}\t{}", v, query_residues, index_path.clone().unwrap());
                    }
                }
            });
        },
        _ => {
            eprintln!("{}", HELP_QUERY);
            std::process::exit(1);
        }
    }
}

pub fn parse_multiple_queries(
    queries: &Vec<(String, String, String)>,
    hash_type: HashType, num_bin_dist: usize, num_bin_angle: usize, 
    dist_thresholds: &Vec<f32>, angle_thresholds: &Vec<f32>,
) -> Vec<(HashMap<GeometricHash, ((usize, usize), bool)>, Vec<GeometricHash>, Vec<(u8, u64)>, String, String)> {
    let mut query_map_vec = Vec::new();
    for (pdb_path, query_string, output_path) in queries {
        let pdb_loaded = PDBReader::from_file(&pdb_path).expect(
            &log_msg(FAIL, &format!("Failed to read PDB file: {}", &pdb_path))
        );
        let structure = pdb_loaded.read_structure().expect(
            &log_msg(FAIL, &format!("Failed to read the structure from {}", &pdb_path))
        );
        // Make query with pdb
        let query_residues = parse_query_string(&query_string, structure.chains[0]);
        let pdb_query_map =  measure_time!(make_query_map(
            &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, dist_thresholds.clone(), angle_thresholds.clone()
        ));
        let pdb_query: Vec<GeometricHash> = pdb_query_map.keys().cloned().collect();
        query_map_vec.push((pdb_query_map, pdb_query, query_residues, pdb_path.clone(), output_path.clone()));
    }
    query_map_vec
}

pub fn get_match_count_filter(match_cutoff: Option<String>, query_length: usize, residue_count: usize) -> Vec<usize> {
    let match_cutoffs = parse_threshold_string(match_cutoff.clone());
    let mut match_count_filter: Vec<usize> = Vec::with_capacity(match_cutoffs.len());
    for (i, &m) in match_cutoffs.iter().enumerate() {
        if m > 1.0 {
            match_count_filter.push(m as usize);
        } else if m > 0.0 {
            if i == 0 || i == 3 {// i=0: Total match count, i=3: exact match count
                match_count_filter.push((m * query_length as f32).ceil() as usize);
            } else if i == 1 {// i=1: node match count
                match_count_filter.push((m * residue_count as f32).ceil() as usize);
            } else if i == 2 {// i=2: edge match count
                let num_edges = residue_count * (residue_count - 1);
                match_count_filter.push((m * num_edges as f32).ceil() as usize);
            } else  {
                match_count_filter.push(u32::MAX as usize);
            }
        } else {
            if i < 4 {
                match_count_filter.push(0);
            } else {
                match_count_filter.push(u32::MAX as usize);
            }
        }
    }
    // If length is less than 4, fill with 0
    while match_count_filter.len() < 6 {
        if match_count_filter.len() < 4 {
            match_count_filter.push(0);
        } else {
            match_count_filter.push(u32::MAX as usize);
        }
    }
    match_count_filter
}

pub fn res_chain_to_string(res_chain: &Vec<(u8, u64)>) -> String {
    let mut output = String::new();
    for (i, (chain, res)) in res_chain.iter().enumerate() {
        output.push_str(&format!("{}:{}", *chain as char, res));
        if i < res_chain.len() - 1 {
            output.push(',');
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    #[ignore]
    fn test_query_pdb_workflow() {
        let pdb_path = String::from("data/serine_peptidases_filtered/4cha.pdb");
        let query_string = String::from("B57,B102,C195");
        let threads = 4;
        let index_path = Some(String::from("data/serine_peptidases_pdbtr"));
        let exact_match = false;
        let retrieve = false;
        let dist_threshold = Some(String::from("0.5,1.0"));
        let angle_threshold = Some(String::from("5.0,10.0"));
        let help = false;
        let match_cutoff = Some(String::from("0.0,0.8,0.0,0.0"));
        let score_cutoff = 0.0;
        let num_res_cutoff = 3000;
        let plddt_cutoff = 0.0;
        let verbose = true;
        let env = AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            retrieve,
            amino_acid: 0,
            dist_threshold,
            angle_threshold,
            match_cutoff,
            score_cutoff,
            num_res_cutoff,
            plddt_cutoff,
            verbose,
            help,
        };
        query_pdb(env);
    }
    #[test]
    #[ignore]
    fn test_query_pdb_with_file() {
        let pdb_path = String::from("");
        let query_string = String::from("data/query.tsv");
        let threads = 4;
        let index_path = Some(String::from("analysis/h_sapiens/d16a4/index_id"));
        let retrieve = false;
        let dist_threshold = Some(String::from("0.5,1.0"));
        let angle_threshold = Some(String::from("5.0,10.0"));
        let help = false;
        let match_cutoff = Some(String::from("0.0,0.8,0.0,0.0"));
        let score_cutoff = 0.0;
        let num_res_cutoff = 3000;
        let plddt_cutoff = 0.0;
        let verbose = true;
        let env = AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            retrieve,
            amino_acid: 0,
            dist_threshold,
            angle_threshold,
            match_cutoff,
            score_cutoff,
            num_res_cutoff,
            plddt_cutoff,
            verbose,
            help,
        };
        query_pdb(env);
    }
}