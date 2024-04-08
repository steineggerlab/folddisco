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

use crate::cli::config::{read_index_config_from_file, IndexConfig, QueryConfig};
use crate::controller::mode::IndexMode;
use crate::cli::*;
use crate::controller::io::{
    get_values_with_offset, get_values_with_offset_u16, get_values_with_offset_u24, 
    get_values_with_offset_u8, read_offset_map, read_u16_vector, read_u32_vector, 
    read_u8_vector, read_usize_vector
};
use crate::controller::query::{check_and_get_indices, get_offset_value_lookup_type, make_query_map, parse_threshold_string};
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
            // Get index paths
            let index_paths = check_and_get_indices(index_path, verbose);
            if verbose {
                print_log_msg(INFO, &format!("Found {} index file(s). Querying with {} threads", index_paths.len(), threads));
            }
            // Set thread pool
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();

            let query_map_vec: HashMap<usize, QueryResult> = pool.install(|| {
                index_paths.into_par_iter().map(|index_path| {
                    let (offset_path, value_path, lookup_path, hash_type_path) = get_offset_value_lookup_type(index_path);
                    let config = read_index_config_from_file(&hash_type_path);
                    let hash_type = config.hash_type;
                    let num_bin_dist = config.num_bin_dist;
                    let num_bin_angle = config.num_bin_angle;
                    let mode = config.mode;
                    
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

                    let pdb_query_map =  measure_time!(make_query_map(
                        &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, dist_thresholds, angle_thresholds
                    ));
                    let pdb_query: Vec<GeometricHash> = pdb_query_map.keys().cloned().collect();

                    let match_cutoffs = parse_threshold_string(match_cutoff.clone());
                    let mut match_count_filter: Vec<usize> = Vec::with_capacity(match_cutoffs.len());
                    for (i, &m) in match_cutoffs.iter().enumerate() {
                        if m > 1.0 {
                            match_count_filter.push(m as usize);
                        } else if m > 0.0 {
                            if i == 0 || i == 3 {// i=0: Total match count, i=3: exact match count
                                match_count_filter.push((m * pdb_query.len() as f32).ceil() as usize);
                            } else if i == 1 {// i=1: node match count
                                match_count_filter.push((m * query_residues.len() as f32).ceil() as usize);
                            } else if i == 2 {// i=2: edge match count
                                let num_edges = query_residues.len() * (query_residues.len() - 1);
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

                    // Load index table
                    let query_count_map: HashMap<usize, QueryResult> = match mode {
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
                            count_query_idmode(
                                pdb_query, pdb_query_map, offset_table, value_vec, &lookup
                            )
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
                            count_query_gridmode(
                                pdb_query, pdb_query_map, offset_table, value_vec, &lookup
                            )
                        }
                        IndexMode::Pos => {
                            let (mmap, value_vec) = measure_time!(read_u32_vector(&value_path).expect("[ERROR] Failed to load value vector"));
                            let offset_table = measure_time!(
                                read_offset_map(&offset_path, hash_type).expect("[ERROR] Failed to load offset table")
                            );
                            let lookup = measure_time!(load_lookup_from_file(&lookup_path));
                            let mut query_count_map: HashMap<usize, QueryResult> = HashMap::new();
                            // todo!("IMPLEMENT NEEDED!");
                            query_count_map
                        },
                    };
                    // Filter query count map
                    query_count_map.into_iter().filter(|(k, v)| {
                        v.total_match_count >= match_count_filter[0] && v.node_count >= match_count_filter[1] && 
                        v.edge_count >= match_count_filter[2] && v.exact_match_count >= match_count_filter[3] &&
                        v.overflow_count <= match_count_filter[4] && v.grid_count <= match_count_filter[5] &&
                        v.idf >= score_cutoff && v.nres <= num_res_cutoff && v.plddt >= plddt_cutoff 
                    }).collect()
                }).reduce(|| HashMap::new(), |mut acc, x| {
                    acc.extend(x);
                    acc
                })
            });
            
            let mut query_count_vec: Vec<(usize, QueryResult)> = query_map_vec.into_iter().collect();
            query_count_vec.sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap());
            for (k, v) in query_count_vec.iter() {
                println!("{:?}", v);
            }
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
        let threads = 4;
        let index_path = Some(String::from("data/serine_peptidases_3di"));
        let exact_match = false;
        let retrieve = false;
        let dist_threshold = Some(String::from("0.5,1.0"));
        let angle_threshold = Some(String::from("5.0,10.0"));
        let help = false;
        let match_cutoff = Some(String::from("0.0,0.0,0.0,0.0"));
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