// File: query_pdb.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for querying PDB files
// When querying PDB files, we need index table and query file.

use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::PathBuf;

use memmap2::{Mmap, MmapMut};
use rayon::prelude::*;

use crate::cli::config::{read_index_config_from_file, IndexConfig};
use crate::controller::map::SimpleHashMap;
use crate::controller::mode::IndexMode;
use crate::cli::*;
use crate::controller::io::{read_compact_structure, read_u16_vector};
use crate::controller::query::{check_and_get_indices, get_offset_value_lookup_type, make_query_map, parse_threshold_string};
use crate::controller::rank::{count_query_bigmode, count_query_idmode, QueryResult};
use crate::controller::retrieve::retrieval_wrapper;
use crate::index::indextable::load_big_index;
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::*;

#[cfg(feature = "foldcomp")]
use crate::controller::retrieve::retrieval_wrapper_for_foldcompdb;
#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;
#[cfg(feature = "foldcomp")]
use crate::structure::io::StructureFileFormat;

pub const HELP_QUERY: &str = "\
USAGE: folddisco query [OPTIONS] <QUERY_PDB> <CHAIN1><RES1>,<CHAIN2><RES2>,<CHAIN3><RES3>...
Example: folddisco query -i index_table.index -t 6 1aq2.pdb A250,A232,A269
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
    --header                         Print header in output
    --node <NODE_COUNT>              Number of nodes to retrieve (default 2)
";
// TODO: need to think about the column name
pub const QUERY_RESULT_HEADER: &str = "id\tidf_score\ttotal_match_count\tnode_count\tedge_count\texact_match_count\toverflow_count\tnres\tplddt\tmatching_residues\tresidues_rescued\tquery_residues\tquery_file\tindex_path";

pub const DEFAULT_NODE_COUNT: usize = 2;

pub fn query_pdb(env: AppArgs) {
    match env {
        AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            retrieve,
            amino_acid: _, // TODO:" Implement amino acid mode"
            dist_threshold,
            angle_threshold,
            match_cutoff,
            score_cutoff,
            num_res_cutoff,
            plddt_cutoff,
            node_count,
            header,
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
            
            let mut big_index = false;
            if index_paths.len() == 1 {
                // Check index mode is Big
                let (_, _, _, hash_type_path) = get_offset_value_lookup_type(index_paths[0].clone());
                let config = read_index_config_from_file(&hash_type_path);
                if config.mode == IndexMode::Big {
                    big_index = true;
                }
            }
            // Set thread pool
            let _pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            
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
            let dist_thresholds = parse_threshold_string(dist_threshold.clone());
            let angle_thresholds = parse_threshold_string(angle_threshold.clone());
            
            let loaded_index_vec = index_paths.into_par_iter().map(|index_path| {
                let (offset_path, value_path, lookup_path, hash_type_path) = get_offset_value_lookup_type(index_path);
                let config = read_index_config_from_file(&hash_type_path);
                let (offset_table, offset_mmap) = if verbose && !big_index { measure_time!(
                    SimpleHashMap::load_from_disk(&PathBuf::from(&offset_path))
                ) } else if !big_index {
                    SimpleHashMap::load_from_disk(&PathBuf::from(&offset_path))
                } else {
                    let anon_mmap: Mmap = MmapMut::map_anon(0).unwrap().make_read_only().unwrap();
                    (Ok(SimpleHashMap::new(0)), anon_mmap)
                };
                let offset_table = offset_table.unwrap();
                let lookup = if verbose {
                    measure_time!(load_lookup_from_file(&lookup_path))
                } else {
                    load_lookup_from_file(&lookup_path)
                };
                (offset_table, offset_mmap, lookup, config, value_path)
            }).collect::<Vec<(SimpleHashMap, Mmap, (Vec<String>, Vec<usize>, Vec<usize>, Vec<f32>), IndexConfig, String)>>();

            // Load foldcomp db 
            #[cfg(feature = "foldcomp")]
            let config = &loaded_index_vec[0].3;
            #[cfg(feature = "foldcomp")]
            let using_foldcomp = config.foldcomp_db.is_some() && config.input_format == StructureFileFormat::FCZDB;

            #[cfg(feature = "foldcomp")]
            let foldcomp_db_reader = match config.input_format {
                StructureFileFormat::FCZDB => {
                    if retrieve {
                        let foldcomp_db_path = config.foldcomp_db.clone().unwrap();
                        FoldcompDbReader::new(foldcomp_db_path.as_str())
                    } else {
                        FoldcompDbReader::empty()
                    }
                },
                _ => FoldcompDbReader::empty(),
            };

            // #[cfg(not(feature = "foldcomp"))]
            // let using_foldcomp = false;

            // Iterate over queries
            queries.into_par_iter().for_each(|(pdb_path, query_string, output_path)| {
                let (query_structure, _) = read_compact_structure(&pdb_path).expect(
                    &log_msg(FAIL, &format!("Failed to read structure: {}", &pdb_path))
                );
                
                let (query_residues, aa_substitutions) = parse_query_string(&query_string, query_structure.chains[0]);
                
                let _residue_count = if query_residues.is_empty() {
                    query_structure.num_residues
                } else {
                    query_residues.len()
                };
                let query_string = if query_residues.is_empty() {
                    query_string
                } else {
                    let mut query_residues = query_residues.clone();
                    query_residues.sort();
                    res_chain_to_string(&query_residues)
                };

                // Get query map for each query in all indices
                let mut queried_from_indices: Vec<(usize, QueryResult)> = loaded_index_vec.par_iter().map(
                    |(offset_table, _offset_mmap, lookup, config, value_path)| {
                        let hash_type = config.hash_type;
                        let num_bin_dist = config.num_bin_dist;
                        let num_bin_angle = config.num_bin_angle;
                        let mode = config.mode;
                        let dist_cutoff = config.grid_width;
                        let (pdb_query_map, query_indices, aa_dist_map ) = if verbose { measure_time!(make_query_map(
                            &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, &dist_thresholds, &angle_thresholds,  &aa_substitutions, dist_cutoff,
                        )) } else {
                            make_query_map(
                                &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, 
                                &dist_thresholds, &angle_thresholds, &aa_substitutions, dist_cutoff
                            )
                        };
                        let pdb_query = pdb_query_map.keys().cloned().collect::<Vec<_>>();
                        match mode {
                            IndexMode::Id => {
                                let (value_mmap, value_vec) = if verbose { measure_time!(read_u16_vector(&value_path).expect(
                                    &log_msg(FAIL, &format!("Failed to load value vector: {}", &value_path))
                                )) } else {
                                    read_u16_vector(&value_path).expect(&log_msg(FAIL, &format!("Failed to load value vector: {}", &value_path)))
                                };
                                let query_count_map = if verbose { measure_time!(count_query_idmode(
                                    &pdb_query, &pdb_query_map, &offset_table, value_vec, &lookup
                                ))} else {
                                    count_query_idmode(&pdb_query, &pdb_query_map, &offset_table, value_vec, &lookup)
                                };
                                let mut match_count_filter = get_match_count_filter(
                                    match_cutoff.clone(), pdb_query.len(), query_residues.len()
                                );
                                if retrieve {
                                    match_count_filter[1] = if node_count > match_count_filter[1] {node_count} else {match_count_filter[1]};
                                    if verbose {
                                        print_log_msg(INFO, &format!("Filtering result with residue >= {}", match_count_filter[1]));
                                    }
                                }
                                let mut query_count_vec: Vec<(usize, QueryResult)> = query_count_map.into_par_iter().filter(|(_k, v)| {
                                    v.total_match_count >= match_count_filter[0] && v.node_count >= match_count_filter[1] && 
                                    v.edge_count >= match_count_filter[2] && v.exact_match_count >= match_count_filter[3] &&
                                    v.overflow_count <= match_count_filter[4] &&
                                    v.idf >= score_cutoff && v.nres <= num_res_cutoff && v.plddt >= plddt_cutoff 
                                }).collect();
                                if verbose {
                                    print_log_msg(INFO, &format!("Found {} structures from inverted index", query_count_vec.len()));
                                }

                                // IF retrieve is true, retrieve matching residues
                                if retrieve {
                                    measure_time!(query_count_vec.par_iter_mut().for_each(|(_, v)| {
                                        #[cfg(not(feature = "foldcomp"))]
                                        let retrieval_result = retrieval_wrapper(
                                            &v.id, node_count, &pdb_query,
                                            hash_type, num_bin_dist, num_bin_angle, dist_cutoff,
                                            &pdb_query_map, &query_structure, &query_indices,
                                            &aa_dist_map
                                        );
                                        #[cfg(feature = "foldcomp")]
                                        let retrieval_result = if using_foldcomp {
                                            retrieval_wrapper_for_foldcompdb(
                                                &v.id, node_count, &pdb_query,
                                                hash_type, num_bin_dist, num_bin_angle, dist_cutoff,
                                                &pdb_query_map, &query_structure, &query_indices,
                                                &aa_dist_map, &foldcomp_db_reader
                                            )
                                        } else {
                                            retrieval_wrapper(
                                                &v.id, node_count, &pdb_query,
                                                hash_type, num_bin_dist, num_bin_angle, dist_cutoff,
                                                &pdb_query_map, &query_structure, &query_indices,
                                                &aa_dist_map
                                            )
                                        };
                                        v.matching_residues = retrieval_result.0;
                                        v.matching_residues_processed = retrieval_result.1;
                                    }));
                                    // Filter query_count_vec with reasonable retrieval results
                                    query_count_vec.retain(|(_, v)| v.matching_residues.len() > 0);
                                    drop(value_mmap);
                                    drop(match_count_filter);
                                    return query_count_vec;
                                }
                                drop(value_mmap);
                                drop(match_count_filter);
                                query_count_vec
                            },
                            IndexMode::Big => {
                                // Get index prefix from value path. index_prefix is just without .value
                                let index_prefix = value_path.rsplitn(2, '.').collect::<Vec<&str>>()[1];
                                let (big_index, big_offset_mmap) = if verbose {
                                    load_big_index(index_prefix)
                                } else {
                                    load_big_index(index_prefix)
                                };
                                let query_count_map = measure_time!(count_query_bigmode(
                                    &pdb_query, &pdb_query_map, &big_index, &lookup
                                ));
                                
                                let mut match_count_filter = get_match_count_filter(
                                    match_cutoff.clone(), pdb_query.len(), query_residues.len()
                                );
                                if retrieve {
                                    match_count_filter[1] = if node_count > match_count_filter[1] {node_count} else {match_count_filter[1]};
                                    if verbose {
                                        print_log_msg(INFO, &format!("Filtering result with residue >= {}", match_count_filter[1]));
                                    }
                                }
                                let mut query_count_vec: Vec<(usize, QueryResult)> = query_count_map.into_par_iter().filter(|(_k, v)| {
                                    v.total_match_count >= match_count_filter[0] && v.node_count >= match_count_filter[1] && 
                                    v.edge_count >= match_count_filter[2] && v.exact_match_count >= match_count_filter[3] &&
                                    v.overflow_count <= match_count_filter[4] && 
                                    v.idf >= score_cutoff && v.nres <= num_res_cutoff && v.plddt >= plddt_cutoff 
                                }).collect();
                                if verbose {
                                    print_log_msg(INFO, &format!("Found {} structures from inverted index", query_count_vec.len()));
                                }

                                // IF retrieve is true, retrieve matching residues
                                if retrieve {
                                    measure_time!(query_count_vec.par_iter_mut().for_each(|(_, v)| {
                                        #[cfg(not(feature = "foldcomp"))]
                                        let retrieval_result = retrieval_wrapper(
                                            &v.id, node_count, &pdb_query,
                                            hash_type, num_bin_dist, num_bin_angle, dist_cutoff,
                                            &pdb_query_map, &query_structure, &query_indices,
                                            &aa_dist_map
                                        );
                                        #[cfg(feature = "foldcomp")]
                                        let retrieval_result = if using_foldcomp {
                                            retrieval_wrapper_for_foldcompdb(
                                                &v.id, node_count, &pdb_query,
                                                hash_type, num_bin_dist, num_bin_angle, dist_cutoff,
                                                &pdb_query_map, &query_structure, &query_indices,
                                                &aa_dist_map, &foldcomp_db_reader
                                            )
                                        } else {
                                            retrieval_wrapper(
                                                &v.id, node_count, &pdb_query,
                                                hash_type, num_bin_dist, num_bin_angle, dist_cutoff,
                                                &pdb_query_map, &query_structure, &query_indices,
                                                &aa_dist_map
                                            )
                                        };
                                        v.matching_residues = retrieval_result.0;
                                        v.matching_residues_processed = retrieval_result.1;
                                    }));
                                    // Filter query_count_vec with reasonable retrieval results
                                    query_count_vec.retain(|(_, v)| v.matching_residues.len() > 0);
                                    drop(big_offset_mmap);
                                    drop(match_count_filter);
                                    return query_count_vec;
                                }
                                drop(big_offset_mmap);
                                drop(match_count_filter);
                                query_count_vec
                            },
                        } // match mode
                    }
                ).reduce(|| Vec::new(), |mut acc, x| {
                    acc.extend(x);
                    acc
                });
                drop(query_residues);
                // Sort query_count_vec by idf
                // TODO: Change feature to sort
                if verbose {
                    measure_time!(queried_from_indices.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap()));
                } else {
                    queried_from_indices.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap());
                }
                    // If output path is not empty, write to file
                if !output_path.is_empty() {
                    // Create file
                    let file = std::fs::File::create(&output_path).expect(
                        &log_msg(FAIL, &format!("Failed to create file: {}", &output_path))
                    );
                    let mut writer = std::io::BufWriter::new(file);
                    if header {
                        writer.write_all(format!("{}\n", QUERY_RESULT_HEADER).as_bytes()).expect(
                            &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
                        );
                    }
                    for (_k, v) in queried_from_indices.iter() {
                        writer.write_all(format!("{:?}\t{}\t{}\t{}\n", v, query_string, pdb_path, index_path.clone().unwrap()).as_bytes()).expect(
                            &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
                        );
                    }
                } else {
                    if header {
                        println!("{}", QUERY_RESULT_HEADER);
                    }
                    for (_k, v) in queried_from_indices.iter() {
                        println!("{:?}\t{}\t{}\t{}", v, query_string, pdb_path, index_path.clone().unwrap());
                    }
                }
            }); // queries
        }, // AppArgs::Query
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
        let (query_residues, aa_substitions) = parse_query_string(&query_string, structure.chains[0]);
        let (pdb_query_map, _query_indices, _) =  measure_time!(make_query_map(
            &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, &dist_thresholds, &angle_thresholds, &aa_substitions, 20.0
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
        output.push_str(&format!("{}{}", *chain as char, res));
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
        let threads = 1;
        let index_path = Some(String::from("data/serine_peptidases_pdbtr_small"));
        // let index_path = Some(String::from("analysis/e_coli/test_index"));
        let retrieve = true;
        let dist_threshold = Some(String::from("0.5,1.0"));
        let angle_threshold = Some(String::from("5.0,10.0"));
        let help = false;
        let match_cutoff = Some(String::from("0.0"));
        let score_cutoff = 0.0;
        let num_res_cutoff = 3000;
        let plddt_cutoff = 0.0;
        let node_count = 2;
        let header = false;
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
            node_count,
            header,
            verbose,
            help,
        };
        query_pdb(env);
    }
    #[test]
    #[ignore]
    fn test_query_with_foldcompdb() {
        #[cfg(feature = "foldcomp")] {
            let pdb_path = String::from("data/foldcomp/example_db:d1asha_");
            let query_string = String::from("1,2,3,4");
            let threads = 1;
            let index_path = Some(String::from("data/example_db_folddisco_db"));
            // let index_path = Some(String::from("analysis/e_coli/test_index"));
            let retrieve = true;
            let dist_threshold = Some(String::from("0.5,1.0"));
            let angle_threshold = Some(String::from("5.0,10.0"));
            let help = false;
            let match_cutoff = Some(String::from("0.0"));
            let score_cutoff = 0.0;
            let num_res_cutoff = 3000;
            let plddt_cutoff = 0.0;
            let node_count = 2;
            let header = false;
            let verbose = false;
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
                node_count,
                header,
                verbose,
                help,
            };
            query_pdb(env);
        }
    }
    #[test]
    #[ignore]
    fn test_query_pdb_with_file() {
        let pdb_path = String::from("");
        let query_string = String::from("data/query.tsv");
        let threads = 4;
        let index_path = Some(String::from("analysis/e_coli/test_index"));
        let retrieve = true;
        let dist_threshold = Some(String::from("0.5,1.0"));
        let angle_threshold = Some(String::from("5.0,10.0"));
        let help = false;
        let match_cutoff = Some(String::from("0.0,0.8,0.0,0.0"));
        let score_cutoff = 0.0;
        let num_res_cutoff = 3000;
        let plddt_cutoff = 0.0;
        let node_count = 2;
        let header = true;
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
            node_count,
            header,
            verbose,
            help,
        };
        query_pdb(env);
    }
}