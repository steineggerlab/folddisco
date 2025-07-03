// File: query_pdb.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for querying PDB files
// When querying PDB files, we need index table and query file.

use std::io::BufRead;
use std::path::PathBuf;

use memmap2::{Mmap, MmapMut};
use rayon::prelude::*;

use crate::cli::config::{read_index_config_from_file, IndexConfig};
use crate::controller::filter::{MatchFilter, StructureFilter};
use crate::controller::map::SimpleHashMap;
use crate::controller::mode::{IndexMode, QueryMode};
use crate::cli::*;
use crate::controller::io::{read_compact_structure, read_u16_vector};
use crate::controller::query::{
    check_and_get_indices, get_offset_value_lookup_type, make_query_map, parse_threshold_string
};
use crate::controller::count_query::{count_query_bigmode, count_query_idmode};
use crate::controller::result::{
    convert_structure_query_result_to_match_query_results, 
    sort_and_print_match_query_result, sort_and_print_structure_query_result, StructureResult
};
use crate::controller::retrieve::retrieval_wrapper;
use crate::index::indextable::{load_big_index, FolddiscoIndex};
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::*;

#[cfg(feature = "foldcomp")]
use crate::controller::retrieve::retrieval_wrapper_for_foldcompdb;
#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;
#[cfg(feature = "foldcomp")]
use crate::structure::io::StructureFileFormat;

pub const HELP_QUERY: &str = "\
usage: folddisco query -p <i:PDB> -q <QUERY> -i <i:INDEX> [OPTIONS] 

input/output:
 -p, --pdb <PATH>                 Path of PDB file to query
 -q, --query <STR>                Query string that specifies residues or a text file containing query
 -i, --index <PATH>               Path of index table to load [REQUIRED]
 -o, --output <PATH>              Output file path [stdout]
 
search parameters:
 -t, --threads <INT>              Number of threads [1]
 -d, --distance <FLOAT>           Distance threshold in Angstroms. Multiple values can be separated by comma [0.0]
 -a, --angle <FLOAT>              Angle threshold. Multiple values can be separated by comma [0.0]
 --ca-distance <FLOAT>            C-alpha distance threshold in matching residues [1.5]
 --sampling-count <INT>           Number of sampled hashes to search [all]
 --sampling-ratio <FLOAT>         Sampling ratio for hashes used in searching. For long queries, smaller ratio is recommended [1.0]
 --freq-filter <FLOAT>            Skip queries with hash frequency higher than given ratio [0.0]
 --length-penalty <FLOAT>         Length penalty for searching. Zero means no penalty and higher value gives more penalty to longer structures [0.5]
 --skip-match                     Skip matching residues
 --serial-index                   Handle residue indices serially

filtering options:
 --total-match <INT>              Filter out structures with less than total match count [0]
 --covered-node <INT>             Filter out structures not covered by given number of nodes with hashes [0]
 --covered-node-ratio <FLOAT>     Filter out structures not covered by given ratio of nodes with hashes [0.0]
 --covered-edge <INT>             Filter out structures not covered by given number of edges with hashes [0]
 --covered-edge-ratio <FLOAT>     Filter out structures not covered by given ratio of edges with hashes [0.0]
 --max-node <INT>                 Filter out structures of maximum matching node size smaller than given value [0]
 --max-node-ratio <FLOAT>         Filter out structures of maximum matching node size smaller than given ratio [0.0]
 --score <FLOAT>                  IDF score cutoff [0.0]
 --connected-node <INT>           Filter out structures/matches with connected node count smaller than given value [0]
 --connected-node-ratio <FLOAT>   Filter out structures/matches with connected node count smaller than given ratio [0.0]
 --num-residue <INT>              Number of residues cutoff [50000]
 --plddt <FLOAT>                  pLDDT cutoff [0.0]
 --top <INT>                      Limit output to top N structures [all]

display options:
 --header                         Print header in output
 --web                            Print output for web
 --per-structure                  Print output per structure
 --per-match                      Print output per match. Not working with --skip-match
 --sort-by-score                  Sort output by score
 --sort-by-rmsd                   Sort output by RMSD. Not working with --skip-match
 --skip-ca-match                  Print matching residues before C-alpha distance check
 --superpose                      Print U, T, CA of matching residues

general options:
 -v, --verbose                    Print verbose messages
 -h, --help                       Print this help menu

examples:
# Search with default settings. This will print out matching motifs with sorting by RMSD.
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6

# Print top 100 structures with sorting by score
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --top 100 --per-structure --sort-by-score

# Query file given as separate text file
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 -d 0.5 -a 5

# Query with amino-acid substitutions and range. 
# Alternative amino acids can be given after colon. Range can be given with dash.
# This will query first 10 residues and 11th residue with subsitution to any amino acid.
folddisco query -p query/4CHA.pdb -q 1-10,11:X -i index/h_sapiens_folddisco -t 6 --serial-index

# Filtering
## Based on connected node and rmsd
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --connected-node 0.75 --rmsd 1.0

## Coverage based filtering & top N filtering without residue matching
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --covered-node 3 --top 1000 --per-structure --skip-match
";

pub const MIN_CONNECTED_COMPONENT_SIZE: usize = 2;
pub const MAX_NUM_LINES_FOR_WEB: usize = 1000;

pub fn query_pdb(env: AppArgs) {
    match env {
        AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            skip_match,
            dist_threshold,
            angle_threshold,
            ca_dist_threshold,
            total_match_count,
            covered_node_count,
            covered_node_ratio,
            covered_edge_count,
            covered_edge_ratio,
            max_matching_node_count,
            max_matching_node_ratio,
            idf_score_cutoff,
            connected_node_count,
            connected_node_ratio,
            num_res_cutoff,
            plddt_cutoff,
            rmsd_cutoff,
            top_n,
            web_mode,
            sampling_count,
            sampling_ratio,
            freq_filter,
            length_penalty,
            sort_by_rmsd,
            sort_by_score,
            output_per_structure,
            output_per_match,
            output_with_superpose,
            skip_ca_match,
            header,
            serial_query,
            output,
            verbose,
            help: _,
        } => {
            if verbose { print_logo(); }
            // help is already handled in main.rs
            // Check if arguments are valid
            if index_path.is_none() {
                eprintln!("{}", HELP_QUERY);
                std::process::exit(1);
            }
            
            let query_mode = QueryMode::from_flags(
                skip_match, web_mode, output_per_structure, output_per_match,
                sort_by_rmsd, sort_by_score
            );
            // Error handling
            match query_mode {
                QueryMode::ContradictoryPrintError => {
                    print_log_msg(FAIL, 
                        "Cannot print output per structure and per match at the same time. Use either --per-structure or --per-match"
                    );
                    std::process::exit(1);
                }
                QueryMode::ContradictorySortError => {
                    print_log_msg(FAIL, 
                        "Cannot sort result by RMSD and score at the same time. Use either --sort-by-rmsd or --sort-by-score"
                    );
                }
                _ => {}
            }

            // Print query information
            if verbose {
                // If pdb_path is empty
                if pdb_path.is_empty() {
                    print_log_msg(INFO, &format!("Querying {} to {}", &query_string, &index_path.clone().unwrap()));
                } else {
                    print_log_msg(INFO, &format!("Querying {}:{} to {}", &pdb_path, &query_string, &index_path.clone().unwrap()));
                }
                // NOTE: If needed, print filter information
            }
            // Get index paths
            let index_paths = check_and_get_indices(index_path.clone(), verbose);
            if verbose {
                print_log_msg(INFO, &format!("Found {} index file(s). Querying with {} threads", index_paths.len(), threads));
            }
            
            let mut use_big_index = false;
            if index_paths.len() == 1 {
                // Check index mode is Big
                let (_, _, _, hash_type_path) = get_offset_value_lookup_type(index_paths[0].clone());
                let config = read_index_config_from_file(&hash_type_path);
                if config.mode == IndexMode::Big {
                    use_big_index = true;
                }
            }
            // Load big index here
            
            let index_prefix = index_paths[0].clone();
            let (big_index, big_offset_mmap) = if use_big_index {
                load_big_index(&index_prefix)
            } else {
                (FolddiscoIndex::new(0, "".to_string(), false),
                 MmapMut::map_anon(0).unwrap().make_read_only().unwrap())
            };

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
                vec![(pdb_path.clone(), query_string.clone(), output.clone())]
            };

            let dist_thresholds = parse_threshold_string(Some(dist_threshold.clone()));
            let angle_thresholds = parse_threshold_string(Some(angle_threshold.clone()));
    
            let loaded_index_vec = index_paths.into_par_iter().map(|index_path| {
                let (offset_path, value_path, lookup_path, hash_type_path) = get_offset_value_lookup_type(index_path);
                let config = read_index_config_from_file(&hash_type_path);
                let (offset_table, offset_mmap) = if verbose && !use_big_index { measure_time!(
                    SimpleHashMap::load_from_disk(&PathBuf::from(&offset_path))
                ) } else if !use_big_index {
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
            }).collect::<Vec<(SimpleHashMap, Mmap, Vec<(String, usize, usize, f32, usize)>, IndexConfig, String)>>();
            
            // Load foldcomp db 
            #[cfg(feature = "foldcomp")]
            let config = &loaded_index_vec[0].3;
            #[cfg(feature = "foldcomp")]
            let using_foldcomp = config.foldcomp_db.is_some() && config.input_format == StructureFileFormat::FCZDB;

            #[cfg(feature = "foldcomp")]
            let foldcomp_db_reader = match config.input_format {
                StructureFileFormat::FCZDB => {
                    if !skip_match {
                        let foldcomp_db_path = config.foldcomp_db.clone().unwrap();
                        if verbose {
                            measure_time!(FoldcompDbReader::new(foldcomp_db_path.as_str()))
                        } else {
                            FoldcompDbReader::new(foldcomp_db_path.as_str())
                        }
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
                let mut queried_from_indices: Vec<(usize, StructureResult)> = loaded_index_vec.par_iter().map(
                    |(offset_table, _offset_mmap, lookup, config, value_path)| {
                        let hash_type = config.hash_type;
                        let num_bin_dist = config.num_bin_dist;
                        let num_bin_angle = config.num_bin_angle;
                        let mode = config.mode;
                        let dist_cutoff = config.grid_width;
                        let multiple_bin = &config.multiple_bin;
                        let (pdb_query_map, query_indices, aa_dist_map ) = if verbose { 
                            measure_time!(make_query_map(
                                &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, multiple_bin,
                                &dist_thresholds, &angle_thresholds, &aa_substitutions, dist_cutoff, serial_query,
                            ))
                        } else {
                            make_query_map(
                                &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, multiple_bin,
                                &dist_thresholds, &angle_thresholds, &aa_substitutions, dist_cutoff, serial_query,
                            )
                        };
                        let pdb_query = pdb_query_map.keys().cloned().collect::<Vec<_>>();
                        // Make filters out of filtering parameters
                        let structure_filter = StructureFilter::new(
                            total_match_count, covered_node_count, covered_node_ratio, covered_edge_count, covered_edge_ratio,
                            idf_score_cutoff, num_res_cutoff, plddt_cutoff, 
                            max_matching_node_count, max_matching_node_ratio, rmsd_cutoff,
                            _residue_count, _residue_count * (_residue_count - 1)
                        );

                        match mode {
                            IndexMode::Id => {
                                let (value_mmap, value_vec) = if verbose { measure_time!(read_u16_vector(&value_path).expect(
                                    &log_msg(FAIL, &format!("Failed to load value vector: {}", &value_path))
                                )) } else {
                                    read_u16_vector(&value_path).expect(&log_msg(FAIL, &format!("Failed to load value vector: {}", &value_path)))
                                };
                                let query_count_map = if verbose { measure_time!(count_query_idmode(
                                    &pdb_query, &pdb_query_map, &offset_table, value_vec, &lookup, 
                                    sampling_ratio, sampling_count, freq_filter, length_penalty
                                ))} else {
                                    count_query_idmode(
                                        &pdb_query, &pdb_query_map, &offset_table, value_vec, &lookup,
                                        sampling_ratio, sampling_count, freq_filter, length_penalty
                                    )
                                };

                                let mut query_count_vec: Vec<(usize, StructureResult)> = query_count_map.into_par_iter().filter(|(_k, v)| {
                                    structure_filter.filter_before_matching(v)
                                }).collect();
                                if verbose {
                                    print_log_msg(INFO, &format!("Found {} structures from inverted index", query_count_vec.len()));
                                }
                                // 
                                if verbose {
                                    measure_time!(query_count_vec.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap()));
                                } else {
                                    query_count_vec.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap());
                                }
                                // Apply top N filter if top_n is not usize::MAX
                                if top_n != usize::MAX {
                                    if verbose {
                                        print_log_msg(INFO, &format!("Limiting result to top {} structures", top_n));
                                    }
                                    query_count_vec.truncate(top_n);
                                }
                                // IF retrieve is true, retrieve matching residues
                                if !skip_match {
                                    if verbose {
                                        measure_time!(query_count_vec.par_iter_mut().for_each(|(_, v)| {
                                            #[cfg(not(feature = "foldcomp"))]
                                            let retrieval_result = retrieval_wrapper(
                                                &v.id, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                &pdb_query_map, &query_structure, &query_indices,
                                                &aa_dist_map, ca_dist_threshold,
                                            );
                                            #[cfg(feature = "foldcomp")]
                                            let retrieval_result = if using_foldcomp {
                                                retrieval_wrapper_for_foldcompdb(
                                                    v.db_key, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                    hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                    &pdb_query_map, &query_structure, &query_indices,
                                                    &aa_dist_map, ca_dist_threshold, &foldcomp_db_reader
                                                )
                                            } else {
                                                retrieval_wrapper(
                                                    &v.id, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                    hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                    &pdb_query_map, &query_structure, &query_indices,
                                                    &aa_dist_map, ca_dist_threshold,
                                                )
                                            };
                                            v.matching_residues = retrieval_result.0;
                                            v.matching_residues_processed = retrieval_result.1;
                                            v.max_matching_node_count = retrieval_result.2;
                                            v.min_rmsd_with_max_match = retrieval_result.3;
                                        }));
                                    } else {
                                        query_count_vec.par_iter_mut().for_each(|(_, v)| {
                                            #[cfg(not(feature = "foldcomp"))]
                                            let retrieval_result = retrieval_wrapper(
                                                &v.id, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                &pdb_query_map, &query_structure, &query_indices,
                                                &aa_dist_map, ca_dist_threshold,
                                            );
                                            #[cfg(feature = "foldcomp")]
                                            let retrieval_result = if using_foldcomp {
                                                retrieval_wrapper_for_foldcompdb(
                                                    v.db_key, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                    hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                    &pdb_query_map, &query_structure, &query_indices,
                                                    &aa_dist_map, ca_dist_threshold, &foldcomp_db_reader
                                                )
                                            } else {
                                                retrieval_wrapper(
                                                    &v.id, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                    hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                    &pdb_query_map, &query_structure, &query_indices,
                                                    &aa_dist_map, ca_dist_threshold,
                                                )
                                            };
                                            v.matching_residues = retrieval_result.0;
                                            v.matching_residues_processed = retrieval_result.1;
                                            v.max_matching_node_count = retrieval_result.2;
                                            v.min_rmsd_with_max_match = retrieval_result.3;
                                        });
                                    }
                                    
                                    // Filter query_count_vec with reasonable retrieval results
                                    query_count_vec.retain(|(_, v)| structure_filter.filter_after_matching(v));
                                    drop(value_mmap);
                                    return query_count_vec;
                                }
                                drop(value_mmap);
                                query_count_vec
                            },
                            IndexMode::Big => {

                                let query_count_map = if verbose { measure_time!(count_query_bigmode(
                                    &pdb_query, &pdb_query_map, &big_index, &lookup, 
                                    sampling_ratio, sampling_count, freq_filter, length_penalty
                                )) } else {
                                    count_query_bigmode(
                                        &pdb_query, &pdb_query_map, &big_index, &lookup,
                                        sampling_ratio, sampling_count, freq_filter, length_penalty
                                    )
                                };
                                let mut query_count_vec: Vec<(usize, StructureResult)> = query_count_map.into_par_iter().filter(|(_k, v)| {
                                    structure_filter.filter_before_matching(v)
                                }).collect();

                                if verbose {
                                    print_log_msg(INFO, &format!("Found {} structures from inverted index", query_count_vec.len()));
                                }

                                // 
                                if verbose {
                                    measure_time!(query_count_vec.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap()));
                                } else {
                                    query_count_vec.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap());
                                }
                                // Apply top N filter if top_n is not usize::MAX
                                if top_n != usize::MAX {
                                    if verbose {
                                        print_log_msg(INFO, &format!("Limiting result to top {} structures", top_n));
                                    }
                                    query_count_vec.truncate(top_n);
                                }
                                
                                // IF retrieve is true, retrieve matching residues
                                if !skip_match {
                                    if verbose {
                                        measure_time!(query_count_vec.par_iter_mut().for_each(|(_, v)| {
                                            #[cfg(not(feature = "foldcomp"))]
                                            let retrieval_result = retrieval_wrapper(
                                                &v.id, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                &pdb_query_map, &query_structure, &query_indices,
                                                &aa_dist_map, ca_dist_threshold,
                                            );
                                            #[cfg(feature = "foldcomp")]
                                            let retrieval_result = if using_foldcomp {
                                                retrieval_wrapper_for_foldcompdb(
                                                    v.db_key, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                    hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                    &pdb_query_map, &query_structure, &query_indices,
                                                    &aa_dist_map, ca_dist_threshold, &foldcomp_db_reader
                                                )
                                            } else {
                                                retrieval_wrapper(
                                                    &v.id, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                    hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                    &pdb_query_map, &query_structure, &query_indices,
                                                    &aa_dist_map, ca_dist_threshold,
                                                )
                                            };
                                            v.matching_residues = retrieval_result.0;
                                            v.matching_residues_processed = retrieval_result.1;
                                            v.max_matching_node_count = retrieval_result.2;
                                            v.min_rmsd_with_max_match = retrieval_result.3;
                                        }));
                                    } else {
                                        query_count_vec.par_iter_mut().for_each(|(_, v)| {
                                            #[cfg(not(feature = "foldcomp"))]
                                            let retrieval_result = retrieval_wrapper(
                                                &v.id, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                &pdb_query_map, &query_structure, &query_indices,
                                                &aa_dist_map, ca_dist_threshold,
                                            );
                                            #[cfg(feature = "foldcomp")]
                                            let retrieval_result = if using_foldcomp {
                                                retrieval_wrapper_for_foldcompdb(
                                                    v.db_key, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                    hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                    &pdb_query_map, &query_structure, &query_indices,
                                                    &aa_dist_map, ca_dist_threshold, &foldcomp_db_reader
                                                )
                                            } else {
                                                retrieval_wrapper(
                                                    &v.id, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                                    hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                                    &pdb_query_map, &query_structure, &query_indices,
                                                    &aa_dist_map, ca_dist_threshold,
                                                )
                                            };
                                            v.matching_residues = retrieval_result.0;
                                            v.matching_residues_processed = retrieval_result.1;
                                            v.max_matching_node_count = retrieval_result.2;
                                            v.min_rmsd_with_max_match = retrieval_result.3;
                                        });
                                    }

                                    // Filter query_count_vec with reasonable retrieval results
                                    query_count_vec.retain(|(_, v)| v.matching_residues.len() > 0);
                                    return query_count_vec;
                                }
                                query_count_vec
                            },
                        } // match mode
                    }
                ).reduce(|| Vec::new(), |mut acc, x| {
                    acc.extend(x);
                    acc
                });
                drop(query_residues);
                let match_filter= MatchFilter::new(
                    connected_node_count, connected_node_ratio, idf_score_cutoff,
                    rmsd_cutoff, _residue_count,
                );

                match query_mode {
                    QueryMode::PerMatchDefault | QueryMode::PerMatchSortByScore => {
                        let mut match_results = convert_structure_query_result_to_match_query_results(
                            &queried_from_indices, skip_ca_match
                        );
                        match_results.retain(|(_, v)| match_filter.filter(v));
                        sort_and_print_match_query_result(
                            &mut match_results, top_n, 
                            &output_path, &query_string, output_with_superpose, header, verbose
                        );
                    }
                    QueryMode::Web => {
                        let mut match_results = convert_structure_query_result_to_match_query_results(
                            &queried_from_indices, skip_ca_match
                        );
                        match_results.retain(|(_, v)| match_filter.filter(v));
                        // If web, set superpose to true.
                        sort_and_print_match_query_result(
                            &mut match_results, MAX_NUM_LINES_FOR_WEB,
                            &output_path, &query_string, true, header, verbose
                        );
                    }
                    QueryMode::PerStructureSortByRmsd => {
                        sort_and_print_structure_query_result(
                            &mut queried_from_indices,  true, &output_path, 
                            &query_string, header, verbose
                        );
                    }
                    QueryMode::PerStructureSortByScore | QueryMode::SkipMatch => {
                        sort_and_print_structure_query_result(
                            &mut queried_from_indices, false, &output_path, 
                            &query_string, header, verbose
                        );
                    }
                    _ => {}
                }
                drop(queried_from_indices);
                drop(query_structure);
            }); // queries
            drop(loaded_index_vec);
            drop(big_offset_mmap);
            drop(big_index);
        }, // AppArgs::Query
        _ => {
            eprintln!("{}", HELP_QUERY);
            std::process::exit(1);
        }
    }
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
        let env = AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            skip_match: false,
            dist_threshold: String::from("0.5"),
            angle_threshold: String::from("5.0"),
            ca_dist_threshold: 1.0,
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.0,
            covered_edge_count: 0,
            covered_edge_ratio: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            idf_score_cutoff: 0.0,
            connected_node_count: 0,
            connected_node_ratio: 0.0,
            num_res_cutoff: 3000,
            plddt_cutoff: 0.0,
            rmsd_cutoff: 1.0,
            top_n: 1000,
            web_mode: false,
            sampling_count: None,
            sampling_ratio: None,
            freq_filter: None,
            length_penalty: None,
            sort_by_rmsd: true,
            sort_by_score: false,
            output_per_structure: false,
            output_per_match: true,
            output_with_superpose: false,
            skip_ca_match: false,
            header: true,
            serial_query: false,
            output: String::from(""),
            verbose: true,
            help: false,
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
            let env = AppArgs::Query {
                pdb_path,
                query_string,
                threads,
                index_path,
                skip_match: false,
                dist_threshold: String::from("0.5"),
                angle_threshold: String::from("5.0"),
                ca_dist_threshold: 1.0,
                total_match_count: 0,
                covered_node_count: 0,
                covered_node_ratio: 0.0,
                covered_edge_count: 0,
                covered_edge_ratio: 0.0,
                idf_score_cutoff: 0.0,
                connected_node_count: 0,
                connected_node_ratio: 0.0,
                max_matching_node_count: 0,
                max_matching_node_ratio: 0.0,
                num_res_cutoff: 3000,
                plddt_cutoff: 0.0,
                rmsd_cutoff: 1.0,
                top_n: 1000,
                web_mode: false,
                sampling_count: None,
                sampling_ratio: None,
                freq_filter: None,
                length_penalty: None,
                sort_by_rmsd: false,
                sort_by_score: true,
                output_per_structure: true,
                output_per_match: false,
                output_with_superpose: true,
                skip_ca_match: false,
                header: true,
                serial_query: false,
                output: String::from(""),
                verbose: true,
                help: false,
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
        let index_path = Some(String::from("analysis/e_coli/test"));
        let env = AppArgs::Query {
            pdb_path,
            query_string,
            threads, 
            index_path,
            skip_match: true,
            dist_threshold: String::from("0.5"),
            angle_threshold: String::from("5.0"),
            ca_dist_threshold: 1.0,
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.0,
            covered_edge_count: 0,
            covered_edge_ratio: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            idf_score_cutoff: 0.0,
            connected_node_count: 0,
            connected_node_ratio: 0.0,
            num_res_cutoff: 3000,
            plddt_cutoff: 0.0,
            rmsd_cutoff: 1.0,
            top_n: 1000,
            web_mode: false,
            sampling_count: None,
            sampling_ratio: None,
            freq_filter: None,
            length_penalty: None,
            sort_by_rmsd: false,
            sort_by_score: true,
            output_per_structure: true,
            output_per_match: false,
            output_with_superpose: true,
            skip_ca_match: false,
            header: true,
            serial_query: false,
            output: String::from(""),
            verbose: true,
            help: false,
        };
        query_pdb(env);
    }
}