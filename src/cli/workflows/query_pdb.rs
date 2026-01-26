// File: query_pdb.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for querying PDB files
// When querying PDB files, we need index table and query file.

use std::io::BufRead;

use rayon::prelude::*;

use crate::cli::config::read_index_config_from_file;
use crate::controller::filter::{MatchFilter, StructureFilter};
use crate::controller::mode::QueryMode;
use crate::controller::sort::{MatchSortStrategy, StructureSortStrategy};
use crate::cli::*;
use crate::controller::io::{
    read_compact_structure, 
    check_and_get_indices, get_lookup_and_type,
};
use crate::controller::query::{
    make_query_map, parse_threshold_string
};
use crate::controller::count_query::count_query;
use crate::controller::result::{
    convert_structure_query_result_to_match_query_results, 
    sort_and_print_match_query_result, sort_and_print_structure_query_result, StructureResult
};
use crate::controller::retrieve::retrieval_wrapper;
use crate::index::indextable::load_folddisco_index;
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::*;

#[cfg(feature = "foldcomp")]
use crate::controller::retrieve::retrieval_wrapper_for_foldcompdb;
#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;
#[cfg(feature = "foldcomp")]
use crate::structure::io::StructureFileFormat;
#[cfg(feature = "foldcomp")]
use crate::controller::io::get_foldcomp_db_path_with_prefix;

pub const HELP_QUERY: &str = "\
usage: folddisco query -p <i:PDB> -q <QUERY> -i <i:INDEX> [OPTIONS] 

input/output:
 -p, --pdb <PATH>                 Path of PDB file to query
 -q, --query <STR>                Query string that specifies residues or a text file containing query
 -i, --index <PATH>               Path of index table to load [REQUIRED]
 -o, --output <PATH>              Output file path [stdout]
 
search parameters:
 -t, --threads <INT>              Number of threads [1]
 -d, --distance <FLOAT>           Distance threshold in Angstroms. Multiple values can be separated by comma [0.5]
 -a, --angle <FLOAT>              Angle threshold. Multiple values can be separated by comma [5.0]
 --ca-distance <FLOAT>            C-alpha distance threshold in matching residues [1.0]
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
 --max-node <INT>                 Filter out structures of maximum matching node size smaller than given value [0]
 --max-node-ratio <FLOAT>         Filter out structures of maximum matching node size smaller than given ratio [0.0]
 --score <FLOAT>                  IDF score cutoff [0.0]
 --connected-node <INT>           Filter out structures/matches with connected node count smaller than given value [0]
 --connected-node-ratio <FLOAT>   Filter out structures/matches with connected node count smaller than given ratio [0.0]
 --num-residue <INT>              Number of residues cutoff [50000]
 --plddt <FLOAT>                  pLDDT cutoff [0.0]
 --rmsd <FLOAT>                   Maximum RMSD cutoff [no limit]
 --tm-score <FLOAT>               Minimum TM-score cutoff [0.0]
 --gdt-ts <FLOAT>                 Minimum GDT-TS cutoff. Thresholds: 1.0Å, 2.0Å, 4.0Å, 8.0Å [0.0]
 --gdt-ha <FLOAT>                 Minimum GDT-HA cutoff. Thresholds: 0.5Å, 1.0Å, 2.0Å, 4.0Å [0.0]
 --chamfer <FLOAT>                Maximum Chamfer distance cutoff. Chamfer distance is mean of nearest neighbor distances between two point clouds [no limit]
 --hausdorff <FLOAT>              Maximum Hausdorff distance cutoff. Hausdorff distance is maximum of nearest neighbor distances between two point clouds [no limit]
 --top <INT>                      Limit output to top N structures based on IDF score [all]

display options:
 --header                         Print header in output
 --web                            Print output for web
 --per-structure                  Print output per structure
 --per-match                      Print output per match. Not working with --skip-match
 --format-output <KEYS>           Comma-separated column names to output
                                  - Per-match: qid, tid, nid, db_key, node_count, idf, rmsd, matching_residues, u_matrix, t_vector,
                                    matching_coordinates, query_residues, tm_score, gdt_ts, gdt_ha, chamfer_distance, hausdorff_distance
                                  - Per-structure: qid, tid, nid, db_key, total_match_count, node_count, edge_count, idf, nres, plddt,
                                    max_node_cov, min_rmsd, matching_residues, query_residues
                                  - Example: --format-output tid,idf,rmsd,tm_score
 --sort-by <KEYS>                 Comma-separated sort keys with optional :asc or :desc [default: node_count:desc,rmsd:asc]
                                  - Per-match: node_count, idf, rmsd, tm_score, gdt_ts, gdt_ha, chamfer_distance, hausdorff_distance
                                  - Per-structure: max_node_count, node_count, idf, min_rmsd, total_match_count, edge_count, nres, plddt
                                  - Example: --sort-by tm_score,rmsd or --sort-by idf:desc
 --skip-ca-match                  Print matching residues before C-alpha distance check
 --partial-fit                    Superposition will find the best aligning substructure using LMS (Least Median of Squares)
 --superpose                      Print U, T, CA of matching residues

general options:
 -v, --verbose                    Print verbose messages
 -h, --help                       Print this help menu

examples:
# Search with default settings. This will print out matching motifs sorted by node count then RMSD.
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6

# Print custom columns (tid, idf, RMSD, and TM-score only)
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --format-output tid,idf,rmsd,tm_score

# Print matches sorted by node count and TM-score
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --sort-by node_count,tm_score

# Print per-structure results sorted by max node count and IDF score
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --per-structure --sort-by max_node_count,idf

# Print per-structure results sorted by IDF only
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --per-structure --sort-by idf

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
            max_matching_node_count,
            max_matching_node_ratio,
            idf_score_cutoff,
            connected_node_count,
            connected_node_ratio,
            num_res_cutoff,
            plddt_cutoff,
            rmsd_cutoff,
            tm_score_cutoff,
            gdt_ts_cutoff,
            gdt_ha_cutoff,
            chamfer_distance_cutoff,
            hausdorff_distance_cutoff,
            top_n,
            web_mode,
            sampling_count,
            sampling_ratio,
            freq_filter,
            length_penalty,
            sort_by,
            format_output,
            output_per_structure,
            output_per_match,
            output_with_superpose,
            skip_ca_match,
            partial_fit,
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
            
            // Determine query mode first to decide which sorting strategy to use
            let query_mode = QueryMode::from_flags(
                skip_match, web_mode, output_per_structure, output_per_match
            );
            
            // Error handling
            if query_mode == QueryMode::ContradictoryPrintError {
                print_log_msg(FAIL, 
                    "Cannot print output per structure and per match at the same time. Use either --per-structure or --per-match"
                );
                std::process::exit(1);
            }
            
            
            // Parse the appropriate sorting strategy based on query mode
            // For per-structure output modes, use StructureSortStrategy
            // For per-match output modes, use MatchSortStrategy
            let use_structure_sort = matches!(
                query_mode,
                QueryMode::PerStructure | QueryMode::SkipMatch
            );
            
            let match_sort_strategy = if !use_structure_sort {
                // If sort_by is given, parse it
                if sort_by.is_empty() {
                    MatchSortStrategy::default()
                } else {
                    MatchSortStrategy::from_str(&sort_by)
                        .unwrap_or_else(|e| {
                            print_log_msg(FAIL, &format!("Error parsing --sort-by: {}", e));
                            std::process::exit(1);
                        })
                }
            } else {
                MatchSortStrategy::default()
            };
            
            let structure_sort_strategy = if use_structure_sort {
                // If sort_by is given, parse it
                if sort_by.is_empty() {
                    StructureSortStrategy::default()
                } else {
                    StructureSortStrategy::from_str(&sort_by)
                        .unwrap_or_else(|e| {
                            print_log_msg(FAIL, &format!("Error parsing --sort-by: {}", e));
                            std::process::exit(1);
                        })
                }
            } else {
                StructureSortStrategy::default()
            };
            
            // Parse format_output if provided
            let parsed_columns: Option<Vec<String>> = format_output.map(|cols| {
                cols.split(',').map(|s| s.trim().to_string()).collect()
            });
            let column_refs: Option<Vec<&str>> = parsed_columns.as_ref().map(|cols| {
                cols.iter().map(|s| s.as_str()).collect()
            });

            // Print query mode and sorting strategy
            if verbose  {
                if use_structure_sort {
                    print_log_msg(INFO, &format!("Printing results {} sorting with {}", query_mode, structure_sort_strategy));
                } else {
                    print_log_msg(INFO, &format!("Printing results {} sorting with {}", query_mode, match_sort_strategy));
                }
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
            
            // Set thread pool if there's no global thread pool yet
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            
            // Get index paths
            let index_paths = check_and_get_indices(index_path.clone(), verbose);
            if verbose {
                print_log_msg(INFO, &format!("Found {} index file(s). Querying with {} threads", index_paths.len(), threads));
            }
            
            // Always use big index
            let index_prefix = index_paths[0].clone();
            let (index, offset_mmap) = measure_time!(load_folddisco_index(&index_prefix), verbose);
            
            // Load lookup and config
            let (lookup_path, hash_type_path) = get_lookup_and_type(&index_prefix);
            let config = read_index_config_from_file(&hash_type_path);
            let lookup = measure_time!(load_lookup_from_file(&lookup_path), verbose);

            
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
    
            // Load foldcomp db 
            #[cfg(feature = "foldcomp")]
            let using_foldcomp = config.foldcomp_db.is_some() && config.input_format == StructureFileFormat::FCZDB;

            #[cfg(feature = "foldcomp")]
            let foldcomp_db_reader = match config.input_format {
                StructureFileFormat::FCZDB => {
                    if !skip_match {
                        let mut foldcomp_db_path = config.foldcomp_db.clone().unwrap();
                        // If foldcomp_db_path is not a valid path, check foldcomp db with index prefix
                        if !std::path::PathBuf::from(&foldcomp_db_path).is_file() {
                            let local_foldcomp_db_path = get_foldcomp_db_path_with_prefix(&index_prefix);
                            if local_foldcomp_db_path.is_some() {
                                foldcomp_db_path = local_foldcomp_db_path.unwrap();
                            }
                        }
                        measure_time!(FoldcompDbReader::new(foldcomp_db_path.as_str()), verbose)
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
                
                let residue_count = if query_residues.is_empty() {
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

                // Get query map for the index
                let hash_type = config.hash_type;
                let num_bin_dist = config.num_bin_dist;
                let num_bin_angle = config.num_bin_angle;
                let dist_cutoff = config.grid_width;
                let multiple_bin = &config.multiple_bin;
                let total_structures = lookup.len() as f32;
                        
                let (pdb_query_map, query_indices, aa_dist_map ) = measure_time!(make_query_map(
                    &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, multiple_bin,
                    &dist_thresholds, &angle_thresholds, &aa_substitutions, dist_cutoff, serial_query,
                    &Some(&index), total_structures
                ), verbose);

                let pdb_query = pdb_query_map.keys().cloned().collect::<Vec<_>>();
                // Make filters out of filtering parameters
                let structure_filter = StructureFilter::new(
                    total_match_count, covered_node_count, covered_node_ratio,
                    idf_score_cutoff, num_res_cutoff, plddt_cutoff, 
                    max_matching_node_count, max_matching_node_ratio, rmsd_cutoff,
                    residue_count,
                );

                let query_count_map = measure_time!(count_query(
                    &pdb_query, &pdb_query_map, &index, &lookup,
                    sampling_ratio, sampling_count, freq_filter, length_penalty
                ), verbose);
                let mut query_count_vec: Vec<(usize, StructureResult)> = query_count_map.into_par_iter().filter(|(_k, v)| {
                    structure_filter.filter_before_matching(v)
                }).collect();

                if verbose {
                    print_log_msg(INFO, &format!("Found {} structures from inverted index", query_count_vec.len()));
                }

                        // 
                measure_time!(query_count_vec.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap()), verbose);
                // Apply top N filter if top_n is not usize::MAX
                if top_n != usize::MAX {
                    if verbose {
                        print_log_msg(INFO, &format!("Limiting result to top {} structures", top_n));
                    }
                    query_count_vec.truncate(top_n);
                }
                        
                // IF retrieve is true, retrieve matching residues
                if !skip_match {
                    measure_time!(query_count_vec.par_iter_mut().for_each(|(_, v)| {
                        #[cfg(not(feature = "foldcomp"))]
                        let retrieval_result = retrieval_wrapper(
                            &v.tid, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                            hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                            &pdb_query_map, &query_structure, &query_indices,
                            &aa_dist_map, ca_dist_threshold, partial_fit
                        );
                        #[cfg(feature = "foldcomp")]
                        let retrieval_result = if using_foldcomp {
                            retrieval_wrapper_for_foldcompdb(
                                v.db_key, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                &pdb_query_map, &query_structure, &query_indices,
                                &aa_dist_map, ca_dist_threshold, partial_fit,
                                &foldcomp_db_reader
                            )
                        } else {
                            retrieval_wrapper(
                                &v.tid, MIN_CONNECTED_COMPONENT_SIZE, &pdb_query,
                                hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                &pdb_query_map, &query_structure, &query_indices,
                                &aa_dist_map, ca_dist_threshold, partial_fit,
                            )
                        };
                        v.matching_residues = retrieval_result.0;
                        v.matching_residues_processed = retrieval_result.1;
                        v.max_matching_node_count = retrieval_result.2;
                        v.min_rmsd_with_max_match = retrieval_result.3;
                    }), verbose);

                    // Filter query_count_vec with reasonable retrieval results
                    query_count_vec.retain(|(_, v)| structure_filter.filter_after_matching(v));
                }
                let mut queried_from_indices = query_count_vec;
                drop(query_residues);
                let evalue_cutoff = f64::MAX; // Currently not used
                let match_filter= MatchFilter::new(
                    connected_node_count, connected_node_ratio, idf_score_cutoff, evalue_cutoff, 
                    rmsd_cutoff, tm_score_cutoff, gdt_ts_cutoff, gdt_ha_cutoff, 
                    chamfer_distance_cutoff, hausdorff_distance_cutoff,
                    residue_count,
                );

                match query_mode {
                    QueryMode::PerMatch => {
                        let mut match_results = convert_structure_query_result_to_match_query_results(
                            &queried_from_indices, skip_ca_match, 
                            total_structures as usize, residue_count
                        );
                        match_results.retain(|(_, v)| match_filter.filter(v));
                        sort_and_print_match_query_result(
                            &mut match_results, top_n, 
                            &output_path, &pdb_path, &query_string, 
                            column_refs.as_deref(), output_with_superpose, header, verbose,
                            match_sort_strategy.clone(),
                        );
                    }
                    QueryMode::Web => {
                        let mut match_results = convert_structure_query_result_to_match_query_results(
                            &queried_from_indices, skip_ca_match, 
                            total_structures as usize, residue_count
                        );
                        match_results.retain(|(_, v)| match_filter.filter(v));
                        // If web, set superpose to true.
                        sort_and_print_match_query_result(
                            &mut match_results, MAX_NUM_LINES_FOR_WEB,
                            &output_path, &pdb_path, &query_string, 
                            column_refs.as_deref(), true, header, verbose,
                            match_sort_strategy.clone(),
                        );
                    }
                    QueryMode::PerStructure | QueryMode::SkipMatch => {
                        sort_and_print_structure_query_result(
                            &mut queried_from_indices, &output_path, 
                            &pdb_path, &query_string, column_refs.as_deref(), header, verbose, structure_sort_strategy.clone()
                        );
                    }
                    QueryMode::ContradictoryPrintError => {
                        // This should have been caught earlier, but handle it just in case
                        print_log_msg(FAIL, "Invalid query mode");
                        std::process::exit(1);
                    }
                }
                drop(queried_from_indices);
                drop(query_structure);
            }); // queries
            drop(lookup);
            drop(offset_mmap);
            drop(index);
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
        let pdb_path = String::from("data/serine_peptidase/4cha.pdb");
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
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            idf_score_cutoff: 0.0,
            connected_node_count: 0,
            connected_node_ratio: 0.0,
            num_res_cutoff: 3000,
            plddt_cutoff: 0.0,
            rmsd_cutoff: 1.0,
            tm_score_cutoff: 0.0,
            gdt_ts_cutoff: 0.0,
            gdt_ha_cutoff: 0.0,
            chamfer_distance_cutoff: 0.0,
            hausdorff_distance_cutoff: 0.0,
            top_n: 1000,
            web_mode: false,
            sampling_count: None,
            sampling_ratio: None,
            freq_filter: None,
            length_penalty: None,
            sort_by: String::from("node_count,rmsd"),
            format_output: None,
            output_per_structure: false,
            output_per_match: true,
            output_with_superpose: false,
            skip_ca_match: false,
            partial_fit: false,
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
                idf_score_cutoff: 0.0,
                connected_node_count: 0,
                connected_node_ratio: 0.0,
                max_matching_node_count: 0,
                max_matching_node_ratio: 0.0,
                num_res_cutoff: 3000,
                plddt_cutoff: 0.0,
                rmsd_cutoff: 1.0,
                tm_score_cutoff: 0.0,
                gdt_ts_cutoff: 0.0,
                gdt_ha_cutoff: 0.0,
                chamfer_distance_cutoff: 0.0,
                hausdorff_distance_cutoff: 0.0,
                top_n: 1000,
                web_mode: false,
                sampling_count: None,
                sampling_ratio: None,
                freq_filter: None,
                length_penalty: None,
                sort_by: String::from("node_count,rmsd"),
                format_output: None,
                output_per_structure: true,
                output_per_match: false,
                output_with_superpose: true,
                skip_ca_match: false,
                header: true,
                serial_query: false,
                output: String::from(""),
                verbose: true,
                partial_fit: false,
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
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            idf_score_cutoff: 0.0,
            connected_node_count: 0,
            connected_node_ratio: 0.0,
            num_res_cutoff: 3000,
            plddt_cutoff: 0.0,
            rmsd_cutoff: 1.0,
            tm_score_cutoff: 0.0,
            gdt_ts_cutoff: 0.0,
            gdt_ha_cutoff: 0.0,
            chamfer_distance_cutoff: 0.0,
            hausdorff_distance_cutoff: 0.0,
            top_n: 1000,
            web_mode: false,
            sampling_count: None,
            sampling_ratio: None,
            freq_filter: None,
            length_penalty: None,
            sort_by: String::from("node_count,rmsd"),
            format_output: None,
            output_per_structure: true,
            output_per_match: false,
            output_with_superpose: true,
            skip_ca_match: false,
            partial_fit: false,
            header: true,
            serial_query: false,
            output: String::from(""),
            verbose: true,
            help: false,
        };
        query_pdb(env);
    }
}
