// File: analyze.rs
// Created: 2025-12-26 20:33:48
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2025 Hyunbin Kim, All rights reserved

//! This workflow handles distribution analysis of given index.
//! If only index is given, it summarizes the overall distribution of the hash values in the index.
//! If PDB files or query motifs are also given, it analyzes the distribution of matches for


use crate::cli::config::read_index_config_from_file;
use crate::cli::{AppArgs, print_logo};
use crate::controller::io::get_lookup_and_type;
use crate::controller::summary::{count_encodings, save_summary};
use crate::prelude::*;

pub const HELP_ANALYZE: &str = "\
usage: folddisco analyze -i <index> [options]

input:
    -i, --index <INDEX_PREFIX>   Path to index prefix
    -p, --pdbs <PDB_PATH>        Directory of Foldcomp DB containing PDB files (optional)

output:
    -o, --output <OUTPUT_PATH>   Path to save analysis results (optional)

options:
    --top <INT>                  Limit output to top N encodings by frequency [default: 10]

general options:
    -t, --threads <INT>          Number of threads to use [default: 1]
    -v, --verbose                Print verbose messages
    -h, --help                   Print this help menu
";

pub fn analyze(env: AppArgs) {
    match env {
        AppArgs::Analyze {
            index_path,
            pdb_container,
            output,
            top_n,
            threads,
            verbose,
            help: _,
        } => {
            if verbose { print_logo(); }

            // index_path is required
            if index_path.is_none() {
                eprintln!("{}", HELP_ANALYZE);
                std::process::exit(1);
            }
 
            // Check if pdb_container is provided
            // If provided, compare the hash distributions of the PDBs against the index
            // If not provided, just summarize the index distribution
            let compare_with_pdbs = pdb_container.is_some();
            let index_path = index_path.unwrap();
            let output_prefix = match output {
                Some(p) => p,
                None => {
                    if compare_with_pdbs {
                        format!("{}_vs_{}",
                            index_path,
                            pdb_container.as_ref().unwrap().replace("/", "_")
                        )
                    } else {
                        format!("{}_summary", index_path)
                    }
                }
            };
            
            if verbose {
                if compare_with_pdbs {
                    print_log_msg(INFO, &format!(
                        "Comparing encoding distribution of {:?} against {}",
                        pdb_container.as_ref().unwrap(), &index_path
                    ));
                } else {
                    print_log_msg(INFO, &format!(
                        "Summarizing encoding distribution of index {:?}", &index_path
                    ));
                }
            }
            // Set thread pool
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().expect(
                &log_msg(FAIL, "Failed to build thread pool")
            );
            // Load index
            let (index, _offset_mmap) = measure_time!(
                load_folddisco_index(&index_path), verbose
            );
            // Load lookup and config
            let (_, hash_type_path) = get_lookup_and_type(&index_path);
            let config = read_index_config_from_file(&hash_type_path);
            let hash_type = config.hash_type;
            let nbin_dist = config.num_bin_dist;
            let nbin_angle = config.num_bin_angle;
            
            // Execute workflow
            if compare_with_pdbs {
                // let analysis_result = 
                todo!("analyze_enrichment(&index, &offset_mmap, pdb_container.unwrap(), verbose)");
                // todo!("save_enrichment_result(&analysis_result, &output_prefix, verbose)");
            } else {
                let summary = count_encodings(
                    &index, hash_type, nbin_dist, nbin_angle, verbose
                );
                let saved = save_summary(
                    &summary, &output_prefix, top_n, verbose
                );
                if saved.is_err() {
                    print_log_msg(FAIL, &format!(
                        "Failed to save summary result to {}", &output_prefix
                    ));
                }
            }

        }
        _ => {
            eprintln!("{}", HELP_ANALYZE);
            std::process::exit(1);
        }
    }
}