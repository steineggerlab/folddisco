// File: main.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

//! Main entry point for FoldDisco CLI

// use crate::*;
use folddisco::cli::{workflows::{build_index, benchmark, query_pdb}, *};
use git_version::git_version;

const VERSION_STRING: &str = git_version!(
    args = ["--abbrev=40", "--always"], fallback = "unknown"
);
const HELP: &str = "\
usage: folddisco <command> [<args>]

subcommands:
  index     Create a new index table from multiple protein structures
  query     Query a motif from an index table
  benchmark Benchmark the performance of folddisco
  version   Print version information

options:
  -h, --help                 Print this help menu
";

fn parse_arg() -> Result<AppArgs, Box<dyn std::error::Error>> {
    let mut args = pico_args::Arguments::from_env();
    match args
        .subcommand()
        .expect("Failed to parse subcommand")
        .as_deref()
    {
        Some("index") => Ok(AppArgs::Index {
            pdb_container: args.opt_value_from_str(["-p", "--pdbs"])?,
            hash_type: args.value_from_str(["-y", "--type"]).unwrap_or("default".into()),
            index_path: args.value_from_str(["-i", "--index"]).unwrap_or("".into()),
            num_threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
            mode: args.value_from_str(["-m", "--mode"]).unwrap_or("id".into()),
            num_bin_dist: args.value_from_str(["-d", "--distance"]).unwrap_or(0),
            num_bin_angle: args.value_from_str(["-a", "--angle"]).unwrap_or(0),
            multiple_bins: args.opt_value_from_str("--multiple-bins")?,
            grid_width: args.value_from_str(["-g", "--grid"]).unwrap_or(20.0),
            chunk_size: args.value_from_str(["-c", "--chunk"]).unwrap_or(65536),
            max_residue: args.value_from_str(["-n", "--residue"]).unwrap_or(50000),
            recursive: args.contains(["-r", "--recursive"]),
            mmap_on_disk: args.contains("--mmap-on-disk"),
            id_type: args.value_from_str("--id").unwrap_or("relpath".into()),
            verbose: args.contains(["-v", "--verbose"]),
            help: args.contains(["-h", "--help"]),
        }),
        Some("query") => Ok(AppArgs::Query {
            pdb_path: args.value_from_str(["-p", "--pdb"]).unwrap_or("".into()),
            query_string: args.value_from_str(["-q", "--query"]).unwrap_or("".into()),
            threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
            index_path: args.opt_value_from_str(["-i", "--index"])?,
            skip_match: args.contains("--skip-match"),
            // Filtering parameters
            dist_threshold: args.value_from_str(["-d", "--distance"]).unwrap_or("0.5".into()),
            angle_threshold: args.value_from_str(["-a", "--angle"]).unwrap_or("5".into()),
            ca_dist_threshold: args.value_from_str("--ca-distance").unwrap_or(1.0),
            total_match_count: args.value_from_str("--total-match").unwrap_or(0),
            covered_node_count: args.value_from_str("--covered-node").unwrap_or(0),
            covered_node_ratio: args.value_from_str("--covered-node-ratio").unwrap_or(0.0),
            max_matching_node_count: args.value_from_str("--max-node").unwrap_or(0),
            max_matching_node_ratio: args.value_from_str("--max-node-ratio").unwrap_or(0.0),
            idf_score_cutoff: args.value_from_str("--score").unwrap_or(0.0),
            connected_node_count: args.value_from_str("--connected-node").unwrap_or(0),
            connected_node_ratio: args.value_from_str("--connected-node-ratio").unwrap_or(0.0),
            num_res_cutoff: args.value_from_str("--num-residue").unwrap_or(50000),
            plddt_cutoff: args.value_from_str("--plddt").unwrap_or(0.0),
            rmsd_cutoff: args.value_from_str("--rmsd").unwrap_or(0.0),
            top_n: args.value_from_str("--top").unwrap_or(usize::MAX),
            web_mode: args.contains("--web"), // Web mode for output
            // Query filtering
            sampling_count: args.opt_value_from_str("--sampling-count")?,
            sampling_ratio: args.opt_value_from_str("--sampling-ratio")?,
            freq_filter: args.opt_value_from_str("--freq-filter")?,
            length_penalty: args.opt_value_from_str("--length-penalty")?,
            // Sorting strategy (comma-separated keys)
            sort_by: args.value_from_str("--sort-by").unwrap_or("node_count,rmsd".into()),
            // Output mode
            output_per_structure: args.contains("--per-structure"),
            output_per_match: args.contains("--per-match"),
            output_with_superpose: args.contains("--superpose"), // Print target CA, U, T
            skip_ca_match: args.contains("--skip-ca-match"),
            partial_fit: args.contains("--partial-fit"), // Enable LMS based superposition.
            header: args.contains("--header"),
            serial_query: args.contains("--serial-index"),
            output: args.value_from_str(["-o", "--output"]).unwrap_or("".into()),
            verbose: args.contains(["-v", "--verbose"]),
            help: args.contains(["-h", "--help"]),
        }),
        Some("benchmark") => Ok(AppArgs::Benchmark {
            result: args.opt_value_from_str(["-r", "--result"])?,
            answer: args.opt_value_from_str(["-a", "--answer"])?,
            neutral: args.opt_value_from_str(["-n", "--neutral"])?,
            index: args.opt_value_from_str(["-i", "--index"])?,
            input: args.opt_value_from_str("--input")?,
            format: args.value_from_str(["-f", "--format"]).unwrap_or("tsv".into()),
            fp: args.opt_value_from_str("--fp")?,
            threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
            afdb_to_uniprot: args.contains("--afdb-to-uniprot"),
            column_result: args.value_from_str("--column-result").unwrap_or(0),
            column_answer: args.value_from_str("--column-answer").unwrap_or(0),
            column_neutral: args.value_from_str("--column-neutral").unwrap_or(0),
            header_result: args.contains("--header-result"),
            header_answer: args.contains("--header-answer"),
            header_neutral: args.contains("--header-neutral"),
        }),
        Some("test") => Ok(AppArgs::Test {
            index_path: args.value_from_str(["-i", "--index"])?,
            verbose: args.contains(["-v", "--verbose"]),
        }),
        Some("version") => {
            println!("{}", VERSION_STRING);
            std::process::exit(0);
        },
        Some(_) => Err("Invalid subcommand".into()),
        None => Ok(AppArgs::Global {
            help: args.contains(["-h", "--help"]),
        }),
    }
}

fn main() {


    let parsed_args = parse_arg().unwrap_or_else(|e| {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    });
    match parsed_args {
        AppArgs::Global { help: _ } => {
            print_logo();
            eprintln!("{}", HELP);
        }
        AppArgs::Index { help, .. } => {
            if help {
                print_logo();
                eprintln!("{}", workflows::build_index::HELP_INDEX);
            } else {
                build_index::build_index(parsed_args);
            }
        }
        AppArgs::Query { help, .. } => {
            if help {
                print_logo();
                eprintln!("{}", workflows::query_pdb::HELP_QUERY);
            } else {
                query_pdb::query_pdb(parsed_args);
            }
        }
        AppArgs::Benchmark { .. } => {
            benchmark::benchmark(parsed_args);
        }
        AppArgs::Test { .. } => {
            println!("Testing");
            // temp::query_test_for_swissprot(parsed_args);
        }
    }
}
