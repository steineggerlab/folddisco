// File: main.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

//! Main entry point for FoldDisco CLI

// use crate::*;
use folddisco::cli::{workflows::{build_index, benchmark, query_pdb}, *};
const HELP: &str = "\
USAGE: folddisco index [OPTIONS] -p <PDBS...> -i <INDEX> -y <TYPE>
       folddisco query [OPTIONS] -p <PDB> -i <INDEX> -q <QUERY>
       folddisco benchmark [OPTIONS] -r <RESULT> -a <ANSWER> -i <INDEX>

SUBCOMMANDS:
  index     Create a new index table from multiple protein structures
  query     Query a motif from an index table
  benchmark Benchmark the performance of FoldDisco
OPTIONS:
  -t, --threads <THREADS>    Number of threads to use
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
            index_path: args.value_from_str(["-i", "--index"]).unwrap_or("folddisco_index".into()),
            num_threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
            mode: args.value_from_str(["-m", "--mode"]).unwrap_or("id".into()),
            num_bin_dist: args.value_from_str(["-d", "--distance"]).unwrap_or(0),
            num_bin_angle: args.value_from_str(["-a", "--angle"]).unwrap_or(0),
            grid_width: args.value_from_str(["-g", "--grid"]).unwrap_or(20.0),
            chunk_size: args.value_from_str(["-c", "--chunk"]).unwrap_or(65535), // TODO: 8192
            max_residue: args.value_from_str(["-n", "--residue"]).unwrap_or(50000),
            recursive: args.contains(["-r", "--recursive"]),
            id_type: args.value_from_str("--id").unwrap_or("relpath".into()),
            verbose: args.contains(["-v", "--verbose"]),
            help: args.contains(["-h", "--help"]),
        }),
        Some("query") => Ok(AppArgs::Query {
            pdb_path: args.value_from_str(["-p", "--pdb"]).unwrap_or("".into()),
            query_string: args.value_from_str(["-q", "--query"]).unwrap_or("".into()),
            threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
            index_path: args.opt_value_from_str(["-i", "--index"])?,
            retrieve: args.contains(["-r", "--retrieve"]),
            amino_acid: args.value_from_str("--amino-acid").unwrap_or(0),
            dist_threshold: args.opt_value_from_str(["-d", "--distance"])?,
            angle_threshold: args.opt_value_from_str(["-a", "--angle"])?,
            match_cutoff: args.opt_value_from_str(["-m", "--match"])?, 
            score_cutoff: args.value_from_str(["-s", "--score"]).unwrap_or(0.0),
            num_res_cutoff: args.value_from_str(["-n", "--residue"]).unwrap_or(50000),
            plddt_cutoff: args.value_from_str(["-l", "--plddt"]).unwrap_or(0.0),
            node_count: args.value_from_str("--node").unwrap_or(2),
            header: args.contains("--header"),
            id_type: args.value_from_str("--id").unwrap_or("relpath".into()),
            verbose: args.contains(["-v", "--verbose"]),
            help: args.contains(["-h", "--help"]),
        }),
        Some("benchmark") => Ok(AppArgs::Benchmark {
            result: args.opt_value_from_str(["-r", "--result"])?,
            answer: args.opt_value_from_str(["-a", "--answer"])?,
            index: args.opt_value_from_str(["-i", "--index"])?,
            format: args.value_from_str(["-f", "--format"]).unwrap_or("tsv".into()),
        }),
        Some("test") => Ok(AppArgs::Test {
            index_path: args.value_from_str(["-i", "--index"])?,
            verbose: args.contains(["-v", "--verbose"]),
        }),
        Some(_) => Err("Invalid subcommand".into()),
        None => Ok(AppArgs::Global {
            help: args.contains(["-h", "--help"]),
        }),
    }
}

fn main() {
    // Init
    print_logo();

    let parsed_args = parse_arg().unwrap_or_else(|e| {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    });
    match parsed_args {
        AppArgs::Global { help: _ } => {
            eprintln!("{}", HELP);
        }
        AppArgs::Index { help, .. } => {
            if help {
                eprintln!("{}", workflows::build_index::HELP_INDEX);
            } else {
                build_index::build_index(parsed_args);
            }
        }
        AppArgs::Query { help, .. } => {
            if help {
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
