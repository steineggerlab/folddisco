// File: main.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

// use crate::*;
use motifsearch::cli::{workflows::{build_index, temp, query_pdb}, *};
use motifsearch::prelude::*;
const HELP: &str = "\
USAGE: motifsearch index [OPTIONS] <PDBS...>
       motifsearch query [OPTIONS] <PDB> <INDEX>

SUBCOMMANDS:
  index     Create a new index table from multiple protein structures
  query     Query a motif from an index table
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
            pdb_dir: args.opt_value_from_str(["-d", "--pdb-dir"])?,
            pdb_path_vec: Vec::new(),
            index_path: args.value_from_str(["-i", "--index-path"])?,
            num_threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
            verbose: args.contains(["-v", "--verbose"]),
            help: args.contains(["-h", "--help"]),
        }),
        Some("query") => Ok(AppArgs::Query {
            threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
            index_path: args.opt_value_from_str(["-i", "--index-path"])?,
            help: args.contains(["-h", "--help"]),
        }),
        Some("test") => Ok(AppArgs::Test {
            index_path: args.value_from_str(["-i", "--index-path"])?,
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
    print_fold_disco_logo();

    let parsed_args = parse_arg().unwrap_or_else(|e| {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    });
    match parsed_args {
        AppArgs::Global { help } => {
            if !help { eprintln!("No subcommand specified."); }
            eprintln!("{}", HELP);
        }
        AppArgs::Index { help, .. } => {
            if help {
                eprintln!("{}", HELP);
            } else {
                build_index::build_index(parsed_args);
            }
        }
        AppArgs::Query { help, threads, .. } => {
            if help {
                eprintln!("{}", HELP);
            } else {
                eprintln!("{} Querying with {} threads...",  INFO, threads);
                query_pdb::query_pdb(parsed_args);
                eprintln!("Querying done.");
            }
        }
        AppArgs::Test { .. } => {
            eprintln!("Testing");
            temp::query_test_for_swissprot(parsed_args);
        }
    }
}
