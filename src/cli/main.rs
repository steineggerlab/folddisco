// File: main.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

// use crate::*;
use motifsearch::cli::{workflows::{build_index, benchmark, query_pdb}, *};
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
            // TODO: NOTE: Change to parse positional arguments for input PDBs
            pdb_dir: args.opt_value_from_str(["-d", "--pdbs"])?,
            pdb_path_vec: Vec::new(),
            hash_type: args.value_from_str(["-H", "--hash"]).unwrap_or("default".into()),
            index_path: args.value_from_str(["-i", "--index"]).unwrap_or("folddisco_index".into()),
            num_threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
            verbose: args.contains(["-v", "--verbose"]),
            help: args.contains(["-h", "--help"]),
        }),
        Some("query") => Ok(AppArgs::Query {
            pdb_path: args.value_from_str(["-d", "--pdb"]).unwrap_or("".into()),
            query_string: args.value_from_str(["-q", "--query"]).unwrap_or("".into()),
            threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
            index_path: args.opt_value_from_str(["-i", "--index"])?,
            help: args.contains(["-h", "--help"]),
            retrieve: args.contains(["-R", "--retrieve"]),
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
        AppArgs::Global { help } => {
            eprintln!("{}", HELP);
        }
        AppArgs::Index { help, .. } => {
            if help {
                eprintln!("{}", workflows::build_index::HELP_INDEX);
            } else {
                build_index::build_index(parsed_args);
            }
        }
        AppArgs::Query { help, threads, .. } => {
            if help {
                eprintln!("{}", workflows::query_pdb::HELP_QUERY);
            } else {
                query_pdb::query_pdb(parsed_args);
            }
        }
        AppArgs::Test { .. } => {
            println!("Testing");
            // temp::query_test_for_swissprot(parsed_args);
        }
    }
}
