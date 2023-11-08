// File: main.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

// use crate::*;
use motifsearch::cli::{workflows::build_index, *};

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
        AppArgs::Global { help } => {
            if help {
                println!("{}", HELP);
            } else {
                println!("No subcommand specified. Try `motifsearch --help` for more information.");
            }
        }
        AppArgs::Index { help, .. } => {
            if help {
                println!("{}", HELP);
            } else {
                build_index::build_index(parsed_args);
            }
        }
        AppArgs::Query { help, threads, .. } => {
            if help {
                println!("{}", HELP);
            } else {
                println!("Querying with {} threads...", threads);
                // let mut query = Query::new();
                // query.build(threads);
                println!("Querying done.");
            }
        }
    }
}
