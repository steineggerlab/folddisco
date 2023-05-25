//

use motifsearch::*;
use pico_args::Arguments;

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

#[derive(Debug)]
enum Subcommand {
    Index,
    Query,
    // Add subcommands here
}

enum AppArgs {
    Global {
        help: bool,
    },
    Index {
        threads: usize,
        help: bool,
    },
    Query {
        threads: usize,
        help: bool,
    },
}

fn parse_arg() -> Result<AppArgs, Box<dyn std::error::Error>> {
    let mut args = pico_args::Arguments::from_env();
    match args.subcommand().expect("Failed to parse subcommand").as_deref() {
        Some("index") => {
            Ok(AppArgs::Index {
                threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
                help: args.contains(["-h", "--help"]),
            })
        }
        Some("query") => {
            Ok(AppArgs::Query {
                threads: args.value_from_str(["-t", "--threads"]).unwrap_or(1),
                help: args.contains(["-h", "--help"]),
            })
        }
        Some(_) => {
            Err("Invalid subcommand".into())
        }
        None => {
            Ok(AppArgs::Global {
                help: args.contains(["-h", "--help"]),
            })
        }
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
        AppArgs::Index { threads, help } => {
            if help {
                println!("{}", HELP);
            } else {
                println!("Indexing with {} threads...", threads);
                // let mut index = Index::new();
                // index.build(threads);
                println!("Indexing done.");
            }
        }
        AppArgs::Query { threads, help } => {
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
