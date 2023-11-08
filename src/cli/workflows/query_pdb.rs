// File: query_pdb.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for querying PDB files
// When querying PDB files, we need index table and query file.

use crate::cli::*;
use crate::index::lookup::{load_lookup_from_file};
use rayon::prelude::*;
use crate::prelude::*;


const HELP_QUERY: &str = "\
USAGE: motifsearch query [OPTIONS]
Options:
    -i, --index-path <INDEX_PATH>   Path of index table to load
    -t, --threads <THREADS>     Number of threads to use
    -v, --verbose               Print verbose messages
    -h, --help                  Print this help menu
";

pub fn query_pdb(env: AppArgs) {
    todo!("Implement query_pdb");
}