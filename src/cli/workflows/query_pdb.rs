// File: query_pdb.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for querying PDB files
// When querying PDB files, we need index table and query file.

use std::collections::HashSet;

use crate::cli::*;
use crate::controller::query::make_query;
use crate::index::io::{read_offset_map, read_u64_vector, get_values_with_offset};
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
    match env {
        AppArgs::Query {
            threads,
            index_path,
            help,
        } => {
            // Check if arguments are valid
            if index_path.is_none() {
                eprintln!("{}", HELP_QUERY);
                std::process::exit(1);
            }
            if help {
                eprintln!("{}", HELP_QUERY);
                std::process::exit(0);
            } 
            
            // Get path. formatting without quotation marks
            let index_path = index_path.unwrap();
            let offset_path = format!("{}.offset", index_path.clone()); 
            let value_path = format!("{}.value", index_path.clone());
            let lookup_path = format!("{}.lookup", index_path.clone());
            println!("offset_path: {}", offset_path);
            println!("value_path: {}", value_path);
            println!("lookup_path: {}", lookup_path);

            // Load index table
            let value_vec = measure_time!(read_u64_vector(&value_path).expect("[ERROR] Failed to load value vector"));
            let offset_table = measure_time!(read_offset_map(&offset_path).expect("[ERROR] Failed to load offset table"));
            let (path_vec, numeric_id_vec, optional_vec) = measure_time!(load_lookup_from_file(&lookup_path));
            let lookup = measure_time!(load_lookup_from_file(&lookup_path));

            // Make query with pdb
            // P00776,149,171,253
            let pdb_dir = "/fast/hyunbin/motif/swissprot_benchmark/swissprot_v4_raw";
            let pdb_path = "AF-P00766-F1-model_v4.pdb";
            let pdb_path = format!("{}/{}", pdb_dir, pdb_path);
            
            let pdb_query =  make_query(&pdb_path, &vec![149u64, 171u64, 253u64]);

            // Get values with offset
            let offset_to_query = pdb_query.iter().map(|&x| offset_table.get(&x).unwrap()).collect::<Vec<_>>();
            let mut intersection = Vec::new();
            for i in 0..offset_to_query.len() {
                println!("offset_to_query: {:?}", offset_to_query[i]);
                let single_queried_values = get_values_with_offset(&value_vec, offset_to_query[i].0, offset_to_query[i].1);
                if i == 0 {
                    intersection = single_queried_values.to_vec();
                } else {
                    // Filter intersection with single_queried_values
                    intersection = intersection.par_iter().filter(|&&x| single_queried_values.contains(&x)).map(|&x| x).collect::<Vec<_>>();

                }
            }
            println!("queried size: {:?}", intersection.len());
            for i in 0..intersection.len() {
                println!("{:?}", lookup.0[intersection[i] as usize]);
            }
        },
        _ => {
            eprintln!("{}", HELP_QUERY);
            std::process::exit(1);
        }
    }
}