// File: query_pdb.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
// Description
// This file contains the workflow for querying PDB files
// When querying PDB files, we need index table and query file.

use crate::cli::*;
use crate::controller::io::{read_offset_map, read_usize_vector, get_values_with_offset};
use crate::index::lookup::{load_lookup_from_file};
use rayon::prelude::*;
use crate::prelude::*;


pub const HELP_QUERY: &str = "\
USAGE: motifsearch query [OPTIONS] <QUERY_PDB> <CHAIN1><RES1>,<CHAIN2><RES2>,<CHAIN3><RES3>...
Example: motifsearch query -i index_table.index -t 6 1aq2.pdb A250,A232,A269
Options:
    -d, --pdb <PDB_PATH>          Path of PDB file to query
    -q, --query <QUERY_STRING>  Query string
    -i, --index <INDEX_PATH>      Path of index table to load
    -t, --threads <THREADS>            Number of threads to use
    -v, --verbose                      Print verbose messages
    -h, --help                         Print this help menu
";

pub fn query_pdb(env: AppArgs) {
    match env {
        AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            help,
        } => {
            if help {
                println!("THIS is query1");
                eprintln!("{}", HELP_QUERY);
                std::process::exit(0);
            }
            // Check if arguments are valid
            if index_path.is_none() {
                                println!("THIS is query2");
                eprintln!("{}", HELP_QUERY);
                std::process::exit(1);
            }
            
            // Get path. formatting without quotation marks
            let index_path = index_path.unwrap();
            let offset_path = format!("{}.offset", index_path.clone()); 
            let value_path = format!("{}.value", index_path.clone());
            let lookup_path = format!("{}.lookup", index_path.clone());
            let hash_type_path = format!("{}.type", index_path.clone());
            let hash_type = HashType::load_from_file(&hash_type_path);
            
            // Load index table
            let (mmap, value_vec) = measure_time!(read_usize_vector(&value_path).expect("[ERROR] Failed to load value vector"));
            let offset_table = measure_time!(
                read_offset_map(&offset_path, hash_type).expect("[ERROR] Failed to load offset table")
            );
            let (path_vec, numeric_id_vec, optional_vec) = measure_time!(load_lookup_from_file(&lookup_path));
            let lookup = measure_time!(load_lookup_from_file(&lookup_path));

            // Make query with pdb
            let query_residues = parse_query_string(&query_string);
            let pdb_query =  make_query(
                &pdb_path, &query_residues, hash_type
            );

            // Get values with offset
            // let offset_to_query = pdb_query.iter().map(|&x| offset_table.get(&x).unwrap()).collect::<Vec<_>>();
            
            // Get offset from offset_table with query
            let mut offset_to_query = Vec::new();
            for i in 0..pdb_query.len() {
                let offset = offset_table.get(&pdb_query[i]).unwrap();
                // Get offset map values
                offset_to_query.push(*offset);
            }
            
            println!("{:?}", offset_to_query);
            let mut intersection = Vec::new();
            for i in 0..offset_to_query.len() {
                let single_queried_values = get_values_with_offset(&value_vec, offset_to_query[i].0, offset_to_query[i].1);
                println!("{:?}", single_queried_values);
                if i == 0 {
                    intersection = single_queried_values.to_vec();
                } else {
                    // Filter intersection with single_queried_values
                    intersection = intersection.par_iter().filter(
                        |&&x| single_queried_values.contains(&x)
                    ).map(|&x| x).collect::<Vec<_>>();
                }
            }

            for i in 0..intersection.len() {
                let nid = lookup.1[intersection[i] as usize];
                println!("{:?}", lookup.0[nid]);
            }
        },
        _ => {
            eprintln!("{}", HELP_QUERY);
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_query_pdb_workflow() {
        let pdb_path = String::from("data/serine_peptidases_filtered/1aq2.pdb");
        let query_string = String::from("A250,A232,A269");
        let threads = 6;
        let index_path = Some(String::from("data/index/serine_peptidases_filtered"));
        let help = false;
        let env = AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            help,
        };
        query_pdb(env);
    }
}