// File: query.rs
// Created: 2023-12-22 17:00:50
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use crate::geometry::core::{GeometricHash, HashType};
use crate::prelude::{print_log_msg, PDBReader, INFO};
use crate::utils::log::{log_msg, FAIL};
use super::feature::get_single_feature;
// TODO: Make this to handle multiple hash types
// TODO: Query should be able to consider both chain & residue index


// Query is expected to be given as a path, and a list of tuples of chain and residue index

pub fn make_query(
    path: &String, query_residues: &Vec<(u8, u64)>, hash_type: HashType, 
    nbin_dist: usize, nbin_angle: usize, check_nearby: bool
) -> Vec<GeometricHash> {
    let pdb_reader = PDBReader::from_file(path).expect("PDB file not found");
    let compact = pdb_reader.read_structure().expect("Failed to read PDB file");
    let compact = compact.to_compact();
    
    let mut hash_collection = Vec::new();
    
    // Convert residue indices to vector indices
    let mut indices = Vec::new();
    let mut query_residues = query_residues.clone();
    if query_residues.is_empty() {
        // Iterate over all residues and set to query_residues
        for i in 0..compact.num_residues {
            let chain = compact.chain_per_residue[i];
            let residue_index = compact.residue_serial[i];
            query_residues.push((chain, residue_index));
        }
    }

    for (chain, ri) in query_residues {
        let index = compact.get_index(&chain, &ri);
        if let Some(index) = index {
            // convert u8 array to string
            let residue: String = compact.get_res_name(index).iter().map(|&c| c as char).collect();
            print_log_msg(INFO, &format!("Found residue: {}:{:?}({})", chain as char, ri, residue));
            indices.push(index);
        }
    }
    // Make combinations
    for i in 0..indices.len() {
        for j in i+1..indices.len() {
            let feature = get_single_feature(
                indices[i], indices[j], &compact, hash_type
            );
            if feature.is_some() {
                let feature = feature.unwrap();
                if check_nearby {
                    let mut feature_near = feature.clone();
                    let mut feature_far = feature.clone();
                    feature_near[2] -= 0.5; // TODO: need to be improved
                    feature_far[2] += 0.5;
                    if nbin_dist == 0 || nbin_angle == 0 {
                        let hash_near = GeometricHash::perfect_hash_default(feature_near, hash_type);
                        let hash_far = GeometricHash::perfect_hash_default(feature_far, hash_type);
                        hash_collection.push(hash_near);
                        hash_collection.push(hash_far);
                    } else {
                        let hash_near = GeometricHash::perfect_hash(feature_near, hash_type, nbin_dist, nbin_angle);
                        let hash_far = GeometricHash::perfect_hash(feature_far, hash_type, nbin_dist, nbin_angle);
                        hash_collection.push(hash_near);
                        hash_collection.push(hash_far);
                    }
                } else {
                    if nbin_dist == 0 || nbin_angle == 0 {
                        let hash_value = GeometricHash::perfect_hash_default(feature, hash_type);
                        hash_collection.push(hash_value);
                    } else {
                        let hash_value = GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle);
                        hash_collection.push(hash_value);
                    }
                }
            }
        }
    }
    let mut hash_collection = hash_collection;
    hash_collection.sort_unstable();
    hash_collection.dedup();
    hash_collection
}

pub fn parse_query_string(query_string: &str) -> Vec<(u8, u64)> {
    let mut query_residues: Vec<(u8, u64)> = Vec::new();
    if query_string.is_empty() {
        return query_residues;
    }
    let mut chain = b' ';
    let mut residue = String::new();
    for c in query_string.chars() {
        if c.is_ascii_alphabetic() {
            chain = c as u8;
        } else if c.is_ascii_digit() {
            residue.push(c);
        } else if c == ',' {
            let res_u64 = residue.parse::<u64>().expect(
                &log_msg(FAIL,  "Failed to parse residue")
            );
            if chain == b' ' {
                chain = b'A';
            }
            query_residues.push((chain, res_u64));
            // Reset
            chain = b' ';
            residue.clear();
        } else if c == ' ' {
            continue;
        } else {
            panic!("Invalid character in query string");
        }
    }
    // Push the last residue
    if chain == b' ' {
        chain = b'A';
    }
    let res_u64 = residue.parse::<u64>().expect("Failed to parse residue");
    query_residues.push((chain, res_u64));
    query_residues
}

// ADD TEST
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_make_query() {
        let path = String::from("data/serine_peptidases_filtered/1aq2.pdb");
        let query_residues = vec![
            (b'A', 250), (b'A', 232), (b'A', 269)
        ];
        let hash_type = HashType::FoldDiscoDefault;
        let hash_collection = make_query(&path, &query_residues, hash_type, 0, 0, false);
        println!("{:?}", hash_collection);
    }
    #[test]
    fn test_parse_query_string() {
        let query_string = "A250,B232,C269";
        let query_residues = parse_query_string(query_string);
        assert_eq!(query_residues, vec![(b'A', 250), (b'B', 232), (b'C', 269)]);
    }
    #[test]
    fn test_parse_query_string_with_space() {
        let query_string = "A250, A232, A269";
        let query_residues = parse_query_string(query_string);
        assert_eq!(query_residues, vec![(b'A', 250), (b'A', 232), (b'A', 269)]);
    }
    
    #[test]
    fn test_parse_query_string_with_space_and_no_chain() {
        let query_string = "250, 232, 269";
        let query_residues = parse_query_string(query_string);
        assert_eq!(query_residues, vec![(b'A', 250), (b'A', 232), (b'A', 269)]);
    }
}