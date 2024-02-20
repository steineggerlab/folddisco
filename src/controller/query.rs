// File: query.rs
// Created: 2023-12-22 17:00:50
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use crate::{prelude::{GeometricHash, HashType}, PDBReader};
use super::feature::get_single_feature;
// TODO: Make this to handle multiple hash types
// TODO: Query should be able to consider both chain & residue index


// Query is expected to be given as a path, and a list of tuples of chain and residue index

pub fn make_query(path: &String, query_residues: &Vec<(u8, u64)>, hash_type: HashType) -> Vec<GeometricHash> {
    let pdb_reader = PDBReader::from_file(path).expect("PDB file not found");
    let compact = pdb_reader.read_structure().expect("Failed to read PDB file");
    let compact = compact.to_compact();
    
    let mut hash_collection = Vec::new();
    
    // Convert residue indices to vector indices
    let mut indices = Vec::new();
    for (chain, ri) in query_residues {
        let index = compact.get_index(chain, ri);
        if let Some(index) = index {
            indices.push(index);
        }
    }
    // Make combinations
    for i in 0..indices.len() {
        for j in i+1..indices.len() {
            let feature = get_single_feature(
                indices[i], indices[j], &compact, hash_type
            );
            let hash_value = GeometricHash::perfect_hash(feature, hash_type);
            hash_collection.push(hash_value);
        }
    }
    hash_collection
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
        let hash_collection = make_query(&path, &query_residues, hash_type);
        println!("{:?}", hash_collection);
    }
}