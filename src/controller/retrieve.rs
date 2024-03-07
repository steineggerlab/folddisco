// 
use std::{collections::HashMap, fs::File, hash::Hash};
use crate::prelude::*;

use super::feature::get_single_feature;

pub fn retrieve_residue_with_hash(
    hash_vec: &Vec<GeometricHash>, path: &str, nbin_dist: usize, nbin_angle: usize
) -> Option<Vec<((u8, u8), (u64, u64))>> {
    let file = File::open(path).expect("File not found");
    let pdb = PDBReader::new(file);
    let compact = pdb.read_structure().expect("Error reading structure");
    let compact = compact.to_compact();
    let comb = get_all_combination(compact.num_residues, false);
    let hash_type = hash_vec[0].hash_type();
    let mut output: Vec<((u8, u8), (u64, u64))> = Vec::new();
    for (i, j) in comb.iter() {
        let feature = get_single_feature(*i, *j, &compact, hash_type);
        if feature.is_some() {
            let feature = feature.unwrap();
            let curr_hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(feature, hash_type)
            } else {
                GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
            };
            for hash in hash_vec {
                if curr_hash == *hash {
                    output.push((
                        (compact.chain_per_residue[*i],
                         compact.chain_per_residue[*j]),
                        (compact.residue_serial[*i],
                         compact.residue_serial[*j])
                    ));
                }
            }
        }
    }
    if output.is_empty() {
        None
    } else {
        Some(output)
    }
}

pub fn connected(res_ind_vec: &Vec<((u8, u8), (u64, u64))>, len: usize) -> usize {
    let mut res_set: HashMap<String, usize> = HashMap::new();
    for (i, j) in res_ind_vec {
        let key1 = format!("{}{}", i.0 as char, j.0);
        let key2 = format!("{}{}", i.1 as char, j.1);
        // if key is not in the set, add it or increment the value
        if !res_set.contains_key(&key1) {
            res_set.insert(key1, 1);
        } else {
            let count = res_set.get_mut(&key1).unwrap();
            *count += 1;
        }
        if !res_set.contains_key(&key2) {
            res_set.insert(key2, 1);
        } else {
            let count = res_set.get_mut(&key2).unwrap();
            *count += 1;
        }
    }
    // Count elements that are greater than 1
    let max = res_set.iter().filter(|(_, &v)| v > 1).count();
    max
}

pub fn res_vec_as_string(res_vec: &Vec<((u8, u8), (u64, u64))>) -> String {
    let mut output = String::new();
    // Merge chain and residue number. Commas separate each pair
    for (k, (i, j)) in res_vec.iter().enumerate() {
        // If last element, don't add comma
        if k == res_vec.len() - 1 {
            output.push_str(&format!("{}{}-{}{}", i.0 as char, j.0, i.1 as char, j.1));
        } else {
            output.push_str(&format!("{}{}-{}{},", i.0 as char, j.0, i.1 as char, j.1));
        }
    }
    output
}