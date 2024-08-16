// File: query.rs
// Created: 2023-12-22 17:00:50
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use std::collections::HashMap;

use crate::geometry::core::{GeometricHash, HashType};
use crate::utils::convert::{is_aa_group_char, map_one_letter_to_u8_vec};
use crate::prelude::{print_log_msg, PDBReader, INFO};
use crate::utils::combination::CombinationIterator;
use crate::utils::log::{log_msg, FAIL};
use super::feature::get_single_feature;
use super::io::read_compact_structure;

// Query is expected to be given as a path, and a list of tuples of chain and residue index
pub fn make_query(
    path: &String, query_residues: &Vec<(u8, u64)>, hash_type: HashType, 
    nbin_dist: usize, nbin_angle: usize, exact_match: bool, dist_thresholds: Vec<f32>, angle_thresholds: Vec<f32>,
    amino_acid_substitutions: &Vec<Option<Vec<u8>>>,
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

    let mut substitution_map: HashMap<usize, Vec<u8>> = HashMap::new();
    
    for (i, (chain, ri)) in query_residues.iter().enumerate() {
        let index = compact.get_index(&chain, &ri);
        if let Some(index) = index {
            // convert u8 array to string
            let _residue: String = compact.get_res_name(index).iter().map(|&c| c as char).collect();
            indices.push(index);
            if let Some(substitution) = amino_acid_substitutions[i].clone() {
                substitution_map.insert(index, substitution);
            }
        }
    }
    let dist_indices = hash_type.dist_index();
    let angle_indices = hash_type.angle_index();
    // Make combinations
    let comb_iter = CombinationIterator::new(indices.len());
    comb_iter.for_each(|(i, j)| {
        if i == j {
            return;
        }
        let feature = get_single_feature(
            indices[i], indices[j], &compact, hash_type
        );
        if feature.is_some() {
            let feature = feature.unwrap();
            if !exact_match & (dist_thresholds.len() > 0 || angle_thresholds.len() > 0) {
                if let Some(dist_indices) = &dist_indices {
                    for dist_threshold in dist_thresholds.iter() {
                        for dist_index in dist_indices.iter() {
                            let mut feature_near = feature.clone();
                            feature_near[*dist_index] -= dist_threshold;
                            let mut feature_far = feature.clone();
                            feature_far[*dist_index] += dist_threshold;
                            append_hash(
                                nbin_dist, nbin_angle, feature_near, hash_type, feature_far, &mut hash_collection
                            );
                        }
                    }
                }
                if angle_indices.is_some() {
                    if let Some(angle_indices) = &angle_indices {
                        for angle_threshold in angle_thresholds.iter() {
                            for angle_index in angle_indices.iter() {
                                let mut feature_near = feature.clone();
                                let angle_threshold = angle_threshold.to_radians();
                                feature_near[*angle_index] -= angle_threshold;
                                let mut feature_far = feature.clone();
                                feature_far[*angle_index] += angle_threshold;
                                append_hash(
                                    nbin_dist, nbin_angle, feature_near, hash_type, feature_far, &mut hash_collection
                                );
                            }
                        }
                    }
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
    });

    let mut hash_collection = hash_collection;
    hash_collection.sort_unstable();
    hash_collection.dedup();
    hash_collection
}

fn append_hash(
    nbin_dist: usize, nbin_angle: usize, feature_near: Vec<f32>, hash_type: HashType, feature_far: Vec<f32>,
    hash_collection: &mut Vec<GeometricHash>
) {
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
}

pub fn parse_threshold_string(threshold_string: Option<String>) -> Vec<f32> {
    if threshold_string.is_none() {
        return Vec::new();
    }
    let threshold_string = threshold_string.unwrap();
    // Remove whitespace
    let threshold_string = threshold_string.replace(" ", "");
    let mut thresholds: Vec<f32> = Vec::new();
    for threshold in threshold_string.split(',') {
        let threshold = threshold.parse::<f32>().expect(
            &log_msg(FAIL, "Failed to parse threshold")
        );
        thresholds.push(threshold);
    }
    thresholds
}


pub fn make_query_map(
    path: &String, query_residues: &Vec<(u8, u64)>, hash_type: HashType, 
    nbin_dist: usize, nbin_angle: usize, dist_thresholds: &Vec<f32>, angle_thresholds: &Vec<f32>,
    amino_acid_substitutions: &Vec<Option<Vec<u8>>>,
) -> (HashMap<GeometricHash, ((usize, usize), bool)>, Vec<usize>, HashMap<(u8, u8), Vec<(f32, usize)>>) {

    let (compact, _) = read_compact_structure(path).expect("Failed to read compact structure");
    
    let mut hash_collection = HashMap::new();
    let mut observed_distance_map: HashMap<(u8, u8), Vec<(f32, usize)>> = HashMap::new();
    
    // Convert residue indices to vector indices
    let mut indices = Vec::new();
    let mut query_residues = query_residues.clone();
    let mut amino_acid_substitutions = amino_acid_substitutions.clone();

    if query_residues.is_empty() {
        // Iterate over all residues and set to query_residues
        for i in 0..compact.num_residues {
            let chain = compact.chain_per_residue[i];
            let residue_index = compact.residue_serial[i];
            query_residues.push((chain, residue_index));
            amino_acid_substitutions.push(None);
        }
    }

    let mut substitution_map: HashMap<usize, Vec<u8>> = HashMap::new();
    
    for (i, (chain, ri)) in query_residues.iter().enumerate() {
        let index = compact.get_index(&chain, &ri);
        if let Some(index) = index {
            // convert u8 array to string
            let _residue: String = compact.get_res_name(index).iter().map(|&c| c as char).collect();
            indices.push(index);
            if let Some(substitution) = amino_acid_substitutions[i].clone() {
                substitution_map.insert(index, substitution);
            }
        }
    }
    let dist_indices = hash_type.dist_index();
    let angle_indices = hash_type.angle_index();
    // Make combinations
    let comb_iter = CombinationIterator::new(indices.len());
    comb_iter.for_each(|(i, j)| {
        if i == j {
            return;
        }
        let feature = get_single_feature(
            indices[i], indices[j], &compact, hash_type
        );
        
        if feature.is_some() {
            // Gather observed distance & aa pairs.
            let aa_dist_info = compact.get_list_amino_acids_and_distances(indices[i], indices[j]);
            if let Some(aa_dist_info) = aa_dist_info {
                // Check if the pair is already in the map
                let aa_pair = (aa_dist_info.0, aa_dist_info.1);
                if observed_distance_map.contains_key(&aa_pair) {
                    observed_distance_map.get_mut(&aa_pair).unwrap().push((aa_dist_info.2, indices[i]));
                } else {
                    observed_distance_map.insert(aa_pair, vec![(aa_dist_info.2, indices[i])]);
                }
            }

            let feature = feature.unwrap();
            let feature_clone = feature.clone();
            if nbin_dist == 0 || nbin_angle == 0 {
                let hash_value = GeometricHash::perfect_hash_default(feature_clone, hash_type);
                hash_collection.insert(hash_value, ((indices[i], indices[j]), true));
            } else {
                let hash_value = GeometricHash::perfect_hash(feature_clone, hash_type, nbin_dist, nbin_angle);
                hash_collection.insert(hash_value, ((indices[i], indices[j]), true));
            }
            // Get zip of substitution of i and j
            if let Some(aa_index) = hash_type.amino_acid_index() {
                if let Some(sub_i) = substitution_map.get(&indices[i]) {
                    // Substitute amino acid of i only
                    for sub_i in sub_i.iter() {
                        let mut feature = feature.clone();
                        feature[aa_index[0]] = *sub_i as f32;
                        let hash_value = if nbin_dist == 0 || nbin_angle == 0 {
                            GeometricHash::perfect_hash_default(feature, hash_type)
                        } else {
                            GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
                        };
                        if hash_collection.contains_key(&hash_value) {
                            continue;
                        } else {
                            hash_collection.insert(hash_value, ((indices[i], indices[j]), false));
                        }
                    }
                    // Substitute amino acid of both i and j
                    if let Some(sub_j) = substitution_map.get(&indices[j]) {
                        for sub_i in sub_i.iter() {
                            for sub_j in sub_j.iter() {
                                let mut feature = feature.clone();
                                feature[aa_index[0]] = *sub_i as f32;
                                feature[aa_index[1]] = *sub_j as f32;
                                let hash_value = if nbin_dist == 0 || nbin_angle == 0 {
                                    GeometricHash::perfect_hash_default(feature, hash_type)
                                } else {
                                    GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
                                };
                                if hash_collection.contains_key(&hash_value) {
                                    continue;
                                } else {
                                    hash_collection.insert(hash_value, ((indices[i], indices[j]), false));
                                }
                            }
                        }
                    }
                } else {
                    if let Some(sub_j) = substitution_map.get(&indices[j]) {
                        // Substitute amino acid of j only
                        for sub_j in sub_j.iter() {
                            let mut feature = feature.clone();
                            feature[aa_index[1]] = *sub_j as f32;
                            let hash_value = if nbin_dist == 0 || nbin_angle == 0 {
                                GeometricHash::perfect_hash_default(feature, hash_type)
                            } else {
                                GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
                            };
                            if hash_collection.contains_key(&hash_value) {
                                continue;
                            } else {
                                hash_collection.insert(hash_value, ((indices[i], indices[j]), false));
                            }
                        }
                    }
                }
            }

            if let Some(dist_indices) = &dist_indices {
                for dist_threshold in dist_thresholds.iter() {
                    for dist_index in dist_indices.iter() {
                        let mut feature_near = feature.clone();
                        feature_near[*dist_index] -= dist_threshold;
                        let mut feature_far = feature.clone();
                        feature_far[*dist_index] += dist_threshold;
                        if nbin_dist == 0 || nbin_angle == 0 {
                            let hash_near = GeometricHash::perfect_hash_default(feature_near, hash_type);
                            let hash_far = GeometricHash::perfect_hash_default(feature_far, hash_type);
                            if hash_collection.contains_key(&hash_near) {
                                continue;
                            } else {
                                hash_collection.insert(hash_near, ((indices[i], indices[j]), false));
                            }
                            if hash_collection.contains_key(&hash_far) {
                                continue;
                            } else {
                                hash_collection.insert(hash_far, ((indices[i], indices[j]), false));
                            }
                        } else {
                            let hash_near = GeometricHash::perfect_hash(feature_near, hash_type, nbin_dist, nbin_angle);
                            let hash_far = GeometricHash::perfect_hash(feature_far, hash_type, nbin_dist, nbin_angle);
                            if hash_collection.contains_key(&hash_near) {
                                continue;
                            } else {
                                hash_collection.insert(hash_near, ((indices[i], indices[j]), false));
                            }
                            if hash_collection.contains_key(&hash_far) {
                                continue;
                            } else {
                                hash_collection.insert(hash_far, ((indices[i], indices[j]), false));
                            }
                        }
                    }
                }
            }
            if let Some(angle_indices) = &angle_indices {
                for angle_threshold in angle_thresholds.iter() {
                    for angle_index in angle_indices.iter() {
                        let mut feature_near = feature.clone();
                        let angle_threshold = angle_threshold.to_radians();
                        feature_near[*angle_index] -= angle_threshold;
                        let mut feature_far = feature.clone();
                        feature_far[*angle_index] += angle_threshold;
                        if nbin_dist == 0 || nbin_angle == 0 {
                            let hash_near = GeometricHash::perfect_hash_default(feature_near, hash_type);
                            let hash_far = GeometricHash::perfect_hash_default(feature_far, hash_type);
                            if hash_collection.contains_key(&hash_near) {
                                continue;
                            } else {
                                hash_collection.insert(hash_near, ((indices[i], indices[j]), false));
                            }
                            if hash_collection.contains_key(&hash_far) {
                                continue;
                            } else {
                                hash_collection.insert(hash_far, ((indices[i], indices[j]), false));
                            }
                        } else {
                            let hash_near = GeometricHash::perfect_hash(feature_near, hash_type, nbin_dist, nbin_angle);
                            let hash_far = GeometricHash::perfect_hash(feature_far, hash_type, nbin_dist, nbin_angle);
                            if hash_collection.contains_key(&hash_near) {
                                continue;
                            } else {
                                hash_collection.insert(hash_near, ((indices[i], indices[j]), false));
                            }
                            if hash_collection.contains_key(&hash_far) {
                                continue;
                            } else {
                                hash_collection.insert(hash_far, ((indices[i], indices[j]), false));
                            }
                        }
                    }
                }
            }

        }
    });
    (hash_collection, indices, observed_distance_map)
}


pub fn parse_query_string(query_string: &str, default_chain: u8) -> (Vec<(u8, u64)>, Vec<Option<Vec<u8>>>) {
    let mut query_residues: Vec<(u8, u64)> = Vec::new();
    let mut amino_acid_substitutions: Vec<Option<Vec<u8>>> = Vec::new();
    if query_string.is_empty() {
        return (query_residues, amino_acid_substitutions);
    }
    let mut default_chain = default_chain;
    // If default chain is not alphabetic, just use A as default
    if !default_chain.is_ascii_alphabetic() {
        default_chain = b'A';
    }

    let mut chain = b' ';
    let mut residue = String::new();
    let mut residue_substitution: Option<Vec<u8>> = None;
    for c in query_string.chars() {
        if is_aa_group_char(c) {
            if residue.is_empty() {
                chain = c as u8;
            } else {
                let curr_sub_aa = map_one_letter_to_u8_vec(c);
                if residue_substitution.is_some() {
                    residue_substitution.as_mut().unwrap().extend(curr_sub_aa);
                } // If None, it's a problem. Don't consider the case. 
            }
        } else if c.is_ascii_digit() {
            residue.push(c);
        } else if c == ':' {
            residue_substitution = Some(Vec::new());
        } else if c == ',' {
            let res_u64 = residue.parse::<u64>().expect(
                &log_msg(FAIL,  "Failed to parse residue")
            );
            if chain == b' ' {
                chain = default_chain;
            }
            query_residues.push((chain, res_u64));
            amino_acid_substitutions.push(residue_substitution);
            // Reset
            chain = b' ';
            residue.clear();
            residue_substitution = None;
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
    amino_acid_substitutions.push(residue_substitution);
    assert_eq!(query_residues.len(), amino_acid_substitutions.len());
    (query_residues, amino_acid_substitutions)
}

pub fn get_offset_value_lookup_type(index_path: String) -> (String, String, String, String) {
    let offset_path = format!("{}.offset", index_path.clone());
    let value_path = format!("{}.value", index_path.clone());
    let lookup_path = format!("{}.lookup", index_path.clone());
    let hash_type_path = format!("{}.type", index_path.clone());
    assert!(std::path::Path::new(&offset_path).is_file());
    assert!(std::path::Path::new(&value_path).is_file());
    assert!(std::path::Path::new(&lookup_path).is_file());
    assert!(std::path::Path::new(&hash_type_path).is_file());
    (offset_path, value_path, lookup_path, hash_type_path)
}

pub fn check_and_get_indices(index_path: Option<String>, verbose: bool) -> Vec<String> {
    // Get path. formatting without quotation marks
    let index_path = index_path.unwrap();
    // Check if index_path_0 is a file.
    let _index_chunk_prefix = format!("{}_0", index_path.clone());
    let index_chunk_path = format!("{}_0.offset", index_path.clone());
    let mut index_paths = Vec::new();
    if std::path::Path::new(&index_chunk_path).is_file() {
        if verbose {
            print_log_msg(INFO, &format!("Index table is chunked"));
        }
        let mut i = 0;
        loop {
            let index_chunk_prefix = format!("{}_{}", index_path.clone(), i);
            let index_chunk_path = format!("{}.offset", index_chunk_prefix);
            if std::path::Path::new(&index_chunk_path).is_file() {
                index_paths.push(index_chunk_prefix);
                i += 1;
            } else {
                break;
            }
        }
    } else {
        index_paths.push(index_path.clone());
    }
    index_paths
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
        let aa_substitution = vec![Some(vec![1]), Some(vec![11]), Some(vec![5])];
        let dist_thresholds: Vec<f32> = vec![0.5];
        let angle_thresholds: Vec<f32> = vec![5.0];
        let hash_type = HashType::PDBMotifSinCos;
        let hash_collection = make_query(
            &path, &query_residues, hash_type, 8, 6, false, dist_thresholds, angle_thresholds, &aa_substitution
        );
        // let hash_collection = make_query(&path, &query_residues, hash_type, 8, 3, true, dist_thresholds, angle_thresholds);
        println!("{:?}", hash_collection);
        println!("{}", hash_collection.len());
    }
    
    #[test]
    fn test_make_query_map() {
        // let path = String::from("data/serine_peptidases_filtered/1aq2.pdb");
        // let query_residues = vec![
        //     (b'A', 250), (b'A', 232), (b'A', 269)
        // ];
        let path = String::from("data/serine_peptidases_filtered/4cha.pdb");
        let query_residues = vec![
            (b'B', 57), (b'B', 102), (b'C', 195)
        ];
        let dist_thresholds: Vec<f32> = vec![0.5];
        let angle_thresholds: Vec<f32> = vec![5.0,10.0,15.0];
        let hash_type = HashType::PDBMotifSinCos;
        let (hash_collection, _index_found, _observed_dist_map) = make_query_map(
            &path, &query_residues, hash_type, 8, 3, &dist_thresholds, &angle_thresholds, &vec![None, None, None]
        );
        println!("{:?}", hash_collection);
        println!("{}", hash_collection.len());
        println!("{:?}", _observed_dist_map);
        // Print the count where value.1 is true
        let mut count = 0;
        hash_collection.iter().for_each(|item| {
            if item.1.1 {
                count += 1;
            }
        });
        println!("Exact: {}", count);
        println!("Not exact: {}", hash_collection.len() - count);
    }

    #[test]
    fn test_parse_query_string() {
        let query_string = "A250,B232,C269";
        let query_residues = parse_query_string(query_string, b'A');
        assert_eq!(query_residues, (vec![(b'A', 250), (b'B', 232), (b'C', 269)], vec![None, None, None]));
    }
    #[test]
    fn test_parse_query_string_with_space() {
        let query_string = "A250, A232, A269";
        let query_residues = parse_query_string(query_string, b'A');
        assert_eq!(query_residues, (vec![(b'A', 250), (b'A', 232), (b'A', 269)], vec![None, None, None]));
    }
    
    #[test]
    fn test_parse_query_string_with_space_and_no_chain() {
        let query_string = "250, 232, 269";
        let query_residues = parse_query_string(query_string, b'A');
        assert_eq!(query_residues, (vec![(b'A', 250), (b'A', 232), (b'A', 269)], vec![None, None, None]));
    }

    #[test]
    fn test_parse_query_string_with_aa_substitution() {
        let query_string = "A250:R,B232:K,C269:Q";
        let query_residues = parse_query_string(query_string, b'A');
        // R = 1, K = 11, Q = 5
        assert_eq!(query_residues, (vec![(b'A', 250), (b'B', 232), (b'C', 269)], vec![Some(vec![1]), Some(vec![11]), Some(vec![5])]));
        let query_string = "250:R,232:K,269:Q";
        let query_residues = parse_query_string(query_string, b'A');
        // R = 1, K = 11, Q = 5
        assert_eq!(query_residues, (vec![(b'A', 250), (b'A', 232), (b'A', 269)], vec![Some(vec![1]), Some(vec![11]), Some(vec![5])]));
    }
}