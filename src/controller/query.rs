// File: query.rs
// Created: 2023-12-22 17:00:50
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use std::collections::HashMap;

use crate::geometry::core::{GeometricHash, HashType};
use crate::utils::convert::{is_aa_group_char, map_one_letter_to_u8_vec};
use crate::utils::combination::CombinationIterator;
use crate::utils::log::{log_msg, FAIL};
use super::feature::get_single_feature;
use super::io::read_compact_structure;

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

// Doesn't support duplicate hash
// If hash is already in the hash_collection, skip.
fn insert_binned_hash(
    hash_collection: &mut HashMap<GeometricHash, ((usize, usize), bool)>,
    feature: &Vec<f32>, indices: (usize, usize), hash_type: HashType,
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
    is_primary: bool,
) {
    if let Some(multiple_bin) = multiple_bin {
        for (nbin_dist, nbin_angle) in multiple_bin.iter() {
            let hash_value = if *nbin_dist == 0 || *nbin_angle == 0 {
                GeometricHash::perfect_hash_default(feature, hash_type)
            } else {
                GeometricHash::perfect_hash(feature, hash_type, *nbin_dist, *nbin_angle)
            };
            if hash_collection.contains_key(&hash_value) {
                continue;
            } else {
                hash_collection.insert(hash_value, (indices, is_primary));
            }
        }
    } else {
        let hash_value = if nbin_dist == 0 || nbin_angle == 0 {
            GeometricHash::perfect_hash_default(feature, hash_type)
        } else {
            GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
        };
        if hash_collection.contains_key(&hash_value) {
            return;
        } else {
            hash_collection.insert(hash_value, (indices, is_primary));
        }
    }
}

fn apply_substitutions(
    i_idx: usize, j_idx: usize, feature: &mut Vec<f32>, aa_index: &[usize],
    nbin_dist: usize, nbin_angle: usize,
    multiple_bin: &Option<Vec<(usize, usize)>>,
    substitution_map: &HashMap<usize, Vec<u8>>,
    hash_collection: &mut HashMap<GeometricHash, ((usize, usize), bool)>,
    hash_type: HashType,
) {
    let orig_aa_value1 = feature[aa_index[0]];
    let orig_aa_value2 = feature[aa_index[1]];

    if let Some(sub_i) = substitution_map.get(&i_idx) {
        // Substitute amino acid of i only
        for aa_sub_i in sub_i {
            let mut temp_feature = feature.to_vec();
            temp_feature[aa_index[0]] = *aa_sub_i as f32;
            insert_binned_hash(
                hash_collection,
                &temp_feature,
                (i_idx, j_idx),
                hash_type,
                nbin_dist,
                nbin_angle,
                multiple_bin,
                false,
            );
        }

        // Substitute amino acid of both i and j
        if let Some(sub_j) = substitution_map.get(&j_idx) {
            for aa_sub_i in sub_i {
                for aa_sub_j in sub_j {
                    feature[aa_index[0]] = *aa_sub_i as f32;
                    feature[aa_index[1]] = *aa_sub_j as f32;
                    insert_binned_hash(
                        hash_collection,
                        feature,
                        (i_idx, j_idx),
                        hash_type,
                        nbin_dist,
                        nbin_angle,
                        multiple_bin,
                        false,
                    );
                    // Reset feature
                    feature[aa_index[0]] = orig_aa_value1;
                    feature[aa_index[1]] = orig_aa_value2;
                }
            }
        }
    } else if let Some(sub_j) = substitution_map.get(&j_idx) {
        // Substitute amino acid of j only
        for aa_sub_j in sub_j {
            let mut temp_feature = feature.to_vec();
            temp_feature[aa_index[1]] = *aa_sub_j as f32;
            insert_binned_hash(
                hash_collection,
                &temp_feature,
                (i_idx, j_idx),
                hash_type,
                nbin_dist,
                nbin_angle,
                multiple_bin,
                false,
            );
        }
    }
}

fn expand_and_insert(
    target_feature_indices: &[usize], threshold_vector: &[f32],
    feature_near: &mut Vec<f32>, feature_far: &mut Vec<f32>,
    hash_collection: &mut HashMap<GeometricHash, ((usize, usize), bool)>,
    i_idx: usize, j_idx: usize, hash_type: HashType,
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
) {
    for &delta in threshold_vector {
        // If we need to convert to radians (for angles), do so now
        for &idx in target_feature_indices {
            feature_near[idx] -= delta;
            feature_far[idx] += delta;

            insert_binned_hash(
                hash_collection, feature_near, (i_idx, j_idx),
                hash_type, nbin_dist, nbin_angle, multiple_bin, false,
            );
            insert_binned_hash(
                hash_collection, feature_far, (i_idx, j_idx),
                hash_type, nbin_dist, nbin_angle, multiple_bin, false,
            );
            // Restore the original value
            feature_near[idx] += delta;
            feature_far[idx] -= delta;
        }
    }
}

pub fn make_query_map(
    path: &String, query_residues: &Vec<(u8, u64)>, hash_type: HashType, 
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
    dist_thresholds: &Vec<f32>, angle_thresholds: &Vec<f32>,
    amino_acid_substitutions: &Vec<Option<Vec<u8>>>, distance_cutoff: f32, serial_query: bool,
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
        let index = if serial_query { Some(*ri as usize) } else { compact.get_index(&chain, &ri) };
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
    let mut feature = vec![0.0; 9];
    let mut feature_near = vec![0.0; 9];
    let mut feature_far = vec![0.0; 9];
    comb_iter.for_each(|(i, j)| {
        if i == j {
            return;
        }
        let is_feature = get_single_feature(
            indices[i], indices[j], &compact, hash_type, distance_cutoff, &mut feature
        );
        
        if is_feature {
            for k in 0..feature.len() {
                feature_near[k] = feature[k].clone();
                feature_far[k] = feature[k].clone();
            }
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

            // Insert observed hash
            insert_binned_hash(
                &mut hash_collection, &feature, (indices[i], indices[j]), 
                hash_type, nbin_dist, nbin_angle, multiple_bin, true
            );

            // Apply substitutions
            if let Some(aa_indices) = hash_type.amino_acid_index() {
                apply_substitutions(
                    indices[i], indices[j], &mut feature_near, aa_indices.as_ref(),
                    nbin_dist, nbin_angle, multiple_bin, &substitution_map, &mut hash_collection, hash_type
                );
            }
            // apply_substitutions(
            //     indices[i], indices[j], &mut feature_near, hash_type.amino_acid_index().as_ref().unwrap(),
            //     nbin_dist, nbin_angle, multiple_bin, &substitution_map, &mut hash_collection, hash_type
            // );

            if let Some(dist_indices) = &dist_indices {
                expand_and_insert(
                    dist_indices, dist_thresholds, &mut feature_near, &mut feature_far,
                    &mut hash_collection, indices[i], indices[j], hash_type, nbin_dist, nbin_angle, multiple_bin
                );
            }
            if let Some(angle_indices) = &angle_indices {
                let angle_thresholds_in_radian = angle_thresholds.iter().map(|&x| x.to_radians()).collect::<Vec<f32>>();
                expand_and_insert(
                    angle_indices, &angle_thresholds_in_radian, &mut feature_near, &mut feature_far,
                    &mut hash_collection, indices[i], indices[j], hash_type, nbin_dist, nbin_angle, multiple_bin
                );
            }
        }
    });
    (hash_collection, indices, observed_distance_map)
}

pub fn parse_query_string(query_string: &str, mut default_chain: u8) -> (Vec<(u8, u64)>, Vec<Option<Vec<u8>>>) {
    let mut query_residues = Vec::new();
    let mut amino_acid_substitutions = Vec::new();

    if query_string.is_empty() {
        return (query_residues, amino_acid_substitutions);
    }
    if !default_chain.is_ascii_alphabetic() {
        default_chain = b'A';
    }
    // Remove whitespace
    let query_string = query_string.replace(" ", "");
    for segment in query_string.split(',') {
        let (chain, rest) = if let Some(first) = segment.chars().next() {
            // NOTE: 2025-01-15 15:55:19
            // Current querying doesn't support chain ID with more than 1 character
            if first.is_ascii_alphabetic() {
                (first as u8, &segment[1..])
            } else {
                (default_chain, segment)
            }
        } else {
            (default_chain, segment)
        };

        let (range_part, subst_part) = match rest.split_once(':') {
            Some((r, s)) => {
                let sub_vec = s
                    .chars()
                    .filter(|c| is_aa_group_char(*c))
                    .flat_map(|c| map_one_letter_to_u8_vec(c))
                    .collect::<Vec<_>>();
                (r, Some(sub_vec))
            }
            None => (rest, None),
        };

        if range_part.contains('-') {
            let (start_str, end_str) = range_part.split_once('-').expect("Invalid range");
            let start = start_str.parse::<u64>().expect("Invalid start residue");
            let end = end_str.parse::<u64>().expect("Invalid end residue");
            for r in start..=end {
                query_residues.push((chain, r));
                amino_acid_substitutions.push(subst_part.clone());
            }
        } else {
            let residue_num = range_part.parse::<u64>().expect("Invalid residue");
            query_residues.push((chain, residue_num));
            amino_acid_substitutions.push(subst_part);
        }
    }

    (query_residues, amino_acid_substitutions)
}



// ADD TEST
#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_make_query_map() {
        let path= String::from("query/1G2F.pdb");
        let query_residues = vec![
            (b'F', 207), (b'F', 212), (b'F', 225)
        ];
        let amino_acid_substitutions = vec![None; query_residues.len()];
        // let path = String::from("data/serine_peptidases/1aq2.pdb");
        // let query_residues = vec![
        //     (b'A', 250), (b'A', 232), (b'A', 269)
        // ];
        let hash_type = HashType::PDBTrRosetta;
        let (hash_collection, _index_found, _observed_dist_map) = make_query_map(
            &path, &query_residues, hash_type, 16, 4, &None,
            &vec![0.0], &vec![0.0], &amino_acid_substitutions, 20.0, false
        );
        let hash_key = hash_collection.keys().cloned().collect::<Vec<GeometricHash>>();
        println!("{}", hash_collection.len());
        println!("{:?}", _observed_dist_map);
        println!("{:?}", hash_key);
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
        let query_string = "A250:R,B232:K,C269:QK";
        let query_residues = parse_query_string(query_string, b'A');
        // R = 1, K = 11, Q = 5
        assert_eq!(query_residues, (vec![(b'A', 250), (b'B', 232), (b'C', 269)], vec![Some(vec![1]), Some(vec![11]), Some(vec![5, 11])]));
        let query_string = "250:R,232:K,269:QK";
        let query_residues = parse_query_string(query_string, b'A');
        // R = 1, K = 11, Q = 5
        assert_eq!(query_residues, (vec![(b'A', 250), (b'A', 232), (b'A', 269)], vec![Some(vec![1]), Some(vec![11]), Some(vec![5, 11])]));
    }
    #[test]
    fn test_parse_query_string_with_range() {
        let query_string = "A250-252,B232-234,C269:Q";
        let query_residues = parse_query_string(query_string, b'A');
        assert_eq!(query_residues, (vec![
            (b'A', 250), (b'A', 251), (b'A', 252), 
            (b'B', 232), (b'B', 233), (b'B', 234), 
            (b'C', 269),
        ], vec![None, None, None, None, None, None, Some(vec![5])]));
    }
}