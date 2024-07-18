// File: feature.rs
// Created: 2024-03-29 19:05:22
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use crate::utils::convert::map_aa_to_u8;
use crate::structure::core::CompactStructure;
use crate::geometry::core::{GeometricHash, HashType};
use crate::structure::grid::{get_grid_index_vector_from_compact, merge_id_with_grid, nearby};
use crate::utils::combination::CombinationIterator;

pub fn get_single_feature(i: usize, j: usize, structure: &CompactStructure, hash_type: HashType) -> Option<Vec<f32>> {
    let res1 = structure.get_res_name(i);
    let res2 = structure.get_res_name(j);
    let res1 = map_aa_to_u8(res1) as f32;
    let res2 = map_aa_to_u8(res2) as f32;
    if res1 == 255.0 || res2 == 255.0 {
        return None;
    }
    match &hash_type {
        HashType::PDBMotif => {
            let ca_dist = structure.get_ca_distance(i, j);
            let cb_dist = structure.get_cb_distance(i, j);
            let ca_cb_angle = structure.get_ca_cb_angle(i, j); // degree
            if ca_dist.is_some() && cb_dist.is_some() && ca_cb_angle.is_some() {
                let feature = vec![
                    res1, res2, ca_dist.unwrap(), cb_dist.unwrap(), ca_cb_angle.unwrap().to_radians()
                ];
                // Ignore distant interactions
                if feature[2] > 20.0 {
                    None    
                } else {
                    Some(feature)
                }
            } else {
                None
            }
        },
        HashType::PDBMotifSinCos | HashType::PDBMotifHalf => {
            let ca_dist = structure.get_ca_distance(i, j);
            let cb_dist = structure.get_cb_distance(i, j);
            let ca_cb_angle = structure.get_ca_cb_angle(i, j);
            if ca_dist.is_some() && cb_dist.is_some() && ca_cb_angle.is_some() {
                let ca_cb_angle_radian = ca_cb_angle.unwrap().to_radians();
                let feature = vec![
                    res1, res2, ca_dist.unwrap(), cb_dist.unwrap(), ca_cb_angle_radian,
                ];
                // Ignore distant interactions
                if feature[2] > 20.0 {
                    None
                } else {
                    Some(feature)
                }
            } else {
                None
            }
        },
        HashType::TrRosetta => {
            let feature = structure.get_default_feature(i, j);
            if feature.is_some() {
                // Concatenate res1 and res2 to the feature
                let mut feature = feature.unwrap();
                if feature[0] > 20.0 {
                    None
                } else {
                    feature.insert(0, res1);
                    feature.insert(1, res2);
                    Some(feature)
                }
            } else {
                None
            }
        },
        HashType::PointPairFeature => {
            let feature = structure.get_ppf(i, j);
            if feature.is_some() {
                let mut feature = feature.unwrap();
                if feature[0] > 20.0 {
                    None
                } else {
                    feature.insert(0, res1);
                    feature.insert(1, res2);
                    Some(feature)
                }
            } else {
                None
            }
        }
        HashType::PDBTrRosetta => {
            let feature = structure.get_pdb_tr_feature(i, j);
            if feature.is_some() {
                let mut feature = feature.unwrap();
                if feature[0] > 20.0 {
                    None
                } else {
                    feature.insert(0, res1);
                    feature.insert(1, res2);
                    Some(feature)
                }
            } else {
                None
            }
        }
        HashType::TertiaryInteraction => {
            if i == 0 || j == 0 || i == structure.num_residues - 1 || j == structure.num_residues - 1 {
                return None;
            }
            let ca_1i = structure.get_ca(i-1);
            let ca_i = structure.get_ca(i);
            let ca_i1 = structure.get_ca(i+1);
            let ca_1j = structure.get_ca(j-1);
            let ca_j = structure.get_ca(j);
            let ca_j1 = structure.get_ca(j+1);
            
            if ca_1i.is_none() || ca_i.is_none() || ca_i1.is_none() || ca_1j.is_none() || ca_j.is_none() || ca_j1.is_none() {
                return None;
            } else {
                let ca_i = ca_i.unwrap();
                let ca_j = ca_j.unwrap();
                let ca_dist = ca_i.distance(&ca_j);
                if ca_dist > 20.0 {
                    return None;
                }
                let ca_1i = ca_1i.unwrap();
                let ca_i1 = ca_i1.unwrap();
                let ca_1j = ca_1j.unwrap();
                let ca_j1 = ca_j1.unwrap();
                let u1 = ca_i.sub(&ca_1i).normalize();
                let u2 = ca_i1.sub(&ca_i).normalize();
                let u3 = ca_j.sub(&ca_1j).normalize();
                let u4 = ca_j1.sub(&ca_j).normalize();
                let u5 = ca_j.sub(&ca_i).normalize();
                let phi_12 = u1.dot(&u2).acos();
                let phi_34 = u3.dot(&u4).acos();
                let phi_15 = u1.dot(&u5).acos();
                let phi_35 = u3.dot(&u5).acos();
                let phi_14 = u1.dot(&u4).acos();
                let phi_23 = u2.dot(&u3).acos();
                let phi_13 = u1.dot(&u3).acos();
                let seq_dist = j as f32 - i as f32;
                let feature = vec![
                    phi_12, phi_34, phi_15, phi_35, phi_14, phi_23, phi_13, ca_dist, seq_dist
                ];
                Some(feature)
            }
        }
        // append new hash type here
        _ => {
            None
        }
    }
}

pub fn _get_geometric_hash_from_structure(structure: &CompactStructure, hash_type: HashType, nbin_dist: usize, nbin_angle: usize) -> Vec<GeometricHash> {
    let res_bound = get_all_combination(
        structure.num_residues, false
    );

    let mut hash_vec = Vec::with_capacity(res_bound.len());

    res_bound.iter().for_each(|(i, j)| {
        let feature = get_single_feature(*i, *j, structure, hash_type);
        if feature.is_some() {
            let feature = feature.unwrap();
            if nbin_dist == 0 || nbin_angle == 0 {
                let hash = GeometricHash::perfect_hash_default(feature, hash_type);
                hash_vec.push(hash);
            } else {
                let hash = GeometricHash::perfect_hash(feature,hash_type, nbin_dist, nbin_angle);
                hash_vec.push(hash);
            }
        }
    });
    // Reduce memory usage
    hash_vec.shrink_to_fit();
    hash_vec
}

fn get_all_combination(n: usize, include_same: bool) -> Vec<(usize, usize)> {
    let mut res = Vec::new();
    for i in 0..n {
        for j in 0..n {
            if i == j && !include_same {
                continue;
            }
            res.push((i, j));
        }
    }
    res
}

pub fn get_geometric_hash_from_structure(structure: &CompactStructure, hash_type: HashType, nbin_dist: usize, nbin_angle: usize) -> Vec<GeometricHash> {
    let res_bound = CombinationIterator::new(structure.num_residues);
    let mut hash_vec = Vec::new();

    res_bound.for_each(|(i, j)| {
        if i == j {
            return;
        }
        let feature = get_single_feature(i, j, structure, hash_type);
        if let Some(feature) = feature {
            let hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(feature, hash_type)
            } else {
                GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
            };
            hash_vec.push(hash);
        }
    });

    // Reduce memory usage
    hash_vec.shrink_to_fit();
    hash_vec
}

pub fn get_geometric_hash_as_u32_from_structure(structure: &CompactStructure, hash_type: HashType, nbin_dist: usize, nbin_angle: usize) -> Vec<u32> {
    let res_bound = CombinationIterator::new(structure.num_residues);
    let mut hash_vec = Vec::new();

    res_bound.for_each(|(i, j)| {
        if i == j {
            return;
        }
        let feature = get_single_feature(i, j, structure, hash_type);
        if let Some(feature) = feature {
            let hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(feature, hash_type)
            } else {
                GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
            };
            hash_vec.push(hash.as_u32());
        }
    });

    // Reduce memory usage
    hash_vec.shrink_to_fit();
    hash_vec
}


pub fn get_geometric_hash_with_grid(structure: &CompactStructure, id: usize, hash_type: HashType, nbin_dist: usize, nbin_angle: usize, grid_width: f32) -> Vec<(GeometricHash, usize)> {
    let grid_indices = get_grid_index_vector_from_compact(structure, grid_width);
    let res_bound = CombinationIterator::new(structure.num_residues);
    let mut hash_vec = Vec::new();

    res_bound.for_each(|(i, j)| {
        if !nearby(grid_indices[i], grid_indices[j]) {
            return;
        }
        let feature = get_single_feature(i, j, structure, hash_type);
        if let Some(feature) = feature {
            let hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(feature, hash_type)
            } else {
                GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
            };
            let id_grid = merge_id_with_grid(id, grid_indices[i]);
            hash_vec.push((hash, id_grid));
        }
    });

    // Reduce memory usage
    hash_vec.shrink_to_fit();
    hash_vec
}

impl HashType {
    pub fn amino_acid_index(&self) -> Option<Vec<usize>> {
        match self {
            HashType::PDBMotif | HashType::PDBMotifSinCos | HashType::PDBMotifHalf |
            HashType::TrRosetta | HashType::PointPairFeature | HashType::PDBTrRosetta => Some(vec![0, 1]), 
            _ => None
        }
    }

    pub fn dist_index(&self) -> Option<Vec<usize>> {
        match self {
            HashType::PDBMotif | HashType::PDBMotifSinCos | HashType::PDBMotifHalf |
            HashType::PDBTrRosetta  => Some(vec![2, 3]),
            HashType::TrRosetta | HashType::PointPairFeature => Some(vec![2]),
            HashType::TertiaryInteraction => Some(vec![7]),
            _ => None
        }
    }
    
    pub fn angle_index(&self) -> Option<Vec<usize>> {
        match self {
            HashType::PDBMotif | HashType::PDBMotifSinCos | HashType::PDBMotifHalf => Some(vec![4]),
            HashType::TrRosetta => Some(vec![3, 4, 5, 6, 7]),
            HashType::PointPairFeature => Some(vec![3, 4, 5]),
            HashType::PDBTrRosetta => Some(vec![4, 5, 6]), 
            HashType::TertiaryInteraction => Some(vec![0, 1, 2, 3, 4, 5, 6]),
            _ => None
        }
    }
}