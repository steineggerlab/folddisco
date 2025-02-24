// File: feature.rs
// Created: 2024-03-29 19:05:22
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use crate::utils::convert::{map_aa_to_u8, map_aa_to_u8_group};
use crate::structure::core::CompactStructure;
use crate::geometry::core::{GeometricHash, HashType};
use crate::utils::combination::CombinationIterator;

pub fn get_single_feature(
    i: usize, j: usize, structure: &CompactStructure, hash_type: HashType, 
    dist_cutoff: f32, feature_container: &mut Vec<f32>
) -> bool {
    if i == j {
        return false;
    }
    let res1 = structure.get_res_name(i);
    let res2 = structure.get_res_name(j);
    let res1 = map_aa_to_u8(res1) as f32;
    let res2 = map_aa_to_u8(res2) as f32;
    if res1 == 255.0 || res2 == 255.0 {
        return false;
    }
    match &hash_type {
        HashType::PDBMotif => {
            let ca_dist = structure.get_ca_distance(i, j);
            if let Some(ca_dist) = ca_dist {
                if ca_dist > dist_cutoff {
                    return false;
                }
            } 
            let cb_dist = structure.get_cb_distance(i, j);
            // Uses angle in degree
            let ca_cb_angle = structure.get_ca_cb_angle(i, j, false); // degree
            if ca_dist.is_some() && cb_dist.is_some() && ca_cb_angle.is_some() {
                feature_container[0] = res1;
                feature_container[1] = res2;
                feature_container[2] = ca_dist.unwrap();
                feature_container[3] = cb_dist.unwrap();
                feature_container[4] = ca_cb_angle.unwrap();
                return true;
            } else {
                return false;
            }
        },
        HashType::PDBMotifSinCos => {
            let ca_dist = structure.get_ca_distance(i, j);
            if let Some(ca_dist) = ca_dist {
                if ca_dist > dist_cutoff {
                    return false;
                }
            }
            let cb_dist = structure.get_cb_distance(i, j);
            // Uses angle in radians
            let ca_cb_angle_radian = structure.get_ca_cb_angle(i, j, true);
            if ca_dist.is_some() && cb_dist.is_some() && ca_cb_angle_radian.is_some() {
                feature_container[0] = res1;
                feature_container[1] = res2;
                feature_container[2] = ca_dist.unwrap();
                feature_container[3] = cb_dist.unwrap();
                feature_container[4] = ca_cb_angle_radian.unwrap();
                return true;
            } else {
                return false;
            }
        },
        HashType::TrRosetta => { 
            let feature = structure.get_trrosetta_feature(i, j, dist_cutoff);
            if let Some(feature) = feature {
                feature_container[0] = res1;
                feature_container[1] = res2;
                feature_container[2] = feature.0;
                feature_container[3] = feature.1;
                feature_container[4] = feature.2;
                feature_container[5] = feature.3;
                feature_container[6] = feature.4;
                feature_container[7] = feature.5;
                return true;
            } else {
                return false;
            }
        },
        HashType::PDBTrRosetta => {
            let feature = structure.get_pdb_tr_feature(i, j, dist_cutoff);
            if feature.is_some() {
                let feature = feature.unwrap();
                feature_container[0] = res1;
                feature_container[1] = res2;
                feature_container[2] = feature.0;
                feature_container[3] = feature.1;
                feature_container[4] = feature.2;
                feature_container[5] = feature.3;
                feature_container[6] = feature.4;
                return true;
            } else {
                return false;
            }
        },
        HashType::PointPairFeature => {
            let feature = structure.get_ppf(i, j, dist_cutoff);
            if let Some(feature) = feature {
                feature_container[0] = res1;
                feature_container[1] = res2;
                feature_container[2] = feature[0];
                feature_container[3] = feature[1];
                feature_container[4] = feature[2];
                feature_container[5] = feature[3];
                return true;
            } else {
                return false;
            }
        },
        HashType::TertiaryInteraction => {
            if i == 0 || j == 0 || i == structure.num_residues - 1 || j == structure.num_residues - 1 {
                return false;
            }
            let ca_1i = structure.get_ca(i-1);
            let ca_i = structure.get_ca(i);
            let ca_i1 = structure.get_ca(i+1);
            let ca_1j = structure.get_ca(j-1);
            let ca_j = structure.get_ca(j);
            let ca_j1 = structure.get_ca(j+1);
            
            if ca_1i.is_none() || ca_i.is_none() || ca_i1.is_none() || ca_1j.is_none() || ca_j.is_none() || ca_j1.is_none() {
                return false;
            } else {
                let ca_i = ca_i.unwrap();
                let ca_j = ca_j.unwrap();
                let ca_dist = ca_i.distance(&ca_j);
                if ca_dist > dist_cutoff {
                    return false;
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
                // let feature = vec![
                //     phi_12, phi_34, phi_15, phi_35, phi_14, phi_23, phi_13, ca_dist, seq_dist
                // ];
                feature_container[0] = phi_12;
                feature_container[1] = phi_34;
                feature_container[2] = phi_15;
                feature_container[3] = phi_35;
                feature_container[4] = phi_14;
                feature_container[5] = phi_23;
                feature_container[6] = phi_13;
                feature_container[7] = ca_dist;
                feature_container[8] = seq_dist;
                return true;
            }
        },
        HashType::Hybrid => {
            if i == 0 || j == 0 || i == structure.num_residues - 1 || j == structure.num_residues - 1 {
                return false;
            }
            let res1 = structure.get_res_name(i);
            let res2 = structure.get_res_name(j);
            let res1_group = map_aa_to_u8_group(res1) as f32;
            let res2_group = map_aa_to_u8_group(res2) as f32;

            let feature = structure.get_hybrid_feature(i, j, dist_cutoff);
            if let Some(feature) = feature {
                feature_container[0] = res1_group;
                feature_container[1] = res2_group;
                feature_container[2] = feature.0;
                feature_container[3] = feature.1;
                feature_container[4] = feature.2;
                feature_container[5] = feature.3;
                feature_container[6] = feature.4;
                feature_container[7] = feature.5;
                feature_container[8] = feature.6;
                return true;
            } else {
                return false;
            }
        },
        // append new hash type here
        _ => {
            return false;
        }
    }
}

pub fn get_geometric_hash_as_u32_from_structure(
    structure: &CompactStructure, hash_type: HashType, 
    nbin_dist: usize, nbin_angle: usize, dist_cutoff: f32,
    multiple_bins: &Option<Vec<(usize, usize)>>,
) -> Vec<u32> {
    let res_bound = CombinationIterator::new(structure.num_residues);
    let mut hash_vec = Vec::with_capacity(res_bound.len());
    let mut feature = vec![0.0; 9];
    res_bound.for_each(|(i, j)| {
        if i == j {
            return;
        }
        let has_feature = get_single_feature(i, j, structure, hash_type, dist_cutoff, &mut feature);
        if has_feature {
            if let Some(multiple_bins) = multiple_bins {
                for (nbin_dist, nbin_angle) in multiple_bins.iter() {
                    let hash = GeometricHash::perfect_hash_as_u32(&feature, hash_type, *nbin_dist, *nbin_angle);
                    hash_vec.push(hash);
                }
            } else {
                if nbin_dist == 0 || nbin_angle == 0 {
                    let hash = GeometricHash::perfect_hash_default_as_u32(&feature, hash_type);
                    hash_vec.push(hash);
                } else {
                    let hash = GeometricHash::perfect_hash_as_u32(&feature, hash_type, nbin_dist, nbin_angle);
                    hash_vec.push(hash);
                }
            }
        }
    });
    // Reduce memory usage
    hash_vec.shrink_to_fit();
    hash_vec
}

impl HashType {
    pub fn amino_acid_index(&self) -> Option<Vec<usize>> {
        match self {
            HashType::PDBMotif | HashType::PDBMotifSinCos | 
            HashType::TrRosetta | HashType::PointPairFeature | HashType::PDBTrRosetta => Some(vec![0, 1]), 
            _ => None
        }
    }

    pub fn dist_index(&self) -> Option<Vec<usize>> {
        match self {
            HashType::PDBMotif | HashType::PDBMotifSinCos | HashType::PDBTrRosetta | HashType::Hybrid  => Some(vec![2, 3]),
            HashType::TrRosetta | HashType::PointPairFeature => Some(vec![2]),
            HashType::TertiaryInteraction => Some(vec![7]),
            _ => None
        }
    }
    
    pub fn angle_index(&self) -> Option<Vec<usize>> {
        match self {
            HashType::PDBMotif | HashType::PDBMotifSinCos => Some(vec![4]),
            HashType::TrRosetta => Some(vec![3, 4, 5, 6, 7]),
            HashType::PointPairFeature => Some(vec![3, 4, 5]),
            HashType::PDBTrRosetta => Some(vec![4, 5, 6]), 
            HashType::TertiaryInteraction => Some(vec![0, 1, 2, 3, 4, 5, 6]),
            HashType::Hybrid => Some(vec![4, 5, 6, 7, 8]),
            _ => None
        }
    }
}

