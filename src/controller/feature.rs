// File: feature.rs
// Created: 2024-03-29 19:05:22
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use std::hash::Hash;

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
        HashType::PDBTrRosetta | HashType::FolddiscoAngle | HashType::FolddiscoDist => {
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

pub fn get_geometric_hash_as_u32_from_structure_with_shifts(
    structure: &CompactStructure, hash_type: HashType, 
    dist_cutoff: f32,
) -> Vec<u32> {
    let res_bound = CombinationIterator::new(structure.num_residues);
    let mut hash_vec = Vec::new();
    let mut feature = vec![0.0; 9];
    
    res_bound.for_each(|(i, j)| {
        if i == j {
            return;
        }
        let has_feature = get_single_feature(i, j, structure, hash_type, dist_cutoff, &mut feature);
        if has_feature {
            let (unique_count, unique_hashes) = GeometricHash::perfect_hash_with_shifts_dedup_inline(&feature, hash_type);
            for k in 0..unique_count as usize {
                hash_vec.push(unique_hashes[k]);
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
            HashType::TrRosetta | HashType::PointPairFeature | HashType::PDBTrRosetta |
            HashType::FolddiscoAngle | HashType::FolddiscoDist => Some(vec![0, 1]), 
            _ => None
        }
    }

    pub fn dist_index(&self) -> Option<Vec<usize>> {
        match self {
            HashType::PDBMotif | HashType::PDBMotifSinCos | HashType::PDBTrRosetta | 
            HashType::Hybrid | HashType::FolddiscoAngle | HashType::FolddiscoDist => Some(vec![2, 3]),
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
            HashType::PDBTrRosetta | HashType::FolddiscoAngle | HashType::FolddiscoDist => Some(vec![4, 5, 6]), 
            HashType::TertiaryInteraction => Some(vec![0, 1, 2, 3, 4, 5, 6]),
            HashType::Hybrid => Some(vec![4, 5, 6, 7, 8]),
            _ => None
        }
    }
    
    pub fn dist_bins(&self, nbin_dist: usize) -> Option<usize> {
        let nbin_dist = if nbin_dist == 0 {
            self.default_dist_bin()
        } else {
            nbin_dist
        };
        match self.dist_index() {
            Some(_) => Some(nbin_dist.pow(self.dist_index().unwrap().len() as u32)),
            None => None,
        }
    }
    
    pub fn angle_bins(&self, nbin_angle: usize) -> Option<usize> {
        let nbin_angle = if nbin_angle == 0 {
            self.default_angle_bin()
        } else {
            nbin_angle
        };
        match self {
            HashType::PDBMotif| HashType::FolddiscoAngle => {
                // Angles are encoded with radian values.
                Some(nbin_angle.pow(self.angle_index().unwrap().len() as u32))
            }
            HashType::FolddiscoDist => {
                if nbin_angle == 16 {
                    // hard code: 8 * 16 * 16 = 2048
                    Some(2048)
                } else {
                    Some(nbin_angle.pow(self.angle_index().unwrap().len() as u32))
                }
            }
            HashType::PDBMotifSinCos | HashType::TrRosetta | HashType::PointPairFeature |
            HashType::PDBTrRosetta | HashType::TertiaryInteraction | HashType::Hybrid => {
                // Angles are encoded with sin-cos values.
                Some((nbin_angle * nbin_angle).pow(self.angle_index().unwrap().len() as u32))
            }
            _ => None
        }
    }
    
    pub fn total_bins(&self, nbin_dist: usize, nbin_angle: usize) -> usize {
        let angle_bins = match self.angle_bins(nbin_angle) {
            Some(num_bins) => num_bins,
            None => 1,
        };
        let dist_bins = match self.dist_bins(nbin_dist) {
            Some(num_bins) => num_bins,
            None => 1,
        };
        let aa_bins = match self.amino_acid_index() {
            Some(_) => 20 * 20,
            None => 1,
        };
        dist_bins * angle_bins * aa_bins
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::controller::io::read_structure_from_path;
    use std::time::Instant;

    #[test]
    fn test_get_geometric_hash_with_shifts() {
        // Load a test structure
        let structure_path = "data/AF-P17538-F1-model_v4.pdb";
        if let Some(structure) = read_structure_from_path(structure_path) {
            let compact_structure = CompactStructure::build(&structure);
            let hash_type = HashType::PDBTrRosetta;
            let dist_cutoff = 20.0;
            
            println!("=== Performance Comparison: Normal vs Shifting ===");
            
            // Time the normal hashing approach
            let start_normal = Instant::now();
            let hash_vec_normal = get_geometric_hash_as_u32_from_structure(
                &compact_structure, hash_type, 16, 4, dist_cutoff, &None
            );
            let duration_normal = start_normal.elapsed();
            
            // Time the shifting approach
            let start_shifts = Instant::now();
            let hash_vec_with_shifts = get_geometric_hash_as_u32_from_structure_with_shifts(
                &compact_structure, hash_type, dist_cutoff
            );
            let duration_shifts = start_shifts.elapsed();
            
            // Display timing results
            println!("Normal approach: {:?} ({} hashes)", duration_normal, hash_vec_normal.len());
            println!("Shifting approach: {:?} ({} hashes)", duration_shifts, hash_vec_with_shifts.len());
            println!("Performance ratio (shifts/normal): {:.2}x", 
                     duration_shifts.as_secs_f64() / duration_normal.as_secs_f64());
            
            // Should have more hashes with shifting than without
            assert!(hash_vec_with_shifts.len() >= hash_vec_normal.len());
            
            // Verify we have unique values due to deduplication
            let mut unique_hashes = std::collections::HashSet::new();
            for hash in &hash_vec_with_shifts {
                unique_hashes.insert(*hash);
            }
            
            // Should have fewer unique hashes than total due to deduplication
            assert!(unique_hashes.len() <= hash_vec_with_shifts.len());
            
            // Display deduplication statistics
            let hash_multiplier = hash_vec_with_shifts.len() as f64 / hash_vec_normal.len() as f64;
            let dedup_efficiency = 100.0 * unique_hashes.len() as f64 / hash_vec_with_shifts.len() as f64;
            
            println!("Hash multiplication factor: {:.1}x", hash_multiplier);
            println!("Total shift hashes: {}", hash_vec_with_shifts.len());
            println!("Unique shift hashes: {}", unique_hashes.len());
            println!("Deduplication efficiency: {:.1}%", dedup_efficiency);
            
            // Calculate performance per hash
            let time_per_hash_normal = duration_normal.as_nanos() as f64 / hash_vec_normal.len() as f64;
            let time_per_hash_shifts = duration_shifts.as_nanos() as f64 / hash_vec_with_shifts.len() as f64;
            
            println!("Time per hash (normal): {:.2} ns", time_per_hash_normal);
            println!("Time per hash (shifts): {:.2} ns", time_per_hash_shifts);
            println!("Per-hash overhead: {:.2}x", time_per_hash_shifts / time_per_hash_normal);
        } else {
            println!("Test structure not found, skipping test");
        }
    }
}

