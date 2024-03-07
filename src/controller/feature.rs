
use crate::geometry::util::map_aa_to_u8;
use crate::structure::core::CompactStructure;
use crate::geometry::core::{GeometricHash, HashType};

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
                    res1, res2, ca_dist.unwrap(), cb_dist.unwrap(), ca_cb_angle.unwrap()
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
        HashType::PDBMotifSinCos => {
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
            let feature = structure.get_trrosetta_feature(i, j);
            if feature.is_some() {
                let feature = feature.unwrap();
                if feature[2] > 20.0 {
                    None
                } else {
                    Some(feature)
                }
            } else {
                None
            }
        },
        HashType::FoldDiscoDefault | HashType::Default32bit => {
            let feature = structure.get_default_feature(i, j);
            if feature.is_some() {
                // Concatenate res1 and res2 to the feature
                let mut feature = feature.unwrap();
                if feature[2] > 20.0 {
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
                if feature[2] > 20.0 {
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
        // append new hash type here
        _ => {
            None
        }
    }
}

pub fn get_geometric_hash_from_structure(structure: &CompactStructure, hash_type: HashType) -> Vec<GeometricHash> {
    let start = std::time::Instant::now();
    let res_bound = get_all_combination(
        structure.num_residues, false
    );

    let mut hash_vec = Vec::with_capacity(res_bound.len());

    res_bound.iter().for_each(|(i, j)| {
        let feature = get_single_feature(*i, *j, structure, hash_type);
        if feature.is_some() {
            let feature = feature.unwrap();
            let hash = GeometricHash::perfect_hash(feature, hash_type);
            hash_vec.push(hash);
        }
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed > 60.0 {
            eprintln!("Elapsed time: {:.2} seconds", elapsed);
            // Print length of hash_vec
            eprintln!("Length of hash_vec: {}", hash_vec.len());
            eprintln!("Residue boundary: {:?}", res_bound);
            eprintln!("Structure info: {:?}", structure);
            eprintln!("Total {} residues", structure.num_residues);
            eprintln!("Residue: {:?}", structure.residue_serial);
            eprintln!("Residue name: {:?}", structure.residue_name);
            eprintln!("Chain: {:?}", structure.chains);
            eprintln!("Chain per residue: {:?}", structure.chain_per_residue);
            // Terminate
            std::process::exit(1);
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

