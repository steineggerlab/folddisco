
use crate::geometry::util::map_aa_to_u8;
use crate::structure::core::CompactStructure;
use crate::geometry::core::{GeometricHash, HashType};

pub fn get_single_feature(i: usize, j: usize, structure: &CompactStructure, hash_type: HashType) -> Vec<f32> {
    let res1 = structure.get_res_name(i);
    let res2 = structure.get_res_name(j);
    let res1 = map_aa_to_u8(res1) as f32;
    let res2 = map_aa_to_u8(res2) as f32;
    let feature = match &hash_type {
        HashType::PDBMotif => {
            let ca_dist = structure.get_ca_distance(i, j);
            let cb_dist = structure.get_cb_distance(i, j);
            let ca_cb_angle = structure.get_ca_cb_angle(i, j); // degree
            if ca_dist.is_some() && cb_dist.is_some() && ca_cb_angle.is_some() {
                let feature = vec![
                    res1, res2, ca_dist.unwrap(), cb_dist.unwrap(), ca_cb_angle.unwrap()
                ];
                feature
            } else {
                vec![0.0; 5]
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
                feature
            } else {
                vec![0.0; 5]
            }
        },
        HashType::TrRosetta => {
            let feature = structure.get_trrosetta_feature(i, j);
            if feature.is_some() {
                feature.unwrap()
            } else {
                vec![0.0; 6]
            }
        },
        HashType::FoldDiscoDefault => {
            let feature = structure.get_default_feature(i, j);
            if feature.is_some() {
                // Concatenate res1 and res2 to the feature
                let mut feature = feature.unwrap();
                feature.insert(0, res1);
                feature.insert(1, res2);
                feature
            } else {
                vec![0.0; 8]
            }
        },
        _ => {
            todo!("Implement feature-collection methods for other hash types here");
        }
    };
    feature
}




pub fn get_geometric_hash_from_structure(structure: &CompactStructure, hash_type: HashType) -> Vec<GeometricHash> {
    let res_bound = get_all_combination(
        structure.num_residues, false
    );
    let mut hash_vec = Vec::new();

    res_bound.iter().for_each(|(i, j)| {
        let feature = get_single_feature(*i, *j, structure, hash_type);
        let hash = GeometricHash::perfect_hash(feature, hash_type);
        hash_vec.push(hash);
    });
    
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

