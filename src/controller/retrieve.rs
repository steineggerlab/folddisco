// 
use std::{collections::{HashMap, HashSet}, fs::File};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{geometry::util::map_u8_to_aa, prelude::*, structure::core::CompactStructure, utils::combination::CombinationIterator};
use crate::controller::{graph::{connected_components_with_given_node_count, create_index_graph}, query::{make_query, parse_query_string}};

use super::feature::{get_single_feature};

pub fn hash_vec_to_aa_pairs(hash_vec: &Vec<GeometricHash>) -> HashSet<(u32, u32)> {
    let mut output: HashSet<(u32, u32)> = HashSet::new();
    for hash in hash_vec {
        let feature = hash.reverse_hash_default(); 
        output.insert((feature[0] as u32, feature[1] as u32));
    }
    output
}

pub fn retrieve_residue_with_hash(
    hash_set: &HashSet<GeometricHash>, _aa_filter: &HashSet<(u32, u32)>, path: &str, hash_type: HashType, nbin_dist: usize, nbin_angle: usize
) -> Option<Vec<((u8, u8), (u64, u64))>> {

    let file = File::open(path).expect("File not found");
    let pdb = PDBReader::new(file);
    let compact = pdb.read_structure().expect("Error reading structure");
    let compact = compact.to_compact();
    let comb = CombinationIterator::new(compact.num_residues);
    let mut output: Vec<((u8, u8), (u64, u64))> = Vec::new();
    comb.for_each(|(i, j)| {
        // let aa_pair = (map_aa_to_u8(&compact.residue_name[i]) as u32, map_aa_to_u8(&compact.residue_name[j]) as u32);
        // if !aa_filter.contains(&aa_pair) {
        //     return;
        // }
        let feature = get_single_feature(i, j, &compact, hash_type);
        if feature.is_some() {
            let feature = feature.unwrap();
            let curr_hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(feature, hash_type)
            } else {
                GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
            };
            if hash_set.get(&curr_hash).is_some() {
                output.push((
                    (compact.chain_per_residue[i], compact.chain_per_residue[j]),
                    (compact.residue_serial[i], compact.residue_serial[j])
                ));
            }
        }
    });

    if output.is_empty() {
        None
    } else {
        Some(output)
    }
}

pub fn connected(res_ind_vec: &Vec<((u8, u8), (u64, u64))>, _len: usize) -> usize {
    let mut res_set: HashMap<String, usize> = HashMap::new();
    for (i, j) in res_ind_vec {
        let key1 = format!("{}{}", i.0 as char, j.0);
        let key2 = format!("{}{}", i.1 as char, j.1);
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

pub fn cycle(_res_ind_vec: &Vec<((u8, u8), (u64, u64))>) -> usize {
    // Use res_ind_vec as a directed graph and count cycles
    todo!();
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


pub fn prefilter_aa_pair(compact: &CompactStructure, hash: &GeometricHash) -> Vec<(usize, usize)> {
    let mut output: Vec<(usize, usize)> = Vec::new();
    let comb = CombinationIterator::new(compact.num_residues);
    let aa_index = hash.hash_type().amino_acid_index();
    if aa_index.is_none() {
        return output;
    }
    let aa_index = aa_index.unwrap();
    let feature = hash.reverse_hash_default();
    let aa1 = map_u8_to_aa(feature[aa_index[0]] as u8).as_bytes();
    let aa2 = map_u8_to_aa(feature[aa_index[1]] as u8).as_bytes();
    comb.for_each(|(i, j)| {
        if compact.residue_name[i] == aa1 && compact.residue_name[j] == aa2 {
            output.push((i, j));
        }
    });
    output
}

pub fn retrieve_with_prefilter(
    compact: &CompactStructure, hash: &GeometricHash, prefilter: &Vec<(usize, usize)>, nbin_dist: usize, nbin_angle: usize
) -> Vec<(usize, usize)> {
    let mut output: Vec<(usize, usize)> = Vec::new();
    for (i, j) in prefilter {
        let feature = get_single_feature(*i, *j, compact, hash.hash_type());
        if feature.is_some() {
            let feature = feature.unwrap();
            let curr_hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(feature, hash.hash_type())
            } else {
                GeometricHash::perfect_hash(feature, hash.hash_type(), nbin_dist, nbin_angle)
            };
            if curr_hash == *hash {
                output.push((*i, *j));
            }
        }
    }
    output
}

pub fn get_chain_and_res_ind(compact: &CompactStructure, i: usize) -> (u8, u64) {
    (compact.chain_per_residue[i], compact.residue_serial[i])
}
pub fn res_index_to_char(chain: u8, res_ind: u64) -> String {
    format!("{}{}", chain as char, res_ind)
}

pub fn retrieval_wrapper(
    path: &str, node_count: usize, query_vector: &Vec<GeometricHash>, _hash_type: HashType, _nbin_dist: usize, _nbin_angle: usize
) -> Vec<String> {
    let pdb_loaded = PDBReader::new(File::open(&path).expect("File not found"));
    let compact = pdb_loaded.read_structure().expect("Error reading structure");
    let compact = compact.to_compact();
    let mut indices_found: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut output: Vec<String> = Vec::new();
    query_vector.iter().for_each(|hash| {
        let prefiltered = prefilter_aa_pair(&compact, hash);
        let retrieved = retrieve_with_prefilter(&compact, hash, &prefiltered, _nbin_dist, _nbin_angle);
        indices_found.push(retrieved);
    });
    let graph = create_index_graph(&indices_found, &query_vector);
    let connected = connected_components_with_given_node_count(&graph, node_count);
    connected.iter().for_each(|component| {
        let mut res_vec: Vec<String> = Vec::new();
        component.iter().for_each(|&node| {
            let (chain, res_ind) = get_chain_and_res_ind(&compact, node);
            res_vec.push(res_index_to_char(chain, res_ind));
        });
        output.push(res_vec.join(","));
    });
    output
}

pub fn rmsd_for_matched(
    compact1: &CompactStructure, compact2: &CompactStructure, 
    index1: &Vec<usize>, index2: &Vec<usize>
) -> f32 {
    todo!()
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prefilter_aa_pair() {
        let path = String::from("data/serine_peptidases_filtered/4cha.pdb");
        let query_string = "B57,B102,C195";
        let query_residues = parse_query_string(query_string, b'A');
        let hash_type = HashType::PDBTrRosetta;
        let nbin_dist = 16;
        let nbin_angle = 4;
        let exact_match = false;
        let dist_thresholds: Vec<f32> = vec![0.5];
        let angle_thresholds: Vec<f32> = vec![5.0];
        let queries = make_query(
            &path, &query_residues, hash_type, nbin_dist, nbin_angle, exact_match, dist_thresholds, angle_thresholds
        );    
        let pdb_loaded = PDBReader::new(File::open(&path).expect("File not found"));
        let compact = pdb_loaded.read_structure().expect("Error reading structure");
        let compact = compact.to_compact();
        // map
        let mut indices_found: Vec<Vec<(usize, usize)>> = Vec::new();
        measure_time!(queries.iter().for_each(|hash| {
            let prefiltered = prefilter_aa_pair(&compact, hash);
            let retrieved = retrieve_with_prefilter(&compact, hash, &prefiltered, nbin_dist, nbin_angle);
            println!("{}: {:?}", hash, retrieved);
            indices_found.push(retrieved);
        }));
        // Convert to graph
        measure_time!({
            let graph = create_index_graph(&indices_found, &queries);
            let connected = connected_components_with_given_node_count(&graph, query_residues.len());
            println!("{:?}", connected);
            let res_vec: Vec<String> = connected.iter().map(|component| {
                let mut res_vec: Vec<String> = Vec::new();
                component.iter().for_each(|&node| {
                    let (chain, res_ind) = get_chain_and_res_ind(&compact, node);
                    res_vec.push(res_index_to_char(chain, res_ind));
                });
                res_vec.join(",")
            }).collect();
            println!("{:?}", res_vec);
        });



    } 
    

    
    
    #[test]
    fn test_retrieve_residue_with_hash() {
        let path = String::from("data/serine_peptidases_filtered/4cha.pdb");
        let query_string = "B57,B102,C195";
        // let path = String::from("analysis/1g91.pdb");
        // let query_string = "A30,A32,A35";
        let query_residues = parse_query_string(query_string, b'A');
        let hash_type = HashType::PDBMotifSinCos;
        let nbin_dist = 16;
        let nbin_angle = 4;
        let exact_match = false;
        let dist_thresholds: Vec<f32> = vec![0.5];
        let angle_thresholds: Vec<f32> = vec![5.0, 10.0];
        let queries = make_query(
            &path, &query_residues, hash_type, nbin_dist, nbin_angle, exact_match, dist_thresholds, angle_thresholds
        );
        let hash_set: HashSet<GeometricHash> = queries.iter().cloned().collect();
    
        let pdb_loaded = PDBReader::new(File::open(&path).expect("File not found"));
        let compact = pdb_loaded.read_structure().expect("Error reading structure");
        let compact = compact.to_compact();
        // print num residues
        println!("Num residues: {}", compact.num_residues);
        let aa_filter = hash_vec_to_aa_pairs(&queries);
        println!("{:?}", aa_filter);
        let retrieved = measure_time!(retrieve_residue_with_hash(
            &hash_set, &aa_filter, &path, hash_type, nbin_dist, nbin_angle
        ));
        println!("{:?}", retrieved);
        let connected = connected(&retrieved.unwrap(), query_residues.len());
        println!("{:?}", connected);
    }
}