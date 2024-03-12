// 
use std::{collections::{HashMap, HashSet}, fs::File, hash::Hash};
use crate::{geometry::util::map_aa_to_u8, prelude::*, utils::combination::CombinationIterator};

use super::feature::{self, get_single_feature};

pub fn hash_vec_to_aa_pairs(hash_vec: &Vec<GeometricHash>) -> HashSet<(u32, u32)> {
    let mut output: HashSet<(u32, u32)> = HashSet::new();
    for hash in hash_vec {
        let feature = hash.reverse_hash_default(); 
        output.insert((feature[0] as u32, feature[1] as u32));
    }
    output
}

pub fn retrieve_residue_with_hash(
    hash_set: &HashSet<GeometricHash>, aa_filter: &HashSet<(u32, u32)>, path: &str, hash_type: HashType, nbin_dist: usize, nbin_angle: usize
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

pub fn connected(res_ind_vec: &Vec<((u8, u8), (u64, u64))>, len: usize) -> usize {
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

pub fn cycle(res_ind_vec: &Vec<((u8, u8), (u64, u64))>) -> usize {
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


#[cfg(test)]
mod tests {
    use super::*;
    use crate::controller::query::{parse_query_string, make_query};

    #[test]
    fn test_retrieve_residue_with_hash() {
        let path = String::from("data/serine_peptidases_filtered/4cha.pdb");
        let query_string = "B57,B102,C195";
        // let path = String::from("analysis/1g91.pdb");
        // let query_string = "A30,A32,A35";
        let query_residues = parse_query_string(query_string);
        let hash_type = HashType::PDBMotifSinCos;
        let nbin_dist = 16;
        let nbin_angle = 3;
        let check_nearby = false;
        let queries = make_query(&path, &query_residues, hash_type, nbin_dist, nbin_angle, check_nearby);
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