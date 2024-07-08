// Functions for ranking queried results

use std::collections::HashMap;
use std::fmt;
use crate::index::indextable::FolddiscoIndex;
use crate::structure::grid::grid_index_to_tuple;
use crate::{prelude::GeometricHash, structure::grid::convert_to_id_grid_vector};


use super::io::{get_values_with_offset_u16, get_values_with_offset_u24};
use super::map::SimpleHashMap;

#[derive(Clone)]
pub struct QueryResult {
    pub id: String,
    pub nid: usize,
    pub total_match_count: usize,
    pub node_count: usize,
    pub edge_count: usize,
    pub exact_match_count: usize,
    pub overflow_count: usize,
    pub grid_count: usize,
    pub idf: f32,
    pub nres: usize,
    pub plddt: f32,
    pub node_set: HashMap<usize, usize>,
    pub edge_set: HashMap<(usize, usize), usize>,
    pub grid_set: HashMap<(u8, u8, u8), usize>,
    pub pos_set: HashMap<(u16, u16), usize>,
    pub matching_residues: Vec<(String, f32)>,
}

impl QueryResult {
    pub fn new(
        id: String, nid: usize, total_match_count: usize, node_count: usize, edge_count: usize,
        exact_match_count: usize, overflow_count: usize, grid_count: usize,
        idf: f32, nres: usize, plddt: f32
    ) -> Self {
        Self {
            id,
            nid,
            total_match_count,
            node_count,
            edge_count,
            exact_match_count,
            overflow_count,
            grid_count,
            idf,
            nres,
            plddt,
            node_set: HashMap::new(),
            edge_set: HashMap::new(),
            grid_set: HashMap::new(),
            pos_set: HashMap::new(),
            matching_residues: Vec::new(),
        }
    }
}

impl fmt::Display for QueryResult {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let matching_residues_with_score = if self.matching_residues.len() == 0 {
            "NA".to_string()
        } else {
            self.matching_residues.iter().map(
                // Only print score with 4 decimal places
                |(x, y)| format!("{}:{:.4}", x, y)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.exact_match_count, self.overflow_count, self.grid_count,
            self.nres, self.plddt, matching_residues_with_score
            // self.pos_set.len(),
            // self.node_set, self.edge_set, self.grid_set, self.pos_set
        )
    }
}

impl fmt::Debug for QueryResult {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let matching_residues_with_score = if self.matching_residues.len() == 0 {
            "NA".to_string()
        } else {
            self.matching_residues.iter().map(
                // Only print score with 4 decimal places
                |(x, y)| format!("{}:{:.4}", x, y)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.exact_match_count, self.overflow_count, self.grid_count,
            self.nres, self.plddt, matching_residues_with_score
            // self.pos_set.len(),
            // self.node_set, self.edge_set, self.grid_set, self.pos_set
        )
    }
}
// write_fmt
impl QueryResult {
    pub fn write_fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let matching_residues_with_score = if self.matching_residues.len() == 0 {
            "NA".to_string()
        } else {
            self.matching_residues.iter().map(
                // Only print score with 4 decimal places
                |(x, y)| format!("{}:{:.4}", x, y)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.exact_match_count, self.overflow_count, self.grid_count,
            self.nres, self.plddt, matching_residues_with_score
            // self.pos_set.len(),
            // self.node_set, self.edge_set, self.grid_set, self.pos_set
        )
    }
}


pub fn count_query_idmode(
    queries: &Vec<GeometricHash>, query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    offset_table: &SimpleHashMap,
    value_vec: &[u16],
    lookup: &(Vec<String>, Vec<usize>, Vec<usize>, Vec<f32>)
) -> HashMap<usize, QueryResult> {
    println!("Counting queries...");
    let mut query_count_map = HashMap::new();
    for (_i, query) in queries.iter().enumerate() {
        println!("{:?}", query); // WARNING: DELETE AFTER DEBUGGING
        let offset = offset_table.get(query);
        if offset.is_none() {
            continue;
        } 
        let offset = offset.unwrap();
        // Print query and offset
        println!("{:?}, {:?}", query, offset); // WARNING: DELETE AFTER DEBUGGING
        let single_queried_values = get_values_with_offset_u16(value_vec, offset.0, offset.1);
        println!("{:?}", single_queried_values); // WARNING: DELETE AFTER DEBUGGING
        let edge_info = query_map.get(query).unwrap();
        let is_exact = edge_info.1;
        let edge = edge_info.0;
        let hash_count = offset.1;
        for j in 0..single_queried_values.len() {
            let id = lookup.0[single_queried_values[j] as usize].clone();
            let nid = lookup.1[single_queried_values[j] as usize];
            let nres = lookup.2[single_queried_values[j] as usize];
            let plddt = lookup.3[single_queried_values[j] as usize];
            
            let result = query_count_map.get_mut(&nid);
            let idf = (lookup.0.len() as f32 / hash_count as f32).log2();
            let nres_norm = (nres as f32).log2() * -1.0 + 12.0;
            
            if result.is_none() {
                let mut node_set = HashMap::new();
                node_set.insert(edge.0, 1);
                node_set.insert(edge.1, 1);
                let mut edge_set = HashMap::new();
                edge_set.insert(edge, 1);
                let exact_match_count = if is_exact { 1usize } else { 0usize };
                let overflow_count = 0usize;
                let total_match_count = 1usize;
                let mut query_result = QueryResult::new(
                    id, nid, total_match_count, 2, 1, exact_match_count,
                    overflow_count, 1usize, idf + nres_norm, nres, plddt
                );
                query_result.node_set = node_set;
                query_result.edge_set = edge_set;
                query_count_map.insert(nid, query_result);
            } else {
                let result = result.unwrap();
                if result.node_set.contains_key(&edge.0) {
                    let count = result.node_set.get_mut(&edge.0).unwrap();
                    *count += 1;
                } else {
                    result.node_set.insert(edge.0, 1);
                    result.node_count += 1;
                }
                if result.node_set.contains_key(&edge.1) {
                    let count = result.node_set.get_mut(&edge.1).unwrap();
                    *count += 1;
                } else {
                    result.node_set.insert(edge.1, 1);
                    result.node_count += 1;
                }
                let is_overflow = result.edge_set.contains_key(&edge);
                if is_overflow {
                    let count = result.edge_set.get_mut(&edge).unwrap();
                    *count += 1;
                    result.overflow_count += 1;
                } else {
                    result.edge_set.insert(edge, 1);
                    result.edge_count += 1;
                }
                result.total_match_count += 1;
                result.exact_match_count += if is_exact { 1usize } else { 0usize };
                result.idf += idf + nres_norm;
            }
        }
    }
    query_count_map
}

pub fn count_query_bigmode(
    queries: &Vec<GeometricHash>, query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    big_index: &FolddiscoIndex,
    lookup: &(Vec<String>, Vec<usize>, Vec<usize>, Vec<f32>)
) -> HashMap<usize, QueryResult> {
    let mut query_count_map = HashMap::new();
    for (_i, query) in queries.iter().enumerate() {
        
        let single_queried_values = big_index.get_entries(query.as_u32());
        let edge_info = query_map.get(query).unwrap();
        let is_exact = edge_info.1;
        let edge = edge_info.0;
        let hash_count = single_queried_values.len();
        for j in 0..single_queried_values.len() {
            let id = lookup.0[single_queried_values[j]].clone();
            let nid = lookup.1[single_queried_values[j]];
            let nres = lookup.2[single_queried_values[j]];
            let plddt = lookup.3[single_queried_values[j]];
            
            let result = query_count_map.get_mut(&nid);
            let idf = (lookup.0.len() as f32 / hash_count as f32).log2();
            let nres_norm = (nres as f32).log2() * -1.0 + 12.0;
            
            if result.is_none() {
                let mut node_set = HashMap::new();
                node_set.insert(edge.0, 1);
                node_set.insert(edge.1, 1);
                let mut edge_set = HashMap::new();
                edge_set.insert(edge, 1);
                let exact_match_count = if is_exact { 1usize } else { 0usize };
                let overflow_count = 0usize;
                let total_match_count = 1usize;
                let mut query_result = QueryResult::new(
                    id, nid, total_match_count, 2, 1, exact_match_count,
                    overflow_count, 1usize, idf + nres_norm, nres, plddt
                );
                query_result.node_set = node_set;
                query_result.edge_set = edge_set;
                query_count_map.insert(nid, query_result);
            } else {
                let result = result.unwrap();
                if result.node_set.contains_key(&edge.0) {
                    let count = result.node_set.get_mut(&edge.0).unwrap();
                    *count += 1;
                } else {
                    result.node_set.insert(edge.0, 1);
                    result.node_count += 1;
                }
                if result.node_set.contains_key(&edge.1) {
                    let count = result.node_set.get_mut(&edge.1).unwrap();
                    *count += 1;
                } else {
                    result.node_set.insert(edge.1, 1);
                    result.node_count += 1;
                }
                let is_overflow = result.edge_set.contains_key(&edge);
                if is_overflow {
                    let count = result.edge_set.get_mut(&edge).unwrap();
                    *count += 1;
                    result.overflow_count += 1;
                } else {
                    result.edge_set.insert(edge, 1);
                    result.edge_count += 1;
                }
                result.total_match_count += 1;
                result.exact_match_count += if is_exact { 1usize } else { 0usize };
                result.idf += idf + nres_norm;
            }
        }
    }
    query_count_map
}

pub fn count_query_gridmode(
    queries: &Vec<GeometricHash>, query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    // offset_table: &DashMap<GeometricHash, (usize, usize)>,
    offset_table: &SimpleHashMap,
    value_vec: &[u8],
    lookup: &(Vec<String>, Vec<usize>, Vec<usize>, Vec<f32>)
) -> HashMap<usize, QueryResult> {
    let mut query_count_map = HashMap::new();
    for (_i, query) in queries.iter().enumerate() {
        let offset = offset_table.get(query);
        if offset.is_none() {
            continue;
        } 
        let offset = offset.unwrap();
        let single_queried_values = get_values_with_offset_u24(value_vec, offset.0, offset.1);
        let single_queried_values = convert_to_id_grid_vector(single_queried_values);
        let edge_info = query_map.get(query).unwrap();
        let is_exact = edge_info.1;
        let edge = edge_info.0;
        let hash_count = offset.1;
        for j in 0..single_queried_values.len() {
            let id = lookup.0[single_queried_values[j].0 as usize].clone();
            let nid = lookup.1[single_queried_values[j].0 as usize];
            let nres = lookup.2[single_queried_values[j].0 as usize];
            let plddt = lookup.3[single_queried_values[j].0 as usize];
            let grid_index = single_queried_values[j].1;
            
            let result = query_count_map.get_mut(&nid);
            let idf = (lookup.0.len() as f32 / hash_count as f32).log2();
            let nres_norm = (nres as f32).log2() * -1.0 + 12.0;
            
            if result.is_none() {
                let mut node_set = HashMap::new();
                node_set.insert(edge.0, 1);
                node_set.insert(edge.1, 1);
                let mut edge_set = HashMap::new();
                edge_set.insert(edge, 1);
                let mut grid_set = HashMap::new();
                grid_set.insert(grid_index_to_tuple(grid_index), 1);
                let exact_match_count = if is_exact { 1usize } else { 0usize };
                let overflow_count = 0usize;
                let total_match_count = 1usize;
                let mut query_result = QueryResult::new(
                    id, nid, total_match_count, 2, 1, exact_match_count,
                    overflow_count, 1usize, idf + nres_norm, nres, plddt
                );
                query_result.node_set = node_set;
                query_result.edge_set = edge_set;
                query_result.grid_set = grid_set;
                query_count_map.insert(nid, query_result);
            } else {
                let result = result.unwrap();
                if result.node_set.contains_key(&edge.0) {
                    let count = result.node_set.get_mut(&edge.0).unwrap();
                    *count += 1;
                } else {
                    result.node_set.insert(edge.0, 1);
                    result.node_count += 1;
                }
                if result.node_set.contains_key(&edge.1) {
                    let count = result.node_set.get_mut(&edge.1).unwrap();
                    *count += 1;
                } else {
                    result.node_set.insert(edge.1, 1);
                    result.node_count += 1;
                }
                let is_overflow = result.edge_set.contains_key(&edge);
                if is_overflow {
                    let count = result.edge_set.get_mut(&edge).unwrap();
                    *count += 1;
                    result.overflow_count += 1;
                } else {
                    result.edge_set.insert(edge, 1);
                    result.edge_count += 1;
                }
                result.total_match_count += 1;
                result.exact_match_count += if is_exact { 1usize } else { 0usize };
                result.idf += idf + nres_norm;
                if result.grid_set.contains_key(&grid_index_to_tuple(grid_index)) {
                    let count = result.grid_set.get_mut(&grid_index_to_tuple(grid_index)).unwrap();
                    *count += 1;
                } else {
                    result.grid_set.insert(grid_index_to_tuple(grid_index), 1);
                    result.grid_count += 1;
                }
            }
        }
    }
    query_count_map
}





// pub fn concat_hashmap<K, V>(mut map1: HashMap<K, V>, map2: HashMap<K, V>) -> HashMap<K, V> 
// where K: Eq + std::hash::Hash + Clone, V: Clone {
//     for (k, v) in map2 {
//         map1.insert(k, v);
//     }
//     map1
// }

// match_cutoff
// 0,0: no cutoff; 0.0 < cutoff < 1.0: relative cutoff; 1.0 < cutoff: count cutoff

// fn rank_and_filter_result(match_cutoff: String, pdb_query: Vec<GeometricHash>, query_count_vec: Vec<(&usize, &(usize, f32, usize, f32, HashSet<usize>, HashSet<(usize, usize)>, usize, usize, HashSet<(u8, u8, u8)>))>, retrieve: bool, lookup: (Vec<String>, Vec<usize>, Vec<usize>, Vec<f32>), hash_type: HashType, num_bin_dist: usize, num_bin_angle: usize, query_residues: Vec<(u8, u64)>, score_cutoff: f32) {
//     let count_cut = match_cutoff.round() as usize;
//     // let count_cut = 0;
//     let hash_set: HashSet<GeometricHash> = pdb_query.iter().cloned().collect();
//     for (nid, count) in query_count_vec {
//         if count.0 >= count_cut {
//             if retrieve {
//                 let aa_filter = hash_vec_to_aa_pairs(&pdb_query);
//                 let retrieved = retrieve_residue_with_hash(
//                     &hash_set, &aa_filter, &lookup.0[*nid], hash_type, num_bin_dist, num_bin_angle
//                 );
//                 if retrieved.is_some() {
//                     let retrieved = retrieved.unwrap();
//                     let connected = connected(&retrieved, query_residues.len());
//                     let total_matches = retrieved.len();
//                     println!(
//                         "{};uniq_matches={};idf={};total_matches={};connected={};{}",
//                         lookup.0[*nid], count.0, count.1, total_matches,
//                         connected, res_vec_as_string(&retrieved)
//                     );
//                 }
//             } else {
//                 if count.1 < score_cutoff {
//                     continue;
//                 }
//                 println!(
//                     "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
//                     lookup.0[*nid], count.0, count.1, count.2, count.3,
//                     count.4.len(), count.5.len(), count.6, count.7, count.8.len(), count.8.iter().map(|&x| format!("{}{}{}", x.0, x.1, x.2)).collect::<Vec<String>>().join(";")
//                 );
//             }
//         }
//     }
// }
