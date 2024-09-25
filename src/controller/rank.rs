// Functions for ranking queried results

use std::collections::HashMap;
use std::fmt;
use crate::index::indextable::FolddiscoIndex;
use crate::prelude::GeometricHash;


use super::io::get_values_with_offset_u16;
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
    pub idf: f32,
    pub nres: usize,
    pub plddt: f32,
    pub node_set: HashMap<usize, usize>,
    pub edge_set: HashMap<(usize, usize), usize>,
    pub pos_set: HashMap<(u16, u16), usize>,
    pub matching_residues: Vec<(String, f32)>,
    pub matching_residues_processed: Vec<(String, f32)>,
}

impl QueryResult {
    pub fn new(
        id: String, nid: usize, total_match_count: usize, node_count: usize, edge_count: usize,
        exact_match_count: usize, overflow_count: usize,  idf: f32, nres: usize, plddt: f32
    ) -> Self {
        Self {
            id,
            nid,
            total_match_count,
            node_count,
            edge_count,
            exact_match_count,
            overflow_count,
            idf,
            nres,
            plddt,
            node_set: HashMap::new(),
            edge_set: HashMap::new(),
            pos_set: HashMap::new(),
            matching_residues: Vec::new(),
            matching_residues_processed: Vec::new(),
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
        let matching_residues_processed_with_score = if self.matching_residues_processed.len() == 0 {
            "NA".to_string()
        } else {
            self.matching_residues_processed.iter().map(
                // Only print score with 4 decimal places
                |(x, y)| format!("{}:{:.4}", x, y)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.exact_match_count, self.overflow_count,
            self.nres, self.plddt, matching_residues_with_score, matching_residues_processed_with_score
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
        let matching_residues_processed_with_score = if self.matching_residues_processed.len() == 0 {
            "NA".to_string()
        } else {
            self.matching_residues_processed.iter().map(
                // Only print score with 4 decimal places
                |(x, y)| format!("{}:{:.4}", x, y)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.exact_match_count, self.overflow_count,
            self.nres, self.plddt, matching_residues_with_score, matching_residues_processed_with_score
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
        let matching_residues_processed_with_score = if self.matching_residues_processed.len() == 0 {
            "NA".to_string()
        } else {
            self.matching_residues_processed.iter().map(
                // Only print score with 4 decimal places
                |(x, y)| format!("{}:{:.4}", x, y)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.exact_match_count, self.overflow_count,
            self.nres, self.plddt, matching_residues_with_score, matching_residues_processed_with_score
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
    let mut query_count_map = HashMap::new();
    for (_i, query) in queries.iter().enumerate() {
        let offset = offset_table.get(query);
        if offset.is_none() {
            continue;
        } 
        let offset = offset.unwrap();
        let single_queried_values = get_values_with_offset_u16(value_vec, offset.0, offset.1);
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
                    overflow_count, idf + nres_norm, nres, plddt
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
                    overflow_count, idf + nres_norm, nres, plddt
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

