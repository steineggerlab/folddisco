// Functions for ranking queried results

use dashmap::DashMap;
use rayon::prelude::*;  // Import rayon for parallel iterators

use std::collections::{HashMap, HashSet};
use std::fmt;
use crate::index::indextable::FolddiscoIndex;
use crate::prelude::GeometricHash;

use super::io::get_values_with_offset_u16;
use super::map::SimpleHashMap;
use super::mode::{parse_path_into_given_id_type, IdType};

pub struct QueryResult<'a> {
    pub id: &'a str,
    pub parsed_id: &'a str,
    pub nid: usize,
    pub total_match_count: usize,
    pub node_count: usize,
    pub edge_count: usize,
    pub exact_match_count: usize,
    pub overflow_count: usize,
    pub idf: f32,
    pub nres: usize,
    pub plddt: f32,
    pub node_set: HashSet<usize>,
    pub edge_set: HashSet<(usize, usize)>,
    pub matching_residues: Vec<(String, f32)>,
    pub matching_residues_processed: Vec<(String, f32)>,
}

impl<'a> QueryResult<'a> {
    pub fn new(
        id: &'a str, parsed_id: &'a str, nid: usize, total_match_count: usize, node_count: usize, edge_count: usize,
        exact_match_count: usize, overflow_count: usize, idf: f32, nres: usize, plddt: f32,
        edge: &(usize, usize), 
    ) -> Self {
        let mut node_set = HashSet::new();
        node_set.insert(edge.0);
        node_set.insert(edge.1);
        let mut edge_set = HashSet::new();
        edge_set.insert(*edge);
        Self {
            id,
            parsed_id,
            nid,
            total_match_count,
            node_count,
            edge_count,
            exact_match_count,
            overflow_count,
            idf,
            nres,
            plddt,
            node_set: node_set,
            edge_set: edge_set,
            matching_residues: Vec::new(),
            matching_residues_processed: Vec::new(),
        }
    }
}

impl<'a> fmt::Display for QueryResult<'a> {
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

impl<'a> fmt::Debug for QueryResult<'a> {
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
impl<'a> QueryResult<'a> {
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

pub fn count_query_idmode<'a>(
    queries: &Vec<GeometricHash>, query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    offset_table: &SimpleHashMap,
    value_vec: &[u16],
    lookup: &'a Vec<(String, usize, usize, f32)>,
    id_type: &IdType,
) -> DashMap<usize, QueryResult<'a>> {
    let query_count_map = DashMap::new();  // Use DashMap instead of HashMap

    queries.par_iter().for_each(|query| {  // Use parallel iterator
        if let Some(offset) = offset_table.get(query) {
            let single_queried_values = get_values_with_offset_u16(value_vec, offset.0, offset.1);
            let edge_info = query_map.get(query).unwrap();
            let is_exact = edge_info.1;
            let edge = edge_info.0;
            let hash_count = offset.1;

            for &value in single_queried_values.iter() {
                let id = &lookup[value as usize].0;
                let parsed_id = parse_path_into_given_id_type(id, id_type);
                let nid = lookup[value as usize].1;
                let nres = lookup[value as usize].2;
                let plddt = lookup[value as usize].3;

                let idf = (lookup.len() as f32 / hash_count as f32).log2();
                let nres_norm = (nres as f32).log2() * -1.0 + 12.0;
                let mut is_new: bool = false;
                let entry = query_count_map.entry(nid);
                // Not consuming the entry, so we can modify it
                let mut ref_mut = entry.or_insert_with(|| {
                    let exact_match_count = if is_exact { 1usize } else { 0usize };
                    let overflow_count = 0usize;
                    let total_match_count = 1usize;
                    is_new = true;
                    QueryResult::new(
                        id, parsed_id,nid, total_match_count, 2, 1, exact_match_count,
                        overflow_count, idf + nres_norm, nres, plddt, &edge
                    )
                });
                
                // Modify with ref_mut
                let result = ref_mut.value_mut();
                if !is_new {
                    // Now modify the `result` directly
                    // Check if node_set has edge.0 and edge.1
                    result.node_set.insert(edge.0);
                    result.node_set.insert(edge.1);
                    result.node_count = result.node_set.len();
                    let has_edge = result.edge_set.contains(&edge);
                    if !has_edge {
                        result.edge_set.insert(edge);
                        result.edge_count += 1;
                    } else {
                        result.overflow_count += 1;
                    }
                    result.total_match_count += 1;
                    result.exact_match_count += if is_exact { 1usize } else { 0usize };
                    result.idf += idf + nres_norm;
                }
            }
        }
    });

    query_count_map
}

pub fn count_query_bigmode<'a>(
    queries: &Vec<GeometricHash>, query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    big_index: &FolddiscoIndex,
    lookup: &'a Vec<(String, usize, usize, f32)>,
    idtype: &IdType,
) -> DashMap<usize, QueryResult<'a>> {
    let query_count_map = DashMap::new();  // Use DashMap instead of HashMap

    queries.par_iter().for_each(|query| {  // Use parallel iterator
        let single_queried_values = big_index.get_entries(query.as_u32());
        let edge_info = query_map.get(query).unwrap();
        let is_exact = edge_info.1;
        let edge = edge_info.0;
        let hash_count = single_queried_values.len();

        for &value in single_queried_values.iter() {
            if value >= lookup.len() {
                println!("Error: {} >= {}", value, lookup.len());
                println!("Error query: {:?}", query);
                continue;
            }
            let id = &lookup[value].0;
            let parsed_id = parse_path_into_given_id_type(id, idtype);
            let nid = lookup[value].1;
            let nres = lookup[value].2;
            let plddt = lookup[value].3;

            let idf = (lookup.len() as f32 / hash_count as f32).log2();
            let nres_norm = (nres as f32).log2() * -1.0 + 12.0;
            let mut is_new: bool = false;
            let entry = query_count_map.entry(nid);
            // Not consuming the entry, so we can modify it
            let mut ref_mut = entry.or_insert_with(|| {
                let exact_match_count = if is_exact { 1usize } else { 0usize };
                let overflow_count = 0usize;
                let total_match_count = 1usize;
                is_new = true;
                QueryResult::new(
                    id, parsed_id, nid, total_match_count, 2, 1, exact_match_count,
                    overflow_count, idf + nres_norm, nres, plddt, &edge
                )
            });
                
            // Modify with ref_mut
            let result = ref_mut.value_mut();
            if !is_new {
                // Now modify the `result` directly
                // Check if node_set has edge.0 and edge.1
                result.node_set.insert(edge.0);
                result.node_set.insert(edge.1);
                result.node_count = result.node_set.len();
                let has_edge = result.edge_set.contains(&edge);
                if !has_edge {
                    result.edge_set.insert(edge);
                    result.edge_count += 1;
                } else {
                    result.overflow_count += 1;
                }
                result.total_match_count += 1;
                result.exact_match_count += if is_exact { 1usize } else { 0usize };
                result.idf += idf + nres_norm;
            }
        }
    });
    query_count_map
}
