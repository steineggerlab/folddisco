// Functions for ranking queried results

use dashmap::DashMap;
use rayon::prelude::*;  // Import rayon for parallel iterators

use std::collections::HashMap;

use crate::index::indextable::FolddiscoIndex;
use crate::prelude::GeometricHash;

use super::io::get_values_with_offset_u16;
use super::map::SimpleHashMap;
use super::result::StructureResult;

pub fn count_query_idmode<'a>(
    queries: &Vec<GeometricHash>, query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    offset_table: &SimpleHashMap, value_vec: &[u16], lookup: &'a Vec<(String, usize, usize, f32)>, 
    sampling_ratio: Option<f32>, sampling_count: Option<usize>,
    freq_filter: Option<f32>, length_penalty_power: Option<f32>,
) -> DashMap<usize, StructureResult<'a>> {
    let query_count_map = DashMap::new();  // Use DashMap instead of HashMap
    // Sampling query
    let queries_to_iter = sample_query_idmode(queries, offset_table, sampling_ratio, sampling_count);
    queries_to_iter.par_iter().for_each(|query| {  // Use parallel iterator
        if let Some(offset) = offset_table.get(query) {
            let single_queried_values = get_values_with_offset_u16(value_vec, offset.0, offset.1);
            let hash_count = offset.1;            
            if let Some(freq_filter) = freq_filter {
                if hash_count as f32 / lookup.len() as f32 > freq_filter {
                    // If the hash count is too low, skip the query
                    return;
                }
            }
            let edge_info = query_map.get(query).unwrap();
            let edge = edge_info.0;

            for &value in single_queried_values.iter() {
                let id = &lookup[value as usize].0;
                let nid = lookup[value as usize].1;
                let nres = lookup[value as usize].2;
                let plddt = lookup[value as usize].3;

                let idf = (lookup.len() as f32 / hash_count as f32).log2();

                let mut is_new: bool = false;
                let entry = query_count_map.entry(nid);
                // Not consuming the entry, so we can modify it
                let mut ref_mut = entry.or_insert_with(|| {
                    let total_match_count = 1usize;
                    is_new = true;
                    StructureResult::new(
                        id, nid, total_match_count, 2, 1, 
                        idf, nres, plddt, &edge
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
                    }
                    result.total_match_count += 1;
                    result.idf += idf;

                }
            }
        }
    });
    
    let length_penalty_power = length_penalty_power.unwrap_or(0.5);
    // Normalize idf by nres
    query_count_map.par_iter_mut().for_each(|mut entry| {
        let result = entry.value_mut();
        let length_penalty = (result.nres as f32).powf(-1.0 * length_penalty_power);
        result.idf *= length_penalty;
    });
    
    query_count_map
}

pub fn count_query_bigmode<'a>(
    queries: &Vec<GeometricHash>, query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    big_index: &FolddiscoIndex, lookup: &'a Vec<(String, usize, usize, f32)>, 
    sampling_ratio: Option<f32>, sampling_count: Option<usize>,
    freq_filter: Option<f32>, length_penalty_power: Option<f32>,
) -> DashMap<usize, StructureResult<'a>> {
    let query_count_map = DashMap::new();  // Use DashMap instead of HashMap
    
    let queries_to_iter = sample_query_bigmode(queries, big_index, sampling_ratio, sampling_count);

    queries_to_iter.par_iter().for_each(|query| {  // Use parallel iterator
        let single_queried_values = big_index.get_entries(query.as_u32());
        let hash_count = single_queried_values.len();
        if let Some(freq_filter) = freq_filter {
            if hash_count as f32 / lookup.len() as f32 > freq_filter {
                // If the hash count is too low, skip the query
                return;
            }
        }        
        let edge_info = query_map.get(query).unwrap();
        let edge = edge_info.0;

        for &value in single_queried_values.iter() {
            if value >= lookup.len() {
                println!("Error: {} >= {}", value, lookup.len());
                println!("Error query: {:?}", query);
                continue;
            }
            let id = &lookup[value].0;
            let nid = lookup[value].1;
            let nres = lookup[value].2;
            let plddt = lookup[value].3;

            let idf = (lookup.len() as f32 / hash_count as f32).log2();
            let mut is_new: bool = false;
            let entry = query_count_map.entry(nid);
            // Not consuming the entry, so we can modify it
            let mut ref_mut = entry.or_insert_with(|| {
                let total_match_count = 1usize;
                is_new = true;
                StructureResult::new(
                    id, nid, total_match_count, 2, 1, 
                    idf, nres, plddt, &edge
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
                }
                result.total_match_count += 1;
                result.idf += idf; 
            }
        }
    });
    
    let length_penalty_power = length_penalty_power.unwrap_or(0.5);
    // Normalize idf by nres
    query_count_map.par_iter_mut().for_each(|mut entry| {
        let result = entry.value_mut();
        let length_penalty = (result.nres as f32).powf(-1.0 * length_penalty_power);
        result.idf *= length_penalty;
    });
    
    query_count_map
}

fn sample_query_idmode(
    queries: &Vec<GeometricHash>, offset_table: &SimpleHashMap,
    sampling_ratio: Option<f32>, sampling_count: Option<usize>,
) -> Vec<GeometricHash> {
    match (sampling_ratio, sampling_count) {
        (None, None) => queries.clone(),
        (Some(sampling_ratio), None) => {
            let mut sampled_queries = queries.par_iter().map(|query| {
                if let Some(offset) = offset_table.get(query) {
                    let hash_count = offset.1;
                    (query, hash_count)
                } else {
                    (query, 0)
                }
            }).filter(|(_, hash_count)| *hash_count > 0).collect::<Vec<_>>();
            sampled_queries.sort_by(|a, b| a.1.cmp(&b.1));
            sampled_queries.truncate((sampling_ratio * sampled_queries.len() as f32).ceil() as usize);
            sampled_queries.into_iter().map(|(query, _)| *query).collect()
        },
        (None, Some(sampling_count)) => {
            let mut sampled_queries = queries.par_iter().map(|query| {
                if let Some(offset) = offset_table.get(query) {
                    let hash_count = offset.1;
                    (query, hash_count)
                } else {
                    (query, 0)
                }
            }).filter(|(_, hash_count)| *hash_count > 0).collect::<Vec<_>>();
            sampled_queries.sort_by(|a, b| a.1.cmp(&b.1));
            sampled_queries.truncate(sampling_count);
            sampled_queries.into_iter().map(|(query, _)| *query).collect()
        },
        (Some(_), Some(_)) => queries.clone(),
    }
}


fn sample_query_bigmode(
    queries: &Vec<GeometricHash>, big_index: &FolddiscoIndex, 
    sampling_ratio: Option<f32>, sampling_count: Option<usize>,
) -> Vec<GeometricHash> {
    match (sampling_ratio, sampling_count) {
        (None, None) => queries.clone(),
        (Some(sampling_ratio), None) => {
            let mut sampled_queries = queries.par_iter().map(|query| {
                let single_queried_values = big_index.get_entries(query.as_u32());
                let hash_count = single_queried_values.len();
                (query, hash_count)
            }).collect::<Vec<_>>();
            sampled_queries.sort_by(|a, b| a.1.cmp(&b.1));
            sampled_queries.truncate((sampling_ratio * sampled_queries.len() as f32).ceil() as usize);
            sampled_queries.into_iter().map(|(query, _)| *query).collect()
        },
        (None, Some(sampling_count)) => {
            let mut sampled_queries = queries.par_iter().map(|query| {
                let single_queried_values = big_index.get_entries(query.as_u32());
                let hash_count = single_queried_values.len();
                (query, hash_count)
            }).collect::<Vec<_>>();
            sampled_queries.sort_by(|a, b| a.1.cmp(&b.1));
            sampled_queries.truncate(sampling_count);
            sampled_queries.into_iter().map(|(query, _)| *query).collect()
        },
        (Some(_), Some(_)) => queries.clone(),
    }
}
