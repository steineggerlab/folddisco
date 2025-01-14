// Functions for ranking queried results

use dashmap::DashMap;
use rayon::prelude::*;  // Import rayon for parallel iterators

use std::collections::HashMap;

use crate::index::indextable::FolddiscoIndex;
use crate::prelude::GeometricHash;

use super::io::get_values_with_offset_u16;
use super::map::SimpleHashMap;
use super::query::subset_query_graph_mst;
use super::result::StructureResult;

pub fn count_query_idmode<'a>(
    queries: &Vec<GeometricHash>, query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    offset_table: &SimpleHashMap,
    value_vec: &[u16],
    lookup: &'a Vec<(String, usize, usize, f32)>
) -> DashMap<usize, StructureResult<'a>> {
    let query_count_map = DashMap::new();  // Use DashMap instead of HashMap

    let mst_filtered_query = subset_query_graph_mst(query_map, None, Some(offset_table));
    mst_filtered_query.par_iter().for_each(|query| {  // Use parallel iterator
    // queries.par_iter().for_each(|query| {  // Use parallel iterator
        if let Some(offset) = offset_table.get(query) {
            let single_queried_values = get_values_with_offset_u16(value_vec, offset.0, offset.1);
            let edge_info = query_map.get(query).unwrap();
            let edge = edge_info.0;
            let hash_count = offset.1;

            for &value in single_queried_values.iter() {
                let id = &lookup[value as usize].0;
                let nid = lookup[value as usize].1;
                let nres = lookup[value as usize].2;
                let plddt = lookup[value as usize].3;

                let idf = (lookup.len() as f32 / hash_count as f32).log2();
                let nres_norm = (nres as f32).log2() * -1.0 + 12.0;
                let mut is_new: bool = false;
                let entry = query_count_map.entry(nid);
                // Not consuming the entry, so we can modify it
                let mut ref_mut = entry.or_insert_with(|| {
                    let total_match_count = 1usize;
                    is_new = true;
                    StructureResult::new(
                        id, nid, total_match_count, 2, 1, idf + nres_norm, nres, plddt, &edge
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
    lookup: &'a Vec<(String, usize, usize, f32)>
) -> DashMap<usize, StructureResult<'a>> {
    let query_count_map = DashMap::new();  // Use DashMap instead of HashMap
    // let mut pruned_queries = queries.par_iter().map(|query| {
    //     let single_queried_values = big_index.get_entries(query.as_u32());
    //     let hash_count = single_queried_values.len();
    //     (query, hash_count)
    // }).collect::<Vec<_>>();
    // pruned_queries.sort_by(|a, b| a.1.cmp(&b.1));
    // pruned_queries.truncate(pruned_queries.len() / 2);
    // pruned_queries.par_iter().for_each(|(query, _)| {  // Use parallel iterator
    let mst_filtered_query = subset_query_graph_mst(query_map, Some(big_index), None);
    mst_filtered_query.par_iter().for_each(|query| {  // Use parallel iterator
    // queries.par_iter().for_each(|query| {  // Use parallel iterator
        let single_queried_values = big_index.get_entries(query.as_u32());
        let edge_info = query_map.get(query).unwrap();
        let edge = edge_info.0;
        let hash_count = single_queried_values.len();

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
            let nres_norm = (nres as f32).log2() * -1.0 + 12.0;
            let mut is_new: bool = false;
            let entry = query_count_map.entry(nid);
            // Not consuming the entry, so we can modify it
            let mut ref_mut = entry.or_insert_with(|| {
                let total_match_count = 1usize;
                is_new = true;
                StructureResult::new(
                    id, nid, total_match_count, 2, 1, 
                    idf + nres_norm, nres, plddt, &edge
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
                result.idf += idf + nres_norm;
            }
        }
    });
    query_count_map
}
