// Functions for ranking queried results

use rayon::prelude::*;  // Import rayon for parallel iterators

use std::collections::HashMap;

use crate::index::indextable::FolddiscoIndex;
use crate::prelude::GeometricHash;

use super::io::get_values_with_offset_u16;
use super::map::SimpleHashMap;
use super::result::StructureResult;

// Optimized: Use compact data structure for HPC with large pre-allocated vectors
#[derive(Debug, Clone)]
struct CompactEntry {
    db_key: usize,
    node_count: u16,
    match_count: u16,
    idf_sum: f32,
    nres: usize,
    plddt: f32,
    // Use bit flags instead of HashSets for efficiency
    nodes_bitmap: u128,  // Support up to 128 unique nodes
    initialized: bool,
}

impl Default for CompactEntry {
    fn default() -> Self {
        Self {
            db_key: 0,
            node_count: 0,
            match_count: 0,
            idf_sum: 0.0,
            nres: 0,
            plddt: 0.0,
            nodes_bitmap: 0,
            initialized: false,
        }
    }
}

pub fn count_query_idmode<'a>(
    queries: &[GeometricHash],
    query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    offset_table: &SimpleHashMap,
    value_vec: &[u16],
    lookup: &'a [(String, usize, usize, f32, usize)],
    sampling_ratio: Option<f32>,
    sampling_count: Option<usize>,
    freq_filter: Option<f32>,
    length_penalty_power: Option<f32>,
) -> Vec<(usize, StructureResult<'a>)> {
    let sampled = sample_query_idmode(queries, offset_table, sampling_ratio, sampling_count);
    let total_hits = lookup.len() as f32;
    let lp = length_penalty_power.unwrap_or(0.5);
    let num_ids = lookup.len();
    
    // Pre-allocate thread results for zero allocation overhead - use only CompactEntry
    let num_threads = rayon::current_num_threads();
    let chunk_size = (sampled.len() + num_threads - 1) / num_threads;
    
    // Pre-allocate Vec<CompactEntry> for multi-threading
    let thread_results: Vec<Vec<CompactEntry>> = sampled
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local_results: Vec<CompactEntry> = vec![CompactEntry::default(); num_ids];
            
            for q in chunk.iter() {
                if let Some((off, hc)) = offset_table.get(q) {
                    if freq_filter.map_or(false, |f| (*hc as f32 / total_hits) > f) {
                        continue;  // Skip queries that do not pass the frequency filter
                    }

                    let edge = query_map.get(q).unwrap().0;

                    for &val in get_values_with_offset_u16(value_vec, *off, *hc).iter() {
                        let lookup_entry = &lookup[val as usize];
                        let nid = lookup_entry.1;
                        let db_key = lookup_entry.4;
                        let nres = lookup_entry.2;
                        let plddt = lookup_entry.3;
                        let idf = (total_hits / (*hc as f32)).log2();

                        let entry = &mut local_results[nid];

                        // Initialize entry if not already initialized
                        if !entry.initialized {
                            entry.db_key = db_key;
                            entry.nres = nres;
                            entry.plddt = plddt;
                            entry.initialized = true;
                        }

                        // Track unique nodes using bitmaps
                        let node0_bit = 1u128 << (edge.0 % 128);
                        let node1_bit = 1u128 << (edge.1 % 128);
                        
                        entry.nodes_bitmap |= node0_bit | node1_bit;
                        
                        entry.match_count += 1;
                        entry.idf_sum += idf;
                    }
                }
            }
            local_results
        })
        .collect();

    // Parallel processing with bitmap aggregation
    let results: Vec<(usize, StructureResult<'a>)> = (0..num_ids)
        .into_par_iter()
        .filter_map(|nid| {
            let mut merged_entry = CompactEntry::default();
            let mut combined_nodes_bitmap = 0u128;
            let mut found_data = false;
            
            // Collect all data from threads for this nid
            for thread_array in &thread_results {
                let entry = &thread_array[nid];
                if entry.initialized {
                    if !found_data {
                        // First thread with data - initialize merged entry
                        merged_entry = entry.clone();
                        combined_nodes_bitmap = entry.nodes_bitmap;
                        found_data = true;
                    } else {
                        // Merge with existing data
                        merged_entry.match_count += entry.match_count;
                        merged_entry.idf_sum += entry.idf_sum;
                        combined_nodes_bitmap |= entry.nodes_bitmap;
                    }
                }
            }
            
            if found_data && merged_entry.match_count > 0 {
                // Count bits in bitmaps for final counts
                merged_entry.node_count = combined_nodes_bitmap.count_ones() as u16;
                merged_entry.idf_sum *= (merged_entry.nres as f32).powf(-lp);
                
                let lookup_entry = &lookup[nid];
                let dummy_edge = (0, 0);
                let sr = StructureResult::new(
                    &lookup_entry.0,
                    nid,
                    merged_entry.match_count as usize,
                    merged_entry.node_count as usize,
                    merged_entry.idf_sum,
                    merged_entry.nres,
                    merged_entry.plddt,
                    &dummy_edge,
                    merged_entry.db_key,
                );
                Some((nid, sr))
            } else {
                None
            }
        })
        .collect();

    results
}

pub fn count_query_bigmode<'a>(
    queries: &Vec<GeometricHash>, query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    big_index: &FolddiscoIndex, lookup: &'a Vec<(String, usize, usize, f32, usize)>, 
    sampling_ratio: Option<f32>, sampling_count: Option<usize>,
    freq_filter: Option<f32>, length_penalty_power: Option<f32>,
) -> Vec<(usize, StructureResult<'a>)> {
    let queries_to_iter = sample_query_bigmode(queries, big_index, sampling_ratio, sampling_count);
    let num_ids = lookup.len();
    let lp = length_penalty_power.unwrap_or(0.5);
    
    // Pre-allocate thread results for zero allocation overhead - use only CompactEntry
    let num_threads = rayon::current_num_threads();
    let chunk_size = (queries_to_iter.len() + num_threads - 1) / num_threads;
    
    // Pre-allocate Vec<CompactEntry> for multi-threading on HPC
    let thread_results: Vec<Vec<CompactEntry>> = queries_to_iter
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local_results: Vec<CompactEntry> = vec![CompactEntry::default(); num_ids];
            
            for query in chunk.iter() {
                let single_queried_values = big_index.get_entries(query.as_u32());
                let hash_count = single_queried_values.len();
                
                if let Some(freq_filter) = freq_filter {
                    if hash_count as f32 / lookup.len() as f32 > freq_filter {
                        continue;  // Skip queries that do not pass the frequency filter
                    }
                }
                
                let edge = query_map.get(query).unwrap().0;

                for &value in single_queried_values.iter() {
                    if value >= lookup.len() {
                        continue;  // Skip invalid values
                    }
                    
                    let lookup_entry = &lookup[value];
                    let nid = lookup_entry.1;
                    let nres = lookup_entry.2;
                    let plddt = lookup_entry.3;
                    let db_key = lookup_entry.4;
                    let idf = (lookup.len() as f32 / hash_count as f32).log2();

                    let entry = &mut local_results[nid];

                    // Initialize entry if not already initialized
                    if !entry.initialized {
                        entry.nres = nres;
                        entry.plddt = plddt;
                        entry.db_key = db_key;
                        entry.initialized = true;
                    }

                    // Track unique nodes using bitmaps
                    let node0_bit = 1u128 << (edge.0 % 128);
                    let node1_bit = 1u128 << (edge.1 % 128);
                    
                    entry.nodes_bitmap |= node0_bit | node1_bit;
                    
                    entry.match_count += 1;
                    entry.idf_sum += idf;
                }
            }
            local_results
        })
        .collect();

    // Lock-free parallel processing with bitmap aggregation
    let results: Vec<(usize, StructureResult<'a>)> = (0..num_ids)
        .into_par_iter()
        .filter_map(|nid| {
            let mut merged_entry = CompactEntry::default();
            let mut combined_nodes_bitmap = 0u128;
            let mut found_data = false;
            
            // Collect all data from threads for this nid
            for thread_array in &thread_results {
                let entry = &thread_array[nid];
                if entry.initialized {
                    if !found_data {
                        // First thread with data - initialize merged entry
                        merged_entry = entry.clone();
                        combined_nodes_bitmap = entry.nodes_bitmap;
                        found_data = true;
                    } else {
                        // Merge with existing data
                        merged_entry.match_count += entry.match_count;
                        merged_entry.idf_sum += entry.idf_sum;
                        combined_nodes_bitmap |= entry.nodes_bitmap;
                    }
                }
            }
            
            if found_data && merged_entry.match_count > 0 {
                // Count bits in bitmaps for final counts
                merged_entry.node_count = combined_nodes_bitmap.count_ones() as u16;
                merged_entry.idf_sum *= (merged_entry.nres as f32).powf(-lp);
                
                let lookup_entry = &lookup[nid];
                let dummy_edge = (0, 0);
                let sr = StructureResult::new(
                    &lookup_entry.0,
                    nid,
                    merged_entry.match_count as usize,
                    merged_entry.node_count as usize,
                    merged_entry.idf_sum,
                    merged_entry.nres,
                    merged_entry.plddt,
                    &dummy_edge,
                    merged_entry.db_key,
                );
                Some((nid, sr))
            } else {
                None
            }
        })
        .collect();

    results
}

fn sample_query_idmode(
    queries: &[GeometricHash], offset_table: &SimpleHashMap,
    sampling_ratio: Option<f32>, sampling_count: Option<usize>,
) -> Vec<GeometricHash> {
    match (sampling_ratio, sampling_count) {
        (None, None) => queries.to_vec(),
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
        (Some(_), Some(_)) => queries.to_vec(),
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
