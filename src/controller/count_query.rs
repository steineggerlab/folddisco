// Functions for ranking queried results

use rayon::prelude::*;  // Import rayon for parallel iterators

// use std::collections::HashMap;
use rustc_hash::FxHashMap as HashMap;

use crate::index::indextable::FolddiscoIndex;
use crate::prelude::GeometricHash;

use super::io::get_values_with_offset_u16;
use super::map::SimpleHashMap;
use super::result::StructureResult;


// Efficient bit vector for tracking sets of IDs
#[derive(Debug, Clone)]
struct BitVector {
    bits: Vec<u64>,
    capacity: usize,
}

impl BitVector {
    fn new(capacity: usize) -> Self {
        let words_needed = (capacity + 63) / 64;
        Self {
            bits: vec![0u64; words_needed],
            capacity,
        }
    }
    #[inline]
    fn set(&mut self, id: usize) {
        if id < self.capacity {
            let word_idx = id / 64;
            let bit_idx = id % 64;
            self.bits[word_idx] |= 1u64 << bit_idx;
        }
    }
    #[inline]
    fn is_set(&self, id: usize) -> bool {
        if id >= self.capacity {
            return false;
        }
        let word_idx = id / 64;
        let bit_idx = id % 64;
        (self.bits[word_idx] & (1u64 << bit_idx)) != 0
    }
    // #[inline]
    // fn count_ones(&self) -> u32 {
    //     self.bits.iter().map(|&word| word.count_ones()).sum()
    // }
    // Clear all bits
    #[inline]
    fn clear(&mut self) {
        self.bits.fill(0);
    }
}


// Optimized: Use compact data structure for HPC with large pre-allocated vectors
#[derive(Debug, Clone)]
struct CompactEntry {
    node_count: u16,
    edge_count: u32,
    match_count: u32,
    idf_sum: f32,
    // idf_max_per_edge: f32,  // Track max IDF per edge and sum.
    initialized: bool,
}

impl Default for CompactEntry {
    fn default() -> Self {
        Self {
            node_count: 0,
            edge_count: 0,
            match_count: 0,
            idf_sum: 0.0,
            // idf_max_per_edge: 0.0,  // Initialize max IDF per edge
            initialized: false,  // Track if this entry has been initialized
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
    
    let node_grouped = build_node_groups(&sampled, query_map);
    
    // Pre-allocate Vec<CompactEntry> for multi-threading
    let thread_results: Vec<Vec<CompactEntry>> = node_grouped
        .par_iter().map(|(_node, chunk)| {
            let mut local_results: Vec<CompactEntry> = vec![CompactEntry::default(); num_ids];
            let mut node_occupancy = BitVector::new(num_ids);
            let mut edge_occupancy = BitVector::new(num_ids);
            let mut prev_edge = None;
            // let mut idf_max_per_edge = 0.0;  // Track max IDF per edge
            
            for (_i, (e, q)) in chunk.iter().enumerate() {
                // Check if we're starting a new edge
                let need_edge_update = prev_edge.map_or(true, |prev| prev != *e);
                
                if need_edge_update {
                    // Finalize previous edge's counts
                    if prev_edge.is_some() {
                        for nid in 0..num_ids {
                            if edge_occupancy.is_set(nid) {
                                local_results[nid].edge_count += 1;
                                // Update max IDF per edge
                                // local_results[nid].idf_max_per_edge += idf_max_per_edge;
                            }
                        }
                    }
                    // Start tracking new edge
                    edge_occupancy.clear();
                    prev_edge = Some(*e);
                    // idf_max_per_edge = 0.0;  // Reset max IDF for new edge
                }
                
                if let Some((off, hc)) = offset_table.get(q) {
                    if freq_filter.map_or(false, |f| (*hc as f32 / total_hits) > f) {
                        continue;  // Skip queries that do not pass the frequency filter
                    }

                    for &val in get_values_with_offset_u16(value_vec, *off, *hc).iter() {
                        let nid = lookup[val as usize].1;
                        let idf = (total_hits / (*hc as f32)).log2();
                        let entry = &mut local_results[nid];

                        // Initialize entry if not already initialized
                        if !entry.initialized {
                            entry.initialized = true;
                            node_occupancy.set(nid);
                        }
                        entry.match_count += 1;
                        entry.idf_sum += idf;
                        
                        // Track that this edge hits this target
                        edge_occupancy.set(nid);
                        // Update max IDF per edge
                        // if idf > idf_max_per_edge {
                        //     idf_max_per_edge = idf;
                        // }
                    }
                }
            }
            // Set node count based on occupancy
            for nid in 0..num_ids {
                if node_occupancy.is_set(nid) {
                    local_results[nid].node_count = 1;  // Initialize node count
                }
            }
            
            // Finalize the last edge's counts
            for nid in 0..num_ids {
                if edge_occupancy.is_set(nid) {
                    local_results[nid].edge_count += 1;  // Increment edge count for last edge
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
            let mut found_data = false;
            
            // Collect all data from threads for this nid
            for thread_array in &thread_results {
                let entry = &thread_array[nid];
                if entry.initialized {
                    if !found_data {
                        // First thread with data - initialize merged entry
                        merged_entry = entry.clone();
                        found_data = true;
                    } else {
                        // Merge with existing data
                        merged_entry.match_count += entry.match_count;
                        merged_entry.idf_sum += entry.idf_sum;
                        // merged_entry.idf_max_per_edge = entry.idf_max_per_edge;  // Merge max IDF per edge  
                        merged_entry.node_count += entry.node_count;
                        merged_entry.edge_count += entry.edge_count;  // Add edge count merging
                    }
                }
            }
            
            if found_data && merged_entry.match_count > 0 {
                let lookup_entry = &lookup[nid];
                // Count bits in bitmaps for final counts
                merged_entry.idf_sum *= (lookup_entry.2 as f32).powf(-lp);
                // No length penalty applied to max IDF per edge
                let sr = StructureResult::new(
                    &lookup_entry.0,
                    nid,
                    merged_entry.match_count as usize,
                    merged_entry.node_count as usize,
                    merged_entry.edge_count as usize,
                    merged_entry.idf_sum,
                    lookup_entry.2,
                    lookup_entry.3,
                    lookup_entry.4,
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
    
    let node_grouped = build_node_groups(&queries_to_iter, query_map);
    
    // Pre-allocate Vec<CompactEntry> for multi-threading
    let thread_results: Vec<Vec<CompactEntry>> = node_grouped
        .par_iter().map(|(_node, chunk)| {
            let mut local_results: Vec<CompactEntry> = vec![CompactEntry::default(); num_ids];
            let mut node_occupancy = BitVector::new(num_ids);
            let mut edge_occupancy = BitVector::new(num_ids);
            let mut prev_edge = None;

            for (_i, (e, query)) in chunk.iter().enumerate() {
                // Check if we're starting a new edge
                let need_edge_update = prev_edge.map_or(true, |prev| prev != *e);
                
                if need_edge_update {
                    // Finalize previous edge's counts
                    if prev_edge.is_some() {
                        for nid in 0..num_ids {
                            if edge_occupancy.is_set(nid) {
                                local_results[nid].edge_count += 1;
                            }
                        }
                    }
                    // Start tracking new edge
                    edge_occupancy.clear();
                    prev_edge = Some(*e);
                }
                
                let single_queried_values = big_index.get_entries(query.as_u32());
                let hash_count = single_queried_values.len();
                
                if let Some(freq_filter) = freq_filter {
                    if hash_count as f32 / lookup.len() as f32 > freq_filter {
                        continue;  // Skip queries that do not pass the frequency filter
                    }
                }

                let idf = (lookup.len() as f32 / hash_count as f32).log2();

                for &value in single_queried_values.iter() {
                    if value >= lookup.len() {
                        continue;  // Skip invalid values
                    }
                
                    let lookup_entry = &lookup[value];
                    let nid = lookup_entry.1;
                    let entry = &mut local_results[nid];

                    // Initialize entry if not already initialized
                    if !entry.initialized {
                        entry.initialized = true;
                        node_occupancy.set(nid);
                    }
                    entry.match_count += 1;
                    entry.idf_sum += idf;
                    
                    // Track that this edge hits this target
                    edge_occupancy.set(nid);
                }
            }
            
            // Set node count based on occupancy
            for nid in 0..num_ids {
                if node_occupancy.is_set(nid) {
                    local_results[nid].node_count = 1;  // Initialize node count
                }
            }
            
            // Finalize the last edge's counts
            for nid in 0..num_ids {
                if edge_occupancy.is_set(nid) {
                    local_results[nid].edge_count += 1;  // Increment edge count for last edge
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
            let mut found_data = false;
            
            // Collect all data from threads for this nid
            for thread_array in &thread_results {
                let entry = &thread_array[nid];
                if entry.initialized {
                    if !found_data {
                        // First thread with data - initialize merged entry
                        merged_entry = entry.clone();
                        found_data = true;
                    } else {
                        // Merge with existing data
                        merged_entry.match_count += entry.match_count;
                        merged_entry.idf_sum += entry.idf_sum;
                        merged_entry.node_count += entry.node_count;
                        merged_entry.edge_count += entry.edge_count;  // Add edge count merging
                    }
                }
            }
            
            if found_data && merged_entry.match_count > 0 {
                let lookup_entry = &lookup[nid];
                // Apply length penalty to the final IDF score
                merged_entry.idf_sum *= (lookup_entry.2 as f32).powf(-lp);
                
                let sr = StructureResult::new(
                    &lookup_entry.0,
                    nid,
                    merged_entry.match_count as usize,
                    merged_entry.node_count as usize,
                    merged_entry.edge_count as usize,  // Use actual edge count instead of node count
                    merged_entry.idf_sum,
                    lookup_entry.2,
                    lookup_entry.3,
                    lookup_entry.4,
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

// Shared function to build node groups from queries
fn build_node_groups(
    sampled_queries: &[GeometricHash],
    query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
) -> HashMap<usize, Vec<((usize, usize), GeometricHash)>> {
    let mut node_groups: HashMap<usize, Vec<((usize, usize), GeometricHash)>> = HashMap::default();

    for &query in sampled_queries {
        if let Some(((e0, e1), _)) = query_map.get(&query) {
            let edge = (*e0, *e1);
            node_groups.entry(edge.0).or_insert_with(Vec::new).push((edge, query));
        }
    }
    for (_, chunk) in node_groups.iter_mut() {
        // Sort the chunk by edge.
        chunk.sort_by_key(|(edge, _)| *edge);
    }
    node_groups
}