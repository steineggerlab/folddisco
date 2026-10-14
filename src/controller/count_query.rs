// First-pass scoring: count index hits per target structure before residue matching.

use rayon::prelude::*;

// use std::collections::HashMap;
use rustc_hash::FxHashMap as HashMap;

use crate::index::indextable::FolddiscoIndex;
use crate::index::lookup::LookupTable;
use crate::prelude::GeometricHash;

use super::query::QueryHashMap;
use super::result::StructureResult;


// Fixed-capacity bit set over structure ids
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
    #[inline]
    fn clear(&mut self) {
        self.bits.fill(0);
    }
}


// Per-target counters, one dense vector per query node
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
            initialized: false,
        }
    }
}

/// Score every structure hit by `queries`.
///
/// Hashes are grouped by the first residue of their query edge (one "node"); per
/// target, `node_count` is the number of nodes with a hit, `edge_count` the number
/// of residue pairs with a hit, and `idf` the summed log2(N / hash frequency)
/// scaled by `nres^-length_penalty` (default 0.5).
/// Substituted hashes add only their weighted IDF, once per edge, and only to targets
/// without an exact hit on that edge.
pub fn count_query<'a>(
    queries: &Vec<GeometricHash>, query_map: &QueryHashMap,
    index: &FolddiscoIndex,
    lookup: &'a LookupTable, 
    sampling_ratio: Option<f32>, sampling_count: Option<usize>,
    freq_filter: Option<f32>, length_penalty_power: Option<f32>,
) -> Vec<(usize, StructureResult<'a>)> {
    let queries_to_iter = sample_query(queries, index, sampling_ratio, sampling_count);
    let num_ids = lookup.len();
    let lp = length_penalty_power.unwrap_or(0.5);
    
    let node_grouped = build_node_groups(&queries_to_iter, query_map);
    
    // One node group per task; each owns a dense result vector
    let thread_results: Vec<Vec<CompactEntry>> = node_grouped
        .par_iter().map(|(_node, chunk)| {
            let mut local_results: Vec<CompactEntry> = vec![CompactEntry::default(); num_ids];
            let mut node_occupancy = BitVector::new(num_ids);
            let mut edge_occupancy = BitVector::new(num_ids);
            let mut edge_exact = BitVector::new(num_ids);
            // Best substituted score per target on the current edge
            let mut edge_substituted: Vec<f32> = Vec::new();
            let mut substituted_hits: Vec<usize> = Vec::new();
            let mut prev_edge = None;

            for (_i, (e, query, weight, capped_idf)) in chunk.iter().enumerate() {
                let need_edge_update = prev_edge.map_or(true, |prev| prev != *e);
                
                if need_edge_update {
                    // Close the previous edge
                    if prev_edge.is_some() {
                        for nid in 0..num_ids {
                            if edge_occupancy.is_set(nid) {
                                local_results[nid].edge_count += 1;
                            }
                        }
                        close_substituted_edge(
                            &mut local_results, &edge_exact, &mut edge_substituted, &mut substituted_hits,
                        );
                    }
                    edge_occupancy.clear();
                    edge_exact.clear();
                    prev_edge = Some(*e);
                }
                
                let single_queried_values = index.get_entries(query.as_u32());
                let hash_count = single_queried_values.len();
                
                if let Some(freq_filter) = freq_filter {
                    if hash_count as f32 / lookup.len() as f32 > freq_filter {
                        continue;
                    }
                }

                let idf = if hash_count > 0 {
                    (lookup.len() as f32 / hash_count as f32).log2()
                } else {
                    continue;  // Hash absent from the index; nothing to count
                };
                let substituted = *weight < 1.0;
                if substituted && edge_substituted.is_empty() {
                    edge_substituted = vec![0.0; num_ids];
                }

                for &value in single_queried_values.iter() {
                    if value >= lookup.len() {
                        continue;
                    }
                
                    let nid = lookup.records()[value].id as usize;
                    let entry = &mut local_results[nid];

                    if !entry.initialized {
                        entry.initialized = true;
                        node_occupancy.set(nid);
                    }
                    entry.match_count += 1;
                    if substituted {
                        let score = weight * capped_idf;
                        if edge_substituted[nid] == 0.0 {
                            substituted_hits.push(nid);
                        }
                        edge_substituted[nid] = edge_substituted[nid].max(score);
                    } else {
                        entry.idf_sum += idf;
                        edge_exact.set(nid);
                    }
                    
                    edge_occupancy.set(nid);
                }
            }
            
            // Each node group contributes at most one node per target
            for nid in 0..num_ids {
                if node_occupancy.is_set(nid) {
                    local_results[nid].node_count = 1;
                }
            }
            
            // Close the last edge
            for nid in 0..num_ids {
                if edge_occupancy.is_set(nid) {
                    local_results[nid].edge_count += 1;
                }
            }
            close_substituted_edge(
                &mut local_results, &edge_exact, &mut edge_substituted, &mut substituted_hits,
            );
            local_results
        })
        .collect();

    // Merge node groups per target
    let results: Vec<(usize, StructureResult<'a>)> = (0..num_ids)
        .into_par_iter()
        .filter_map(|nid| {
            let mut merged_entry = CompactEntry::default();
            let mut found_data = false;
            
            for thread_array in &thread_results {
                let entry = &thread_array[nid];
                if entry.initialized {
                    if !found_data {
                        merged_entry = entry.clone();
                        found_data = true;
                    } else {
                        merged_entry.match_count += entry.match_count;
                        merged_entry.idf_sum += entry.idf_sum;
                        merged_entry.node_count += entry.node_count;
                        merged_entry.edge_count += entry.edge_count;
                    }
                }
            }
            
            if found_data && merged_entry.match_count > 0 {
                let lookup_entry = lookup.entry(nid);
                merged_entry.idf_sum *= (lookup_entry.nres as f32).powf(-lp);

                let sr = StructureResult::new(
                    lookup_entry.name,
                    nid,
                    merged_entry.match_count as usize,
                    merged_entry.node_count as usize,
                    merged_entry.edge_count as usize,
                    merged_entry.idf_sum,
                    lookup_entry.nres,
                    lookup_entry.plddt,
                    lookup_entry.db_key,
                );
                Some((nid, sr))
            } else {
                None
            }
        })
        .collect();

    results
}

/// Add each target's best substituted score on the edge just closed, unless the target
/// also hit the edge exactly, then reset the buffers.
fn close_substituted_edge(
    results: &mut [CompactEntry], edge_exact: &BitVector,
    edge_substituted: &mut [f32], substituted_hits: &mut Vec<usize>,
) {
    for &nid in substituted_hits.iter() {
        if !edge_exact.is_set(nid) {
            results[nid].idf_sum += edge_substituted[nid];
        }
        edge_substituted[nid] = 0.0;
    }
    substituted_hits.clear();
}

/// Keep the rarest hashes: a fraction (`sampling_ratio`) or a count (`sampling_count`).
/// With neither or both set, all hashes are kept.
fn sample_query(
    queries: &Vec<GeometricHash>, 
    index: &FolddiscoIndex,
    sampling_ratio: Option<f32>, 
    sampling_count: Option<usize>,
) -> Vec<GeometricHash> {
    match (sampling_ratio, sampling_count) {
        (None, None) => queries.clone(),
        (Some(sampling_ratio), None) => {
            let mut sampled_queries = queries.par_iter().map(|query| {
                let single_queried_values = index.get_entries(query.as_u32());
                let hash_count = single_queried_values.len();
                (query, hash_count)
            }).collect::<Vec<_>>();
            sampled_queries.sort_by(|a, b| a.1.cmp(&b.1));
            let sample_query_size: usize = sampled_queries.len();
            sampled_queries.truncate((sampling_ratio * sample_query_size as f32).ceil() as usize);
            sampled_queries.into_iter().map(|(query, _)| *query).collect()
        },
        (None, Some(sampling_count)) => {
            let mut sampled_queries = queries.par_iter().map(|query| {
                let single_queried_values = index.get_entries(query.as_u32());
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

/// Group hashes by the first residue of their edge, sorted by edge, with each hash's
/// weight and IDF from `query_map`.
fn build_node_groups(
    sampled_queries: &[GeometricHash],
    query_map: &QueryHashMap,
) -> HashMap<usize, Vec<((usize, usize), GeometricHash, f32, f32)>> {
    let mut node_groups: HashMap<usize, Vec<((usize, usize), GeometricHash, f32, f32)>> = HashMap::default();

    for &query in sampled_queries {
        if let Some(h) = query_map.get(&query) {
            node_groups.entry(h.pair.0).or_insert_with(Vec::new).push((h.pair, query, h.weight, h.idf));
        }
    }
    for (_, chunk) in node_groups.iter_mut() {
        chunk.sort_by_key(|(edge, _, _, _)| *edge);
    }
    node_groups
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::core::HashType;
    use crate::index::lookup::{load_lookup_from_file, save_lookup_to_file};
    use crate::controller::query::QueryHash;

    fn hash(value: u32) -> GeometricHash {
        GeometricHash::from_u32(value, HashType::PDBTrRosetta)
    }

    #[test]
    fn substituted_hits_count_once_per_edge_and_only_without_an_exact_hit() {
        let dir = std::env::temp_dir().join(format!("folddisco_count_query_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx").to_string_lossy().to_string();

        // Target 0: exact h0 and substituted h1, h2 on edge (0, 1)
        // Target 1: substituted h1, h2 on edge (0, 1); target 2: substituted h3 on edge (0, 2)
        let postings: Vec<(usize, Vec<u32>)> = vec![(0, vec![0, 1, 2]), (1, vec![1, 2]), (2, vec![3])];
        let mut index = FolddiscoIndex::new(8, prefix.clone(), false);
        for (id, hashes) in &postings {
            index.count_entries(hashes, *id);
        }
        index.allocate_entries();
        let mut buffer = Vec::new();
        for (id, hashes) in &postings {
            index.add_entries(hashes, *id, &mut buffer);
        }
        index.wrapup_offset_and_save_entries();
        index.prune_to_sparse();

        let lookup_path = format!("{}.lookup", prefix);
        let names: Vec<String> = (0..3).map(|i| format!("t{}", i)).collect();
        save_lookup_to_file(&lookup_path, &names, &vec![0, 1, 2], Some(&vec![1, 1, 1]), Some(&vec![0.0; 3]), None);
        let lookup = load_lookup_from_file(&lookup_path);

        let mut query_map = QueryHashMap::default();
        let entry = |pair, weight, idf| QueryHash { pair, weight, idf, observed: true };
        query_map.insert(hash(0), entry((0, 1), 1.0, 9.0));
        query_map.insert(hash(1), entry((0, 1), 0.5, 1.0));
        query_map.insert(hash(2), entry((0, 1), 0.5, 1.2));
        query_map.insert(hash(3), entry((0, 2), 0.25, 2.0));
        let queries = (0..4).map(hash).collect::<Vec<_>>();

        let mut results = count_query(&queries, &query_map, &index, &lookup, None, None, None, None);
        results.sort_by_key(|(id, _)| *id);
        let idf = |i: usize| results[i].1.idf;
        // Exact hit only; its IDF is recomputed from the index (3 structures, 1 hit)
        assert!((idf(0) - 3f32.log2()).abs() < 1e-5, "{}", idf(0));
        // Best substituted hash on the edge, once
        assert!((idf(1) - 0.6).abs() < 1e-5, "{}", idf(1));
        assert!((idf(2) - 0.5).abs() < 1e-5, "{}", idf(2));
        // Coverage counts every hit
        assert_eq!((results[0].1.edge_count, results[0].1.total_match_count), (1, 3));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
