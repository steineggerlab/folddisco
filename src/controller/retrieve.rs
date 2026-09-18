//! Residue matching: find target residues that realise the query motif, then superpose.

use std::collections::BTreeSet;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use petgraph::Graph;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::structure::chain_id::{format_chain_residue, ChainId};
use crate::structure::lms_qcp::LmsQcpSuperimposer;
use crate::structure::metrics::{deformation_stats_indexed, PrecomputedDistances, StructureSimilarityMetrics};
use crate::utils::convert::{map_aa_to_u8, map_u8_to_aa}; 
use crate::prelude::*; 
use crate::structure::{coordinate::Coordinate, core::CompactStructure, kabsch::KabschSuperimposer}; 
use crate::utils::combination::{CombinationIterator, CombinationVecIterator};
use crate::controller::graph::{connected_components_with_given_node_count, create_index_graph};
use crate::controller::feature::get_single_feature;
use crate::controller::ResidueMatch;
use crate::controller::query::{QueryHash, QueryHashMap};
use crate::controller::io::read_structure_from_path;

#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;

const PREFILTER_AA_SKIPPING_SIZE: usize = 200; // Above this, prefilter by amino acid code instead of name
const RESIDUE_RESCUE_COUNT_CUTOFF: usize = 2; 

/// Distinct (aa1, aa2) codes decoded from `hash_vec` with default binning.
pub fn hash_vec_to_aa_pairs(hash_vec: &Vec<GeometricHash>) -> HashSet<(u32, u32)> {
    let mut output: HashSet<(u32, u32)> = HashSet::default();
    let mut feature = vec![0.0; 9];
    for hash in hash_vec {
        hash.reverse_hash_default(&mut feature);
        output.insert((feature[0] as u32, feature[1] as u32));
    }
    output
}


/// Format residue pairs as `A12-B34`, comma separated.
pub fn res_vec_as_string(res_vec: &Vec<((u8, u8), (u64, u64))>) -> String {
    let mut output = String::new();
    for (k, (i, j)) in res_vec.iter().enumerate() {
        if k == res_vec.len() - 1 {
            output.push_str(&format!("{}{}-{}{}", i.0 as char, j.0, i.1 as char, j.1));
        } else {
            output.push_str(&format!("{}{}-{}{},", i.0 as char, j.0, i.1 as char, j.1));
        }
    }
    output
}

/// Target residue pairs whose hash is in `hash_set`, and (query residue, target pair)
/// candidates whose amino acids and Ca distance agree with `query_aa_dist_map`.
/// Only `prefilter` pairs are scanned when given; `None` scans every pair.
pub fn retrieve_with_prefilter(
    compact: &CompactStructure,
    hash_set: &HashSet<GeometricHash>,
    prefilter: Option<CombinationVecIterator>,
    nbin_dist: usize,
    nbin_angle: usize,
    multiple_bin: &Option<Vec<(usize, usize)>>,
    dist_cutoff: f32,
    ca_distance_cutoff: f32,
    query_aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
) -> (
    Vec<(usize, usize, GeometricHash)>,
    Vec<(usize, (usize, usize))>,
) {
    let mut output: Vec<(usize, usize, GeometricHash)> = Vec::new();
    let mut candidate_pairs: Vec<(usize, (usize, usize))> = Vec::new();

    let hash_type = match hash_set.iter().next() {
        Some(h) => h.hash_type(),
        None => return (output, candidate_pairs),
    };
    
    // Pre-convert all residues to u8
    let residue_types: Vec<u8> = compact.residue_name
        .iter()
        .map(|r| map_aa_to_u8(r))
        .collect();


    let mut temp_candidates: Vec<(usize, (usize, usize))> = Vec::with_capacity(16);
    let mut feature = vec![0.0; 9];

    let mut process_pair = |i: usize, j: usize| {
        // 1. Ca distance within cutoff
        let curr_dist = match compact.get_ca_distance(i, j) {
            Some(d) if d <= dist_cutoff => d,
            _ => return,
        };

        // 2. Amino acid pair present in the query, at a similar Ca distance
        let aa1 = residue_types[i];
        let aa2 = residue_types[j];
        
        let dists = match query_aa_dist_map.get(&(aa1, aa2)) {
            Some(d) => d,
            None => return, // AA pair not found in query
        };

        let mut is_valid_dist_for_query = false;
        temp_candidates.clear();
        
        for (dist, qi) in dists {
            if (curr_dist - dist).abs() < ca_distance_cutoff {
                temp_candidates.push((*qi, (i, j)));
                is_valid_dist_for_query = true;
            }
        }

        // Features are costlier than the distance check
        if !is_valid_dist_for_query {
            return;
        }

        // 3. Features and hashes
        let is_feature = get_single_feature(i, j, compact, hash_type, dist_cutoff, &mut feature);

        if is_feature {
            candidate_pairs.extend_from_slice(&temp_candidates);
            
            if let Some(multiple_bins) = multiple_bin {
                for (nb_d, nb_a) in multiple_bins.iter() {
                    let curr_hash =
                        GeometricHash::perfect_hash(&feature, hash_type, *nb_d, *nb_a);
                    if hash_set.contains(&curr_hash) {
                        output.push((i, j, curr_hash));
                    }
                }
            } else {
                let curr_hash = if nbin_dist == 0 || nbin_angle == 0 {
                    GeometricHash::perfect_hash_default(&feature, hash_type)
                } else {
                    GeometricHash::perfect_hash(&feature, hash_type, nbin_dist, nbin_angle)
                };
                if hash_set.contains(&curr_hash) {
                    output.push((i, j, curr_hash));
                }
            }
        }
    };

    match prefilter {
        None => CombinationIterator::new(compact.num_residues).for_each(|(i, j)| process_pair(i, j)),
        Some(pairs) => pairs.for_each(|(i, j)| process_pair(i, j)),
    }

    (output, candidate_pairs)
}


/// (chain, residue serial) of residue index `i`.
pub fn get_chain_and_res_ind(compact: &CompactStructure, i: usize) -> (ChainId, u64) {
    (compact.chain_per_residue[i], compact.residue_serial[i])
}
/// `A21`, or `AA_21` when the chain needs a separator.
pub fn res_index_to_char(chain: &ChainId, res_ind: u64) -> String {
    format_chain_residue(chain, res_ind, chain.needs_separator())
}

/// `retrieval_wrapper` for a target read from a Foldcomp database by `db_key`.
#[cfg(feature = "foldcomp")]
pub fn retrieval_wrapper_for_foldcompdb(
    db_key: usize, node_count: usize, query_vector: &Vec<GeometricHash>,
    _hash_type: HashType, _nbin_dist: usize, _nbin_angle: usize,
    multiple_bin: &Option<Vec<(usize, usize)>>, dist_cutoff: f32,
    query_map: &QueryHashMap,
    query_structure: &CompactStructure, all_query_indices: &Vec<usize>,
    aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
    ca_distance_cutoff: f32, partial_fit: bool,
    foldcomp_db_reader: &FoldcompDbReader,
) -> (Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>,
      Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>,
      usize, f32, f32) {
    let compact = foldcomp_db_reader.read_single_structure_by_id(db_key).expect("Error reading structure from foldcomp db");
    let compact = compact.to_compact();

    // let mut indices_found: Vec<Vec<(usize, usize)>> = Vec::new();
    let query_set: HashSet<GeometricHash> = HashSet::from_iter(query_vector.clone());
    let query_symmetry_map = get_hash_symmetry_map(&query_set);

    let aa_filter = amino_acid_prefilter(&query_set, _hash_type, &compact, aa_dist_map);
    let (indices_found , candidate_pairs) = retrieve_with_prefilter(
        &compact, &query_set, aa_filter, _nbin_dist, _nbin_angle, multiple_bin,
        dist_cutoff, ca_distance_cutoff, aa_dist_map
    );

    let candidate_pair_map: HashMap<usize, Vec<(usize, usize)>> = candidate_pairs.into_iter().fold(
        HashMap::default(), |mut map, (qi, pair)| {
            map.entry(qi).or_insert_with(Vec::new).push(pair);
            map
        }
    );
    
    // Connected components of the hit graph; naive, partial matches are handled poorly
    let graph = create_index_graph(&indices_found);
    let connected = connected_components_with_given_node_count(&graph, node_count);
    
    let output: Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32,
                     Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)> = connected.par_iter().map(|component| {
        let subgraph: Graph<usize, GeometricHash> = graph.filter_map(
            |node, _| {
                if component.contains(&graph[node]) {
                    Some(graph[node])
                } else {
                    None
                }
            },
            |_, edge| Some(*edge)
        );
        let node_count = subgraph.node_count();
        
        // Find mapping between query residues and retrieved residues
        let (query_indices, retrieved_indices) = map_query_and_retrieved_residues(
            &subgraph, query_map, node_count, &query_symmetry_map,
        );

        let retrieved_indices_set: HashSet<usize> = retrieved_indices.iter().cloned().collect();
        let subgraph_idf = calculate_subgraph_idf(&subgraph, query_map, &retrieved_indices_set);
        let query_to_retrieved: HashMap<usize, usize> = query_indices.iter()
            .zip(retrieved_indices.iter())
            .map(|(&q, &r)| (q, r))
            .collect();

        let mut query_indices_scanned: Vec<usize> = Vec::with_capacity(all_query_indices.len());
        let mut retrieved_indices_scanned: Vec<usize> = Vec::with_capacity(retrieved_indices.len());
        // let mut retrieved_indices_scanned_set: HashSet<usize> = HashSet::with_capacity(retrieved_indices.len());
        let mut retrieved_indices_scanned_set: HashSet<usize> = HashSet::default();
        
        // Sort component to match retrieved indices
        let mut res_vec: Vec<ResidueMatch> = Vec::with_capacity(all_query_indices.len());
        let mut res_vec_from_hash: Vec<ResidueMatch> = Vec::with_capacity(all_query_indices.len());
        let mut count_map: HashMap<usize, usize> = HashMap::default();
        let mut pairs_vec: Vec<(usize, usize)> = Vec::new();
        all_query_indices.iter().for_each(|&i| {
            // If i is in query_indices, get the corresponding retrieved index
            count_map.clear();
            if let Some(&retrieved_index) = query_to_retrieved.get(&i) {
                let (chain, res_ind) = get_chain_and_res_ind(&compact, retrieved_index);
                res_vec_from_hash.push(Some((chain, res_ind)));
                
                if !retrieved_indices_scanned_set.contains(&retrieved_index) {
                    res_vec.push(Some((chain, res_ind)));
                    query_indices_scanned.push(i);
                    retrieved_indices_scanned.push(retrieved_index);
                    retrieved_indices_scanned_set.insert(retrieved_index);
                } else {
                    // Rare: replace the previous entry
                    if let Some(prev_pos) = retrieved_indices_scanned.iter().position(|&x| x == retrieved_index) {
                        res_vec[prev_pos] = None;
                        res_vec.push(Some((chain, res_ind)));
                        query_indices_scanned.remove(prev_pos);
                        retrieved_indices_scanned.remove(prev_pos);
                        retrieved_indices_scanned_set.remove(&retrieved_index);
                        query_indices_scanned.push(i);
                        retrieved_indices_scanned.push(retrieved_index);
                        retrieved_indices_scanned_set.insert(retrieved_index);
                    }
                }
            } else {
                res_vec_from_hash.push(None);
                if candidate_pair_map.contains_key(&i) {
                    let pairs = candidate_pair_map.get(&i).unwrap().clone();
                    pairs_vec.clear();
                    pairs_vec.extend(pairs);
                    pairs_vec.sort_by_key(|&(j, _)| j);
                    let mut max_count = 0usize;
                    for (j, k) in &pairs_vec {
                        if retrieved_indices_set.contains(k) {
                            *count_map.entry(*j).or_insert(0) += 1;
                            // Track maximum count for this residue
                            if *count_map.get(j).unwrap() > max_count {
                                max_count = *count_map.get(j).unwrap();
                            }
                        }
                    }
                }
                if !count_map.is_empty() {
                    let max = count_map.iter().filter(|&(_, &v)| v == *count_map.values().max().unwrap())
                        .map(|(&k, &v)| (k, v))
                        .collect::<Vec<_>>();
                    
                    if max.len() == 1 && max[0].1 >= RESIDUE_RESCUE_COUNT_CUTOFF && !retrieved_indices_scanned_set.contains(&max[0].0) {
                        // If only one max entry and it has count > 1, add it
                        let (chain, res_ind) = get_chain_and_res_ind(&compact, max[0].0);
                        res_vec.push(Some((chain, res_ind)));
                        query_indices_scanned.push(i);
                        retrieved_indices_scanned.push(max[0].0);
                        retrieved_indices_scanned_set.insert(max[0].0);
                    } else {
                        res_vec.push(None);
                    }
                } else {
                    res_vec.push(None);
                }
            }
        });

        let (rmsd_from_hash, u_mat_from_hash, t_mat_from_hash, ca_coords_from_hash, metrics_from_hash) = rmsd_with_calpha_and_rottran(
            query_structure, &compact, &query_indices, &retrieved_indices, partial_fit
        );
        
        let (rmsd, u_mat, t_mat, ca_coords, metrics) = if res_vec == res_vec_from_hash {
            (rmsd_from_hash, u_mat_from_hash, t_mat_from_hash, ca_coords_from_hash.clone(), metrics_from_hash.clone())
        } else {
            rmsd_with_calpha_and_rottran(
                query_structure, &compact, &query_indices_scanned, &retrieved_indices_scanned, partial_fit
            )
        };
        
        (res_vec_from_hash, rmsd_from_hash, u_mat_from_hash, t_mat_from_hash, ca_coords_from_hash, metrics_from_hash, subgraph_idf,
         res_vec, rmsd, u_mat, t_mat, ca_coords, metrics, subgraph_idf)
    }).collect();
    let (result_from_hash, result): (Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>,
        Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>) = output.into_iter().map(|(a, b, c, d, e, f, g, h, i, j, k, l, m, n)| {
        ((a, b, c, d, e, f, g), (h, i, j, k, l, m, n))
    }).unzip();
    let (max_matching_node_count, min_rmsd_with_max_match, min_drmsd_with_max_match) =
        summarize_best_match(&result);
    (result_from_hash, result, max_matching_node_count, min_rmsd_with_max_match, min_drmsd_with_max_match)
}




/// Match the query motif in the structure at `path`.
///
/// Returns matches from hashes alone and after residue rescue, each as (residues, RMSD,
/// U, T, target Ca, metrics, IDF), plus max matched node count and its best RMSD/dRMSD.
pub fn retrieval_wrapper(
    path: &str, node_count: usize, query_vector: &Vec<GeometricHash>,
    _hash_type: HashType, _nbin_dist: usize, _nbin_angle: usize, 
    multiple_bin: &Option<Vec<(usize, usize)>>, dist_cutoff: f32,
    query_map: &QueryHashMap,
    query_structure: &CompactStructure, all_query_indices: &Vec<usize>,
    aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
    ca_distance_cutoff: f32, partial_fit: bool,
) -> (Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>,
      Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>,
      usize, f32, f32) {
    let compact = read_structure_from_path(&path).expect("Error reading structure from path");
    let compact = compact.to_compact();
    // let mut indices_found: Vec<Vec<(usize, usize)>> = Vec::new();
    let query_set: HashSet<GeometricHash> = HashSet::from_iter(query_vector.clone());
    let query_symmetry_map = get_hash_symmetry_map(&query_set);
    
    let aa_filter = amino_acid_prefilter(&query_set, _hash_type, &compact, aa_dist_map);

    let (indices_found , candidate_pairs) = retrieve_with_prefilter(
        &compact, &query_set, aa_filter, _nbin_dist, _nbin_angle,
        multiple_bin, dist_cutoff, ca_distance_cutoff, aa_dist_map
    );

    let candidate_pair_map: HashMap<usize, Vec<(usize, usize)>> = candidate_pairs.into_iter().fold(
        HashMap::default(), |mut map, (qi, pair)| {
            map.entry(qi).or_insert_with(Vec::new).push(pair);
            map
        }
    );
    
    // Connected components of the hit graph; naive, partial matches are handled poorly
    let graph = create_index_graph(&indices_found);
    let connected = connected_components_with_given_node_count(&graph, node_count);

// Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>)>
    let output: Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32,
                     Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)> = connected.par_iter().map(|component| {
        let subgraph: Graph<usize, GeometricHash> = graph.filter_map(
            |node, _| {
                if component.contains(&graph[node]) {
                    Some(graph[node])
                } else {
                    None
                }
            },
            |_, edge| Some(*edge)
        );
        
        let node_count = subgraph.node_count();
        
        // Find mapping between query residues and retrieved residues
        let (query_indices, retrieved_indices) = map_query_and_retrieved_residues(
            &subgraph, query_map, node_count, &query_symmetry_map,
        );

        let retrieved_indices_set: HashSet<usize> = retrieved_indices.iter().cloned().collect();
        let subgraph_idf = calculate_subgraph_idf(&subgraph, query_map, &retrieved_indices_set);
        let query_to_retrieved: HashMap<usize, usize> = query_indices.iter()
            .zip(retrieved_indices.iter())
            .map(|(&q, &r)| (q, r))
            .collect();

        let mut query_indices_scanned: Vec<usize> = Vec::with_capacity(all_query_indices.len());
        let mut retrieved_indices_scanned: Vec<usize> = Vec::with_capacity(all_query_indices.len());
        let mut retrieved_indices_scanned_set: HashSet<usize> = HashSet::default();
        retrieved_indices_scanned_set.reserve(all_query_indices.len());

        // Sort component to match retrieved indices
        let mut res_vec: Vec<ResidueMatch> = Vec::with_capacity(all_query_indices.len());
        let mut res_vec_from_hash: Vec<ResidueMatch> = Vec::with_capacity(all_query_indices.len());
        let mut count_map: HashMap<usize, usize> = HashMap::default();

        all_query_indices.iter().for_each(|&i| {
            // If i is in query_indices, get the corresponding retrieved index
            count_map.clear();
            if let Some(&retrieved_index) = query_to_retrieved.get(&i) {
                let (chain, res_ind) = get_chain_and_res_ind(&compact, retrieved_index);
                res_vec_from_hash.push(Some((chain, res_ind)));
                
                if !retrieved_indices_scanned_set.contains(&retrieved_index) {
                    res_vec.push(Some((chain, res_ind)));
                    query_indices_scanned.push(i);
                    retrieved_indices_scanned.push(retrieved_index);
                    retrieved_indices_scanned_set.insert(retrieved_index);
                } else {
                    // Rare: replace the previous entry
                    if let Some(prev_pos) = retrieved_indices_scanned.iter().position(|&x| x == retrieved_index) {
                        res_vec[prev_pos] = None;
                        res_vec.push(Some((chain, res_ind)));
                        query_indices_scanned.remove(prev_pos);
                        retrieved_indices_scanned.remove(prev_pos);
                        retrieved_indices_scanned_set.remove(&retrieved_index);
                        query_indices_scanned.push(i);
                        retrieved_indices_scanned.push(retrieved_index);
                        retrieved_indices_scanned_set.insert(retrieved_index);
                    }
                }
            } else {
                res_vec_from_hash.push(None);
                if candidate_pair_map.contains_key(&i) {
                    let pairs = candidate_pair_map.get(&i).unwrap();
                    let mut max_count = 0usize;
                    for (j, k) in pairs {
                        if retrieved_indices_set.contains(k) {
                            *count_map.entry(*j).or_insert(0) += 1;
                            // Track maximum count for this residue
                            if *count_map.get(j).unwrap() > max_count {
                                max_count = *count_map.get(j).unwrap();
                            }
                        }
                    }
                }
                if !count_map.is_empty() {
                    // let max = count_map.iter().max_by(|a, b| a.1.cmp(b.1)).unwrap();
                    // Get all max entries as a vector
                    let max = count_map.iter().filter(|&(_, &v)| v == *count_map.values().max().unwrap())
                        .map(|(&k, &v)| (k, v))
                        .collect::<Vec<_>>();

                    if max.len() == 1 && max[0].1 >= RESIDUE_RESCUE_COUNT_CUTOFF && !retrieved_indices_scanned_set.contains(&max[0].0) {
                        // If only one max entry and it has count > 1, add it
                        let (chain, res_ind) = get_chain_and_res_ind(&compact, max[0].0);
                        res_vec.push(Some((chain, res_ind)));
                        query_indices_scanned.push(i);
                        retrieved_indices_scanned.push(max[0].0);
                        retrieved_indices_scanned_set.insert(max[0].0);
                    } else {
                        res_vec.push(None);
                    }
                } else {
                    res_vec.push(None);
                }
            }
        });

        let (rmsd_from_hash, u_mat_from_hash, t_mat_from_hash, ca_coords_from_hash, metrics_from_hash) = rmsd_with_calpha_and_rottran(
            query_structure, &compact, &query_indices, &retrieved_indices, partial_fit
        );
        
        let (rmsd, u_mat, t_mat, ca_coords, metrics) = if res_vec == res_vec_from_hash {
            (rmsd_from_hash, u_mat_from_hash, t_mat_from_hash, ca_coords_from_hash.clone(), metrics_from_hash.clone())
        } else {
            rmsd_with_calpha_and_rottran(
                query_structure, &compact, &query_indices_scanned, &retrieved_indices_scanned, partial_fit
            )
        };

        (res_vec_from_hash, rmsd_from_hash, u_mat_from_hash, t_mat_from_hash, ca_coords_from_hash, metrics_from_hash, subgraph_idf,
         res_vec, rmsd, u_mat, t_mat, ca_coords, metrics, subgraph_idf)
        }).collect();
    let (result_from_hash, result): (Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>, 
        Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>) = output.into_iter().map(|(a, b, c, d, e, f, g, h, i, j, k, l, m, n)| {
        ((a, b, c, d, e, f, g), (h, i, j, k, l, m, n))
    }).unzip();
    let (max_matching_node_count, min_rmsd_with_max_match, min_drmsd_with_max_match) =
        summarize_best_match(&result);
    (result_from_hash, result, max_matching_node_count, min_rmsd_with_max_match, min_drmsd_with_max_match)
}

/// Largest matched residue count, with the lowest RMSD and dRMSD at that count.
fn summarize_best_match(
    result: &[(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)]
) -> (usize, f32, f32) {
    let mut max_matching_node_count = 0;
    let mut min_rmsd = 0.0f32;
    let mut min_drmsd = 0.0f32;
    for (res_vec, rmsd, _, _, _, metrics, _) in result.iter() {
        let count = res_vec.iter().filter(|&x| x.is_some()).count();
        if count > max_matching_node_count {
            max_matching_node_count = count;
            min_rmsd = *rmsd;
            min_drmsd = metrics.drmsd;
        } else if count == max_matching_node_count {
            if *rmsd < min_rmsd {
                min_rmsd = *rmsd;
            }
            if metrics.drmsd < min_drmsd {
                min_drmsd = metrics.drmsd;
            }
        }
    }
    (max_matching_node_count, min_rmsd, min_drmsd)
}

/// Whether each query hash is symmetric under swapping its residue pair.
fn get_hash_symmetry_map(query_set: &HashSet<GeometricHash>) -> HashMap<GeometricHash, bool> {
    let mut query_symmetry_map = HashMap::default();
    query_set.iter().for_each(|hash| {
        let symmetry = hash.is_symmetric();
        query_symmetry_map.insert(hash.clone(), symmetry);
    });
    query_symmetry_map
}

/// Target residue pairs to scan, or `None` for all of them.
///
/// Up to `PREFILTER_AA_SKIPPING_SIZE` hashes, residues are picked by the names in the query
/// hashes (an empty pick scans everything, as before). Above it, by the amino acid codes
/// in `query_aa_dist_map`: the same pairs, in the same order, that a full scan keeps.
fn amino_acid_prefilter(
    query_set: &HashSet<GeometricHash>, hash_type: HashType, compact: &CompactStructure,
    query_aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
) -> Option<CombinationVecIterator> {
    hash_type.amino_acid_index()?;
    if query_set.len() <= PREFILTER_AA_SKIPPING_SIZE {
        let (index_set1, index_set2) = prefilter_amino_acid(query_set, hash_type, compact);
        let pairs = CombinationVecIterator::new_from_btreesets(&index_set1, &index_set2);
        return (!pairs.is_empty()).then_some(pairs);
    }
    let mut first = [false; 256];
    let mut second = [false; 256];
    for &(aa1, aa2) in query_aa_dist_map.keys() {
        first[aa1 as usize] = true;
        second[aa2 as usize] = true;
    }
    let (mut index_vec1, mut index_vec2) = (Vec::new(), Vec::new());
    for (i, name) in compact.residue_name.iter().enumerate() {
        let code = map_aa_to_u8(name) as usize;
        if first[code] {
            index_vec1.push(i);
        }
        if second[code] {
            index_vec2.push(i);
        }
    }
    Some(CombinationVecIterator::new(index_vec1, index_vec2))
}

/// Target residue indices whose amino acid appears first / second in any query hash.
/// Both empty (no prefilter) above `PREFILTER_AA_SKIPPING_SIZE` hashes.
pub fn prefilter_amino_acid(query_set: &HashSet<GeometricHash>, _hash_type: HashType, compact: &CompactStructure) -> (BTreeSet<usize>, BTreeSet<usize>) {
    let mut observed_aa1: HashSet<u8> = HashSet::default();
    let mut observed_aa2: HashSet<u8> = HashSet::default();
    let mut index_vec1 = BTreeSet::new();
    let mut index_vec2 = BTreeSet::new();
    if query_set.len() > PREFILTER_AA_SKIPPING_SIZE {
        return (index_vec1, index_vec2);
    }

    let mut feature_holder = vec![0.0; 9];
    query_set.iter().for_each(|hash| {
        hash.reverse_hash_default(&mut feature_holder);
        let aa1 = feature_holder[_hash_type.amino_acid_index().unwrap()[0]] as u8;
        let aa2 = feature_holder[_hash_type.amino_acid_index().unwrap()[1]] as u8;
        if !observed_aa1.contains(&aa1) {
            observed_aa1.insert(aa1);
            let indices: Vec<usize> = compact.residue_name.iter().enumerate().filter_map(|(i, &res)| {
                if res == map_u8_to_aa(aa1).as_bytes() {
                    Some(i)
                } else {
                    None
                }
            }).collect();
            index_vec1.extend(indices);
        }
        if !observed_aa2.contains(&aa2) {
            observed_aa2.insert(aa2);
            let indices: Vec<usize> = compact.residue_name.iter().enumerate().filter_map(|(i, &res)| {
                if res == map_u8_to_aa(aa2).as_bytes() {
                    Some(i)
                } else {
                    None
                }
            }).collect();
            index_vec2.extend(indices);
        }
    });
    (index_vec1, index_vec2)
}

/// Greedily pair query and target residues by how many hit edges support each pairing,
/// stopping at `node_count` pairs.
pub fn map_query_and_retrieved_residues(
    retrieved: &Graph<usize, GeometricHash>, 
    query_map: &QueryHashMap,
    node_count: usize,
    query_symmetry_map: &HashMap<GeometricHash, bool>,
) -> (Vec<usize>, Vec<usize>) {
    let max_query_idx = query_map.values().map(|h| h.pair.0.max(h.pair.1)).max().unwrap_or(0);
    let max_retrieved_idx = retrieved.node_weights().max().copied().unwrap_or(0);
    
    let q_size = max_query_idx + 1;
    let r_size = max_retrieved_idx + 1;
    
    // Flat array for counts
    let mut counts = vec![0u8; q_size * r_size];
    
    // Track best match per query: (count, retrieved_idx)
    let mut best_match: Vec<(u8, usize)> = vec![(0, 0); q_size];
    
    macro_rules! count_at {
        ($q:expr, $r:expr) => {
            counts[$q * r_size + $r]
        };
    }
    
    // Single pass: count votes and track best per query
    for edge in retrieved.edge_indices() {
        let (i, j) = retrieved.edge_endpoints(edge).unwrap();
        let hash = retrieved[edge];
        
        if let Some(&QueryHash { pair: (query_i, query_j), .. }) = query_map.get(&hash) {
            let is_symmetric = *query_symmetry_map.get(&hash).unwrap();
            
            let pairs = if is_symmetric {
                let (q1, q2, r1, r2) = if query_i < query_j {
                    if retrieved[i] < retrieved[j] {
                        (query_i, query_j, retrieved[i], retrieved[j])
                    } else {
                        (query_i, query_j, retrieved[j], retrieved[i])
                    }
                } else {
                    if retrieved[i] < retrieved[j] {
                        (query_j, query_i, retrieved[i], retrieved[j])
                    } else {
                        (query_j, query_i, retrieved[j], retrieved[i])
                    }
                };
                [(q1, r1), (q2, r2)]
            } else {
                [(query_i, retrieved[i]), (query_j, retrieved[j])]
            };
            
            for (q, r) in pairs {
                count_at!(q, r) = count_at!(q, r).saturating_add(1);
                let new_count = count_at!(q, r);
                
                // Update best match for this query if this is better
                if new_count > best_match[q].0  || (new_count == best_match[q].0 && r < best_match[q].1) {
                    best_match[q] = (new_count, r);
                }
            }
        }
    }
    
    // Bucket sort by count
    const MAX_BUCKETS: usize = 256;
    let mut buckets: Vec<Vec<(usize, usize)>> = vec![Vec::new(); MAX_BUCKETS];
    
    for (q, &(count, r)) in best_match.iter().enumerate() {
        if count > 0 {
            let bucket_idx = (count as usize).min(MAX_BUCKETS - 1);
            buckets[bucket_idx].push((q, r));
        }
    }
    
    // Greedy assignment from highest to lowest count
    let mut query_indices: Vec<usize> = Vec::with_capacity(node_count);
    let mut retrieved_indices: Vec<usize> = Vec::with_capacity(node_count);
    let mut query_used = vec![false; q_size];
    let mut retrieved_used = vec![false; r_size];
    
    for bucket in buckets.iter().rev() {
        for &(q, r) in bucket {
            if !query_used[q] && !retrieved_used[r] {
                query_indices.push(q);
                retrieved_indices.push(r);
                query_used[q] = true;
                retrieved_used[r] = true;
                
                if query_indices.len() == node_count {
                    return (query_indices, retrieved_indices);
                }
            }
        }
    }
    
    (query_indices, retrieved_indices)
}

/// Sum of weighted query IDFs over the subgraph's edges. Edges from substituted or
/// neighbouring-bin hashes count only between `matched` residues, so a large component of
/// similar residues adds nothing.
pub fn calculate_subgraph_idf(
    subgraph: &Graph<usize, GeometricHash>,
    query_map: &QueryHashMap,
    matched: &HashSet<usize>,
) -> f32 {
    let mut total_idf = 0.0f32;
    
    for edge in subgraph.edge_indices() {
        let hash = subgraph[edge];
        if let Some(h) = query_map.get(&hash) {
            if !h.is_exact() {
                let (i, j) = subgraph.edge_endpoints(edge).unwrap();
                if !matched.contains(&subgraph[i]) || !matched.contains(&subgraph[j]) {
                    continue;
                }
            }
            total_idf += h.weight * h.idf;
        }
    }
    
    total_idf
}

/// RMSD of matched Ca+Cb coordinates; `lms` uses LMS for more than 3 residues.
pub fn rmsd_for_matched(
    compact1: &CompactStructure, compact2: &CompactStructure, 
    index1: &Vec<usize>, index2: &Vec<usize>, lms: bool
) -> f32 {
    let coord_vec1: Vec<Coordinate> = index1.iter().map(
        |&i| (compact1.ca_vector.get_coord(i).unwrap(), compact1.cb_vector.get_coord(i).unwrap())
    ).flat_map(|(a, b)| vec![a, b]).collect();
    
    let coord_vec2: Vec<Coordinate> = index2.iter().map(
        |&i| (compact2.ca_vector.get_coord(i).unwrap(), compact2.cb_vector.get_coord(i).unwrap())
    ).flat_map(|(a, b)| vec![a, b]).collect();

    match lms {
        true => {
            if index1.len() <= 3 {
                let mut superposer = KabschSuperimposer::new();
                superposer.set_atoms(&coord_vec1, &coord_vec2);
                superposer.run();
                superposer.get_rms()
            } else {
                let mut superposer = LmsQcpSuperimposer::new();
                superposer.set_atoms(&coord_vec1, &coord_vec2);
                superposer.run();
                superposer.get_rms_inliers()
            }
        }
        false => {
            let mut superposer = KabschSuperimposer::new();
            superposer.set_atoms(&coord_vec1, &coord_vec2);
            superposer.run();
            superposer.get_rms()
        }
    }
}

/// Like `rmsd_for_matched`, also returning rotation, translation, target Ca and metrics.
pub fn rmsd_with_calpha_and_rottran(
    compact1: &CompactStructure, compact2: &CompactStructure, 
    index1: &Vec<usize>, index2: &Vec<usize>, lms: bool
) -> (f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics) {

    let coord_vec1: Vec<Coordinate> = index1.iter().map(
        |&i| (compact1.ca_vector.get_coord(i).unwrap(), compact1.cb_vector.get_coord(i).unwrap())
    ).flat_map(|(a, b)| vec![a, b]).collect();
    
    let coord_vec2: Vec<Coordinate> = index2.iter().map(
        |&i| (compact2.ca_vector.get_coord(i).unwrap(), compact2.cb_vector.get_coord(i).unwrap())
    ).flat_map(|(a, b)| vec![a, b]).collect();

    let target_calpha: Vec<Coordinate> = index2.iter().map(
        |&i| compact2.ca_vector.get_coord(i).unwrap()
    ).collect();

    // Superposition-free deformation from the original coordinates (hinges stay low)
    let deformation = || deformation_stats_indexed(
        coord_vec1.len().min(coord_vec2.len()),
        |i, j| coord_vec1[i].calc_distance(&coord_vec1[j]),
        |i, j| coord_vec2[i].calc_distance(&coord_vec2[j]),
    );

    let metrics_from = |reference: &Option<Vec<[f32; 3]>>, transformed: &Option<Vec<[f32; 3]>>| {
        let mut metrics = match (reference, transformed) {
            (Some(ref_coords), Some(trans_coords)) => {
                let precomputed_distances = PrecomputedDistances::new(
                    &ref_coords, &trans_coords
                );
                let mut metrics = StructureSimilarityMetrics::new();
                metrics.calculate_all(&precomputed_distances);
                metrics
            },
            _ => StructureSimilarityMetrics::new(),
        };
        let (drmsd, max_dist_deviation) = deformation();
        metrics.drmsd = drmsd;
        metrics.max_dist_deviation = max_dist_deviation;
        metrics
    };

    match lms {
        true => {
            if index1.len() <= 3 {
                let mut superposer = KabschSuperimposer::new();
                superposer.set_atoms(&coord_vec1, &coord_vec2);
                superposer.run();
                let target_metrics = metrics_from(
                    &superposer.reference_coords, &superposer.transformed_coords
                );
                (superposer.get_rms(), superposer.rot.unwrap(), superposer.tran.unwrap(), target_calpha, target_metrics)
            } else {
                let mut superposer = LmsQcpSuperimposer::new();
                superposer.set_atoms(&coord_vec1, &coord_vec2);
                superposer.run();
                let target_metrics = metrics_from(
                    &superposer.reference_coords, &superposer.transformed_coords
                );
                (superposer.get_rms_inliers(), superposer.rot.unwrap(), superposer.tran.unwrap(), target_calpha, target_metrics)
            }
        }
        false => {
            let mut superposer = KabschSuperimposer::new();
            superposer.set_atoms(&coord_vec1, &coord_vec2);
            superposer.run();
            let target_metrics = metrics_from(
                &superposer.reference_coords, &superposer.transformed_coords
            );
            (superposer.get_rms(), superposer.rot.unwrap(), superposer.tran.unwrap(), target_calpha, target_metrics)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::controller::query::make_query_map;

    use super::*;

    use crate::controller::expand::ToleranceConfig;

    #[test]
    fn non_exact_edges_score_only_between_matched_residues() {
        let hash = |v: u32| GeometricHash::from_u32(v, HashType::PDBTrRosetta);
        let mut graph: Graph<usize, GeometricHash> = Graph::new();
        let nodes: Vec<_> = [1usize, 2, 3, 9].iter().map(|&r| graph.add_node(r)).collect();
        graph.add_edge(nodes[0], nodes[1], hash(0)); // exact
        graph.add_edge(nodes[1], nodes[2], hash(1)); // substituted, matched
        graph.add_edge(nodes[2], nodes[3], hash(2)); // substituted, residue 9 unmatched
        graph.add_edge(nodes[0], nodes[3], hash(3)); // exact, counted as before
        let entry = |pair, weight, idf, observed| QueryHash { pair, weight, idf, observed };
        let mut query_map = QueryHashMap::default();
        query_map.insert(hash(0), entry((0, 1), 1.0, 2.0, true));
        query_map.insert(hash(1), entry((1, 2), 0.5, 4.0, true));
        query_map.insert(hash(2), entry((1, 2), 0.5, 8.0, true));
        query_map.insert(hash(3), entry((0, 2), 1.0, 1.0, true));
        let matched: HashSet<usize> = [1, 2, 3].into_iter().collect();
        assert_eq!(calculate_subgraph_idf(&graph, &query_map, &matched), 2.0 + 0.5 * 4.0 + 1.0);
        // A neighbouring bin is treated like a substitution
        query_map.insert(hash(3), entry((0, 2), 1.0, 1.0, false));
        assert_eq!(calculate_subgraph_idf(&graph, &query_map, &matched), 2.0 + 0.5 * 4.0);
    }

    #[test]
    fn test_retrieval_wrapper() {
        let path = String::from("data/serine_peptidases/4cha.pdb");
        let query_string = "B57,B102,C195";
        let (query_residues, aa_substitutions) = parse_query_string(query_string, ChainId::from_byte(b'A'));
        let hash_type = HashType::PDBTrRosetta;
        let nbin_dist = 16;
        let nbin_angle = 4;
        let tolerance = ToleranceConfig::new(vec![0.5, 1.0], vec![5.0, 10.0], 1);
        let dist_cutoff = 20.0;
        let (query_map, query_indices, aa_dist_map ) = make_query_map(
            &path, &query_residues, hash_type, nbin_dist, nbin_angle, &None,
            &tolerance, &aa_substitutions, dist_cutoff, false,
            &None, 1000.0
        );
        let queries: Vec<GeometricHash> = query_map.keys().cloned().collect();
        let compact = read_structure_from_path(&path).expect("Error reading structure from path");
        let compact = compact.to_compact();
        let new_path = String::from("data/serine_peptidases/4cha.pdb");
        let output = measure_time!(retrieval_wrapper(
            &new_path, query_residues.len(), &queries, hash_type, nbin_dist, nbin_angle, &None,
            dist_cutoff, &query_map, &compact, &query_indices, &aa_dist_map, 1.5, false,
        ));
        println!("{:?}", output);
    }

    /// Max matched node count of the 4CHA catalytic triad query against `target`.
    fn triad_match_nodes(query_string: &str, target: &str) -> usize {
        use crate::controller::query::resolve_query_substitutions;
        use crate::controller::substitution::SubstitutionScheme;
        let path = String::from("data/serine_peptidases/4cha.pdb");
        let query = read_structure_from_path(&path).unwrap().to_compact();
        let (residues, subs) = parse_query_string(query_string, ChainId::from_byte(b'A'));
        let subs = resolve_query_substitutions(&query, &residues, &subs, SubstitutionScheme::Blosum62, false, false);
        let (query_map, query_indices, aa_dist_map) = make_query_map(
            &path, &residues, HashType::PDBTrRosetta, 16, 4, &None,
            &ToleranceConfig::default_query(), &subs, 20.0, false, &None, 1000.0,
        );
        let hashes: Vec<GeometricHash> = query_map.keys().cloned().collect();
        retrieval_wrapper(
            target, residues.len(), &hashes, HashType::PDBTrRosetta, 16, 4, &None, 20.0,
            &query_map, &query, &query_indices, &aa_dist_map, 1.0, false,
        ).2
    }

    #[test]
    fn substituted_residue_is_matched_not_only_looked_up() {
        // 4CHA (two copies) with Asp102 renamed to Asn: same geometry, different residue
        let pdb = std::fs::read_to_string("data/serine_peptidases/4cha.pdb").unwrap();
        let mutated: String = pdb.lines().map(|line| {
            if line.starts_with("ATOM") && line.len() > 26 && &line[17..20] == "ASP" && &line[22..26] == " 102" {
                format!("{}ASN{}\n", &line[..17], &line[20..])
            } else {
                format!("{}\n", line)
            }
        }).collect();
        let target = std::env::temp_dir().join(format!("folddisco_4cha_d102n_{}.pdb", std::process::id()));
        std::fs::write(&target, mutated).unwrap();
        let target = target.to_str().unwrap().to_string();

        assert_eq!(triad_match_nodes("B57,B102,C195", "data/serine_peptidases/4cha.pdb"), 3);
        assert!(triad_match_nodes("B57,B102,C195", &target) < 3);
        assert_eq!(triad_match_nodes("B57,B102:N,C195", &target), 3, "explicit substitution");
        assert_eq!(triad_match_nodes("B57,B102:*,C195", &target), 3, "BLOSUM62 scheme (D -> N, E)");
        std::fs::remove_file(&target).ok();
    }

    #[test]
    fn self_match_has_zero_deformation() {
        // Self match: RMSD and dRMSD vanish
        let path = String::from("data/serine_peptidases/4cha.pdb");
        let compact = read_structure_from_path(&path)
            .expect("Error reading structure from path").to_compact();
        let indices: Vec<usize> = vec![(ChainId::from(b'B'), 57u64), (b'B'.into(), 102), (b'C'.into(), 195)].iter()
            .map(|(chain, res)| compact.get_index(chain, res).expect("residue not found"))
            .collect();
        let (rmsd, _u, _t, _ca, metrics) = rmsd_with_calpha_and_rottran(
            &compact, &compact, &indices, &indices, false
        );
        assert!(rmsd < 1e-3, "rmsd {}", rmsd);
        assert!(metrics.drmsd < 1e-3, "drmsd {}", metrics.drmsd);
        assert!(metrics.max_dist_deviation < 1e-3, "max dev {}", metrics.max_dist_deviation);
    }

    #[test]
    fn deformation_metrics_ignore_rigid_motion_but_not_bending() {
        let stats = |a: &[[f32; 3]], b: &[[f32; 3]]| {
            let dist = |p: [f32; 3], q: [f32; 3]| {
                ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
            };
            deformation_stats_indexed(a.len(), |i, j| dist(a[i], a[j]), |i, j| dist(b[i], b[j]))
        };
        let reference = [
            [0.0f32, 0.0, 0.0], [5.0, 0.0, 0.0], [5.0, 5.0, 0.0], [0.0, 5.0, 0.0],
        ];
        // Rotate 90 degrees about z and translate: internal distances unchanged
        let moved = [
            [10.0f32, 0.0, 3.0], [10.0, 5.0, 3.0], [5.0, 5.0, 3.0], [5.0, 0.0, 3.0],
        ];
        let (drmsd, worst) = stats(&reference, &moved);
        assert!(drmsd < 1e-3 && worst < 1e-3, "rigid motion showed up: {} {}", drmsd, worst);
        // Pull one point away: the deformation now shows up
        let bent = [
            [0.0f32, 0.0, 0.0], [5.0, 0.0, 0.0], [5.0, 5.0, 0.0], [0.0, 8.0, 0.0],
        ];
        let (drmsd, worst) = stats(&reference, &bent);
        assert!(drmsd > 1.0, "drmsd {}", drmsd);
        assert!(worst >= 3.0 - 1e-3, "worst {}", worst);
    }
}