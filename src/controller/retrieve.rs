///! Module for retrieval functions
///! Contains functions for retrieving target residues containing the motif

use std::collections::BTreeSet;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

use petgraph::Graph;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::structure::lms_qcp::LmsQcpSuperimposer;
use crate::structure::metrics::{PrecomputedDistances, StructureSimilarityMetrics};
use crate::utils::convert::{map_aa_to_u8, map_u8_to_aa}; 
use crate::prelude::*; 
use crate::structure::{coordinate::Coordinate, core::CompactStructure, kabsch::KabschSuperimposer}; 
use crate::utils::combination::{CombinationIterator, CombinationVecIterator};
use crate::controller::graph::{connected_components_with_given_node_count, create_index_graph};
use crate::controller::feature::get_single_feature;
use crate::controller::ResidueMatch;
use crate::controller::io::read_structure_from_path;

#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;

const PREFILTER_AA_SKIPPING_SIZE: usize = 200; // If query vector is larger than this, skip prefiltering amino acids
const RESIDUE_RESCUE_COUNT_CUTOFF: usize = 2; 

pub fn hash_vec_to_aa_pairs(hash_vec: &Vec<GeometricHash>) -> HashSet<(u32, u32)> {
    let mut output: HashSet<(u32, u32)> = HashSet::default();
    let mut feature = vec![0.0; 9];
    for hash in hash_vec {
        hash.reverse_hash_default(&mut feature);
        output.insert((feature[0] as u32, feature[1] as u32));
    }
    output
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


pub fn retrieve_with_prefilter(
    compact: &CompactStructure, hash_set: &HashSet<GeometricHash>, prefilter: CombinationVecIterator,
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
    dist_cutoff: f32, ca_distance_cutoff: f32,
    query_aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
) -> (Vec<(usize, usize, GeometricHash)>, Vec<(usize, (usize, usize))>) {
    let mut output: Vec<(usize, usize, GeometricHash)> = Vec::new();
    let mut candidate_pairs: Vec<(usize, (usize, usize))> = Vec::new();
    let hash = hash_set.iter().next().cloned().unwrap();
    let mut feature = vec![0.0; 9];
    if prefilter.is_empty() {
        let comb = CombinationIterator::new(compact.num_residues);
        comb.for_each(|(i, j)| {
            let is_feature = get_single_feature(i, j, compact, hash.hash_type(), dist_cutoff, &mut feature);
            if is_feature {
                // Check distance & if it is within the threshold, add to candidate_pairs
                let aa1 = map_aa_to_u8(&compact.residue_name[i]);
                let aa2 = map_aa_to_u8(&compact.residue_name[j]);
                if query_aa_dist_map.contains_key(&(aa1, aa2)) {
                    let dists = query_aa_dist_map.get(&(aa1, aa2)).unwrap();
                    let curr_dist = compact.get_ca_distance(i, j);
                    // If curr_dist is Some and within the threshold, add to candidate_pairs
                    if curr_dist.is_some() {
                        let curr_dist = curr_dist.unwrap();
                        for (dist, qi) in dists {
                            if (curr_dist - dist).abs() < ca_distance_cutoff {
                                candidate_pairs.push((*qi, (i, j)));
                            }
                        }
                    }
                }
                if let Some(multiple_bins) = multiple_bin {
                    for (nbin_dist, nbin_angle) in multiple_bins.iter() {
                        let curr_hash = GeometricHash::perfect_hash(&feature, hash.hash_type(), *nbin_dist, *nbin_angle);
                        if hash_set.contains(&curr_hash) {
                            output.push((i, j, curr_hash));
                        }
                    }
                } else {
                    let curr_hash = if nbin_dist == 0 || nbin_angle == 0 {
                        GeometricHash::perfect_hash_default(&feature, hash.hash_type())
                    } else {
                        GeometricHash::perfect_hash(&feature, hash.hash_type(), nbin_dist, nbin_angle)
                    };
                    if hash_set.contains(&curr_hash) {
                        output.push((i, j, curr_hash));
                    }
                }
            }
        });
    } else {
        for (i, j) in prefilter {
            let is_feature = get_single_feature(i, j, compact, hash.hash_type(), dist_cutoff, &mut feature);
            if is_feature {
                // Check distance & if it is within the threshold, add to candidate_pairs
                let aa1 = map_aa_to_u8(&compact.residue_name[i]);
                let aa2 = map_aa_to_u8(&compact.residue_name[j]);
                if query_aa_dist_map.contains_key(&(aa1, aa2)) {
                    let dists = query_aa_dist_map.get(&(aa1, aa2)).unwrap();
                    let curr_dist = compact.get_ca_distance(i, j);
                    // If curr_dist is Some and within the threshold, add to candidate_pairs
                    if curr_dist.is_some() {
                        let curr_dist = curr_dist.unwrap();
                        for (dist, qi) in dists {
                            if (curr_dist - dist).abs() < ca_distance_cutoff{
                                candidate_pairs.push((*qi, (i, j)));
                            }
                        }
                    }
                }
                if let Some(multiple_bins) = multiple_bin {
                    for (nbin_dist, nbin_angle) in multiple_bins.iter() {
                        let curr_hash = GeometricHash::perfect_hash(&feature, hash.hash_type(), *nbin_dist, *nbin_angle);
                        if hash_set.contains(&curr_hash) {
                            output.push((i, j, curr_hash));
                        }
                    }
                } else {
                    let curr_hash = if nbin_dist == 0 || nbin_angle == 0 {
                        GeometricHash::perfect_hash_default(&feature, hash.hash_type())
                    } else {
                        GeometricHash::perfect_hash(&feature, hash.hash_type(), nbin_dist, nbin_angle)
                    };
                    if hash_set.contains(&curr_hash) {
                        // output.push((*i, *j));
                        output.push((i, j, curr_hash));
                    }
                }
            }
        }
    }
    // Sort candidatpar_sort_byfirst element (query index)
    // candidate_pairs.par_sort_by_key(|(qi, _)| *qi);
    (output, candidate_pairs)
}


pub fn get_chain_and_res_ind(compact: &CompactStructure, i: usize) -> (u8, u64) {
    (compact.chain_per_residue[i], compact.residue_serial[i])
}
pub fn res_index_to_char(chain: u8, res_ind: u64) -> String {
    format!("{}{}", chain as char, res_ind)
}

#[cfg(feature = "foldcomp")]
pub fn retrieval_wrapper_for_foldcompdb(
    db_key: usize, node_count: usize, query_vector: &Vec<GeometricHash>,
    _hash_type: HashType, _nbin_dist: usize, _nbin_angle: usize,
    multiple_bin: &Option<Vec<(usize, usize)>>, dist_cutoff: f32,
    query_map: &HashMap<GeometricHash, ((usize, usize), bool, f32)>,
    query_structure: &CompactStructure, all_query_indices: &Vec<usize>,
    aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
    ca_distance_cutoff: f32, partial_fit: bool,
    foldcomp_db_reader: &FoldcompDbReader,
) -> (Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>, 
      Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>, usize, f32) {
    let compact = foldcomp_db_reader.read_single_structure_by_id(db_key).expect("Error reading structure from foldcomp db");
    let compact = compact.to_compact();

    // let mut indices_found: Vec<Vec<(usize, usize)>> = Vec::new();
    // Iterate over query vector and retrieve indices
    // Parallel
    let query_set: HashSet<GeometricHash> = HashSet::from_iter(query_vector.clone());
    let query_symmetry_map = get_hash_symmetry_map(&query_set);

    let aa_filter = if _hash_type.amino_acid_index().is_some() {
        let (index_set1, index_set2) = prefilter_amino_acid(&query_set, _hash_type, &compact);
        CombinationVecIterator::new_from_btreesets(&index_set1, &index_set2)
    } else {
        CombinationVecIterator::new(vec![], vec![])
    };
    
    // let (index_set1, index_set2) = prefilter_amino_acid(&query_set, _hash_type, &compact);
    // let aa_filter = CombinationVecIterator::new_from_btreesets(&index_set1, &index_set2);
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
    
    // Make a graph and find connected components with the same node count
    // NOTE: Naive implementation to find matching components. Need to be improved to handle partial matches
    let graph = create_index_graph(&indices_found);
    let connected = connected_components_with_given_node_count(&graph, node_count);
    
    // Parallel
    let output: Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32,
                     Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)> = connected.par_iter().map(|component| {
        // Filter graph to get subgraph with component
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
        
        // Calculate IDF for this subgraph
        let subgraph_idf = calculate_subgraph_idf(&subgraph, query_map);
        
        // Find mapping between query residues and retrieved residues
        let (query_indices, retrieved_indices) = map_query_and_retrieved_residues(
            &subgraph, query_map, node_count, &query_symmetry_map,
        );

        // Pre-build HashSets for O(1) lookups instead of O(n) vector operations
        let retrieved_indices_set: HashSet<usize> = retrieved_indices.iter().cloned().collect();
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
                    // Find and replace previous entry - more complex but still O(n) in worst case
                    // This case should be rare, so we keep it simple
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
                        // Use HashSet for O(1) lookup instead of O(n) vector contains
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
    // Split
    let (result_from_hash, result): (Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>,
        Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>) = output.into_iter().map(|(a, b, c, d, e, f, g, h, i, j, k, l, m, n)| {
        ((a, b, c, d, e, f, g), (h, i, j, k, l, m, n))
    }).unzip();
    // In result, find the maximum matching node count and minimum RMSD with max match
    let mut max_matching_node_count = 0;
    let mut min_rmsd_with_max_match = 0.0;
    result.iter().for_each(|(res_vec, rmsd, _, _, _, _, _)| {
        // Count number of Some in res_vec
        let count = res_vec.iter().filter(|&x| x.is_some()).count();
        if count > max_matching_node_count {
            max_matching_node_count = count;
            min_rmsd_with_max_match = *rmsd;
        } else if count == max_matching_node_count && *rmsd < min_rmsd_with_max_match {
            min_rmsd_with_max_match = *rmsd;
        }
    });
    (result_from_hash, result, max_matching_node_count, min_rmsd_with_max_match)
}




// Returns a vector of 
// 1) chain+residue index as String
// 2) RMSD value as f32
// TODO: Add U, T matrix and C-alpha coordinates
// 3) U, T matrix and C-alpha coordinates
pub fn retrieval_wrapper(
    path: &str, node_count: usize, query_vector: &Vec<GeometricHash>,
    _hash_type: HashType, _nbin_dist: usize, _nbin_angle: usize, 
    multiple_bin: &Option<Vec<(usize, usize)>>, dist_cutoff: f32,
    query_map: &HashMap<GeometricHash, ((usize, usize), bool, f32)>,
    query_structure: &CompactStructure, all_query_indices: &Vec<usize>,
    aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
    ca_distance_cutoff: f32, partial_fit: bool,
) -> (Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>, 
      Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>, usize, f32) {
    // Load structure to retrieve motif
    let compact = read_structure_from_path(&path).expect("Error reading structure from path");
    let compact = compact.to_compact();
    // let mut indices_found: Vec<Vec<(usize, usize)>> = Vec::new();
    // Iterate over query vector and retrieve indices
    // Parallel
    let query_set: HashSet<GeometricHash> = HashSet::from_iter(query_vector.clone());
    let query_symmetry_map = get_hash_symmetry_map(&query_set);
    
    // let (index_set1, index_set2) = prefilter_amino_acid(&query_set, _hash_type, &compact);

    let aa_filter = if _hash_type.amino_acid_index().is_some() {
        let (index_set1, index_set2) = prefilter_amino_acid(&query_set, _hash_type, &compact);
        CombinationVecIterator::new_from_btreesets(&index_set1, &index_set2)
    } else {
        CombinationVecIterator::new(vec![], vec![])
    };


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
    
    // Make a graph and find connected components with the same node count
    // NOTE: Naive implementation to find matching components. Need to be improved to handle partial matches
    let graph = create_index_graph(&indices_found);
    let connected = connected_components_with_given_node_count(&graph, node_count);

    // Parallel
// Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>)>
    let output: Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32,
                     Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)> = connected.par_iter().map(|component| {
        // Filter graph to get subgraph with component
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
        
        // Calculate IDF for this subgraph
        let subgraph_idf = calculate_subgraph_idf(&subgraph, query_map);
        
        // Find mapping between query residues and retrieved residues
        let (query_indices, retrieved_indices) = map_query_and_retrieved_residues(
            &subgraph, query_map, node_count, &query_symmetry_map,
        );

        // Pre-build HashSets for O(1) lookups instead of O(n) vector operations
        let retrieved_indices_set: HashSet<usize> = retrieved_indices.iter().cloned().collect();
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
                    // Find and replace previous entry - more complex but still O(n) in worst case
                    // This case should be rare, so we keep it simple
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
                        // Use HashSet for O(1) lookup instead of O(n) vector contains
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
    // Split 
    let (result_from_hash, result): (Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>, 
        Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>) = output.into_iter().map(|(a, b, c, d, e, f, g, h, i, j, k, l, m, n)| {
        ((a, b, c, d, e, f, g), (h, i, j, k, l, m, n))
    }).unzip();
    // In result, find the maximum matching node count and minimum RMSD with max match
    let mut max_matching_node_count = 0;
    let mut min_rmsd_with_max_match = 0.0;
    result.iter().for_each(|(res_vec, rmsd, _, _, _, _, _)| {
        // Count number of Some in res_vec
        let count = res_vec.iter().filter(|&x| x.is_some()).count();
        if count > max_matching_node_count {
            max_matching_node_count = count;
            min_rmsd_with_max_match = *rmsd;
        } else if count == max_matching_node_count && *rmsd < min_rmsd_with_max_match {
            min_rmsd_with_max_match = *rmsd;
        }
    });
    (result_from_hash, result, max_matching_node_count, min_rmsd_with_max_match)
}

fn get_hash_symmetry_map(query_set: &HashSet<GeometricHash>) -> HashMap<GeometricHash, bool> {
    let mut query_symmetry_map = HashMap::default();
    query_set.iter().for_each(|hash| {
        let symmetry = hash.is_symmetric();
        query_symmetry_map.insert(hash.clone(), symmetry);
    });
    query_symmetry_map
}

pub fn prefilter_amino_acid(query_set: &HashSet<GeometricHash>, _hash_type: HashType, compact: &CompactStructure) -> (BTreeSet<usize>, BTreeSet<usize>) {
    let mut observed_aa1: HashSet<u8> = HashSet::default();
    let mut observed_aa2: HashSet<u8> = HashSet::default();
    let mut index_vec1 = BTreeSet::new();
    let mut index_vec2 = BTreeSet::new();
    // If query_set is too large, just return empty sets
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

pub fn map_query_and_retrieved_residues(
    retrieved: &Graph<usize, GeometricHash>, 
    query_map: &HashMap<GeometricHash, ((usize, usize), bool, f32)>,
    node_count: usize,
    query_symmetry_map: &HashMap<GeometricHash, bool>,
) -> (Vec<usize>, Vec<usize>) {
    // Find max indices
    let max_query_idx = query_map.values().map(|((i, j), _, _)| (*i).max(*j)).max().unwrap_or(0);
    let max_retrieved_idx = retrieved.node_weights().max().copied().unwrap_or(0);
    
    let q_size = max_query_idx + 1;
    let r_size = max_retrieved_idx + 1;
    
    // Flat array for counts
    let mut counts = vec![0u8; q_size * r_size];
    
    // Track best match per query: (count, retrieved_idx)
    let mut best_match: Vec<(u8, usize)> = vec![(0, 0); q_size];
    
    // Helper macro for flat array access
    macro_rules! count_at {
        ($q:expr, $r:expr) => {
            counts[$q * r_size + $r]
        };
    }
    
    // Single pass: count votes and track best per query
    for edge in retrieved.edge_indices() {
        let (i, j) = retrieved.edge_endpoints(edge).unwrap();
        let hash = retrieved[edge];
        
        if let Some(&((query_i, query_j), _, _)) = query_map.get(&hash) {
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

/// Calculate total IDF score for a subgraph by summing IDFs of all edges
pub fn calculate_subgraph_idf(
    subgraph: &Graph<usize, GeometricHash>,
    query_map: &HashMap<GeometricHash, ((usize, usize), bool, f32)>,
) -> f32 {
    let mut total_idf = 0.0f32;
    
    for edge in subgraph.edge_indices() {
        let hash = subgraph[edge];
        if let Some(&(_, _, idf)) = query_map.get(&hash) {
            total_idf += idf;
        }
    }
    
    total_idf
}

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
    
    match lms {
        true => {
            if index1.len() <= 3 {
                let mut superposer = KabschSuperimposer::new();
                superposer.set_atoms(&coord_vec1, &coord_vec2);
                superposer.run();

                // Calculate target metrics
                let target_metrics = match (&superposer.reference_coords, &superposer.transformed_coords) {
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

                (superposer.get_rms(), superposer.rot.unwrap(), superposer.tran.unwrap(), target_calpha, target_metrics)
            } else {
                let mut superposer = LmsQcpSuperimposer::new();
                superposer.set_atoms(&coord_vec1, &coord_vec2);
                superposer.run();
                
                // Calculate target metrics
                let target_metrics = match (&superposer.reference_coords, &superposer.transformed_coords) {
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
                
                (superposer.get_rms_inliers(), superposer.rot.unwrap(), superposer.tran.unwrap(), target_calpha, target_metrics)
            }
        }
        false => {
            let mut superposer = KabschSuperimposer::new();
            superposer.set_atoms(&coord_vec1, &coord_vec2);
            superposer.run();
            // Calculate target metrics
            let target_metrics = match (&superposer.reference_coords, &superposer.transformed_coords) {
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
            (superposer.get_rms(), superposer.rot.unwrap(), superposer.tran.unwrap(), target_calpha, target_metrics)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::controller::query::make_query_map;

    use super::*;

    #[test]
    fn test_retrieval_wrapper() {
        let path = String::from("data/serine_peptidases/4cha.pdb");
        let query_string = "B57,B102,C195";
        let (query_residues, aa_substitutions) = parse_query_string(query_string, b'A');
        let hash_type = HashType::PDBTrRosetta;
        let nbin_dist = 16;
        let nbin_angle = 4;
        let dist_thresholds: Vec<f32> = vec![0.5,1.0];
        let angle_thresholds: Vec<f32> = vec![5.0,10.0];
        let dist_cutoff = 20.0;
        let (query_map, query_indices, aa_dist_map ) = make_query_map(
            &path, &query_residues, hash_type, nbin_dist, nbin_angle, &None,
            &dist_thresholds, &angle_thresholds, &aa_substitutions, dist_cutoff, false,
            &None, &None, 1000.0
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

}