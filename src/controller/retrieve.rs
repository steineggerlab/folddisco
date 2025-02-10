// 
use std::collections::{BTreeSet, HashMap, HashSet};
use petgraph::Graph;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::utils::convert::{map_aa_to_u8, map_u8_to_aa}; 
use crate::prelude::*; 
use crate::structure::{coordinate::Coordinate, core::CompactStructure, qcp::QCPSuperimposer}; 
use crate::utils::combination::{CombinationIterator, CombinationVecIterator};
use crate::controller::graph::{connected_components_with_given_node_count, create_index_graph};
use crate::controller::feature::get_single_feature;
use crate::controller::ResidueMatch;
use crate::controller::io::read_structure_from_path;

#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;



pub fn hash_vec_to_aa_pairs(hash_vec: &Vec<GeometricHash>) -> HashSet<(u32, u32)> {
    let mut output: HashSet<(u32, u32)> = HashSet::new();
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
    path: &str, node_count: usize, query_vector: &Vec<GeometricHash>,
    _hash_type: HashType, _nbin_dist: usize, _nbin_angle: usize,
    multiple_bin: &Option<Vec<(usize, usize)>>, dist_cutoff: f32,
    query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    query_structure: &CompactStructure, all_query_indices: &Vec<usize>,
    aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
    ca_distance_cutoff: f32, foldcomp_db_reader: &FoldcompDbReader,
) -> (Vec<(Vec<ResidueMatch>, f32)>, Vec<(Vec<ResidueMatch>, f32)>, usize, f32) {
    let compact = foldcomp_db_reader.read_single_structure(path).expect("Error reading structure from foldcomp db");
    let compact = compact.to_compact();

    // let mut indices_found: Vec<Vec<(usize, usize)>> = Vec::new();
    // Iterate over query vector and retrieve indices
    // Parallel
    let query_set: HashSet<GeometricHash> = HashSet::from_iter(query_vector.clone());
    let query_symmetry_map = get_hash_symmetry_map(&query_set);

    let (index_set1, index_set2) = prefilter_amino_acid(&query_set, _hash_type, &compact);
    let aa_filter = CombinationVecIterator::new_from_btreesets(&index_set1, &index_set2);
    let (indices_found , candidate_pairs) = retrieve_with_prefilter(
        &compact, &query_set, aa_filter, _nbin_dist, _nbin_angle, multiple_bin,
        dist_cutoff, ca_distance_cutoff, aa_dist_map
    );

    let candidate_pair_map: HashMap<usize, BTreeSet<(usize, usize)>> = candidate_pairs.into_iter().fold(
        HashMap::new(), |mut map, (qi, pair)| {
            if !map.contains_key(&qi) {
                let mut set = BTreeSet::new();
                set.insert(pair);
                map.insert(qi, set);
            } else {
                map.get_mut(&qi).unwrap().insert(pair);
            }
            map
        }
    );
    
    // Make a graph and find connected components with the same node count
    // NOTE: Naive implementation to find matching components. Need to be improved to handle partial matches
    let graph = create_index_graph(&indices_found);
    let connected = connected_components_with_given_node_count(&graph, node_count);
    
    // Parallel
    let output: Vec<(Vec<ResidueMatch>, f32, Vec<ResidueMatch>, f32)>  = connected.par_iter().map(|component| {
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
        // Find mapping between query residues and retrieved residues
        let (query_indices, retrieved_indices) = map_query_and_retrieved_residues(
            &subgraph, query_map, node_count, &query_symmetry_map,
        );

        let mut query_indices_scanned: Vec<usize> = Vec::new();
        let mut retrieved_indices_scanned: Vec<usize> = Vec::new();
        // Sort component to match retrieved indices
        let mut res_vec: Vec<ResidueMatch> = Vec::new();
        let mut res_vec_from_hash: Vec<ResidueMatch> = Vec::new();
        let mut count_map: HashMap<usize, usize> = HashMap::new();
        all_query_indices.iter().for_each(|&i| {
            // If i is in query_indices, get the corresponding retrieved index
            count_map.clear();
            if query_indices.contains(&i) {
                let index = query_indices.iter().position(|&x| x == i).unwrap();
                let (chain, res_ind) = get_chain_and_res_ind(&compact, retrieved_indices[index]);
                res_vec_from_hash.push(Some((chain, res_ind)));
                if !retrieved_indices_scanned.contains(&retrieved_indices[index]) {
                    res_vec.push(Some((chain, res_ind)));
                    query_indices_scanned.push(i);
                    retrieved_indices_scanned.push(retrieved_indices[index]);
                } else {
                    // Substitute res_vec
                    let prev_index = retrieved_indices_scanned.iter().position(|&x| x == retrieved_indices[index]).unwrap();
                    res_vec[prev_index] = None;
                    res_vec.push(Some((chain, res_ind)));
                    // Delete previous indices in query_indices_scanned and retrieved_indices_scanned
                    query_indices_scanned.remove(prev_index);
                    retrieved_indices_scanned.remove(prev_index);
                    query_indices_scanned.push(i);
                    retrieved_indices_scanned.push(retrieved_indices[index]);
                }
                // res_vec.push(res_index_to_char(chain, res_ind));
                // query_indices_scanned.push(i);
                // retrieved_indices_scanned.push(retrieved_indices[index]);
            } else {
                res_vec_from_hash.push(None);
                if candidate_pair_map.contains_key(&i) {
                    let pairs = candidate_pair_map.get(&i).unwrap().clone();
                    for (j, k) in pairs {
                        // If retrieved_indices contains k, add j to mapping
                        if retrieved_indices.contains(&k) {
                            if !count_map.contains_key(&j) {
                                count_map.insert(j, 1);
                            } else {
                                let count = count_map.get_mut(&j).unwrap();
                                *count += 1;
                            }
                        }
                    }
                }
                if !count_map.is_empty() {
                    let max = count_map.iter().max_by(|a, b| a.1.cmp(b.1)).unwrap();
                    if *max.1 > 1 && !retrieved_indices_scanned.contains(max.0) {
                        let (chain, res_ind) = get_chain_and_res_ind(&compact, *max.0);
                        res_vec.push(Some((chain, res_ind)));
                        query_indices_scanned.push(i);
                        retrieved_indices_scanned.push(*max.0);
                    } else {
                        res_vec.push(None);
                    }
                } else {
                    res_vec.push(None);
                }
            }
        });

        let rmsd_from_hash = rmsd_for_matched(
            query_structure, &compact, &query_indices, &retrieved_indices
        );
        
        let rmsd = if res_vec == res_vec_from_hash {
            rmsd_from_hash
        } else {
            rmsd_for_matched(
                query_structure, &compact, &query_indices_scanned, &retrieved_indices_scanned
            )
        };
        
        (res_vec_from_hash, rmsd_from_hash, res_vec, rmsd)
    }).collect();
    // Split 
    let (result_from_hash, result): (Vec<(Vec<ResidueMatch>, f32)>, Vec<(Vec<ResidueMatch>, f32)>) = output.into_iter().map(|(a, b, c, d)| {
        ((a, b), (c, d))
    }).unzip();
    // In result, find the maximum matching node count and minimum RMSD with max match
    let mut max_matching_node_count = 0;
    let mut min_rmsd_with_max_match = 0.0;
    result.iter().for_each(|(res_vec, rmsd)| {
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




// Returns a vector of 1) chain+residue index as String and 2) RMSD value as f32
// 2025-01-08 10:51:23 
// Return a vector of ResidueMatch and RMSD values
pub fn retrieval_wrapper(
    path: &str, node_count: usize, query_vector: &Vec<GeometricHash>,
    _hash_type: HashType, _nbin_dist: usize, _nbin_angle: usize, 
    multiple_bin: &Option<Vec<(usize, usize)>>, dist_cutoff: f32,
    query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    query_structure: &CompactStructure, all_query_indices: &Vec<usize>,
    aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
    ca_distance_cutoff: f32,
) -> (Vec<(Vec<ResidueMatch>, f32)>, Vec<(Vec<ResidueMatch>, f32)>, usize, f32) {
    // Load structure to retrieve motif
    let compact = read_structure_from_path(&path).expect("Error reading structure from path");
    let compact = compact.to_compact();
    // let mut indices_found: Vec<Vec<(usize, usize)>> = Vec::new();
    // Iterate over query vector and retrieve indices
    // Parallel
    let query_set: HashSet<GeometricHash> = HashSet::from_iter(query_vector.clone());
    let query_symmetry_map = get_hash_symmetry_map(&query_set);
    
    let (index_set1, index_set2) = prefilter_amino_acid(&query_set, _hash_type, &compact);
    let aa_filter = CombinationVecIterator::new_from_btreesets(&index_set1, &index_set2);
    let (indices_found , candidate_pairs) = retrieve_with_prefilter(
        &compact, &query_set, aa_filter, _nbin_dist, _nbin_angle,
        multiple_bin, dist_cutoff, ca_distance_cutoff, aa_dist_map
    );

    let candidate_pair_map: HashMap<usize, BTreeSet<(usize, usize)>> = candidate_pairs.into_iter().fold(
        HashMap::new(), |mut map, (qi, pair)| {
            if !map.contains_key(&qi) {
                let mut set = BTreeSet::new();
                set.insert(pair);
                map.insert(qi, set);
            } else {
                map.get_mut(&qi).unwrap().insert(pair);
            }
            map
        }
    );
    
    // Make a graph and find connected components with the same node count
    // NOTE: Naive implementation to find matching components. Need to be improved to handle partial matches
    let graph = create_index_graph(&indices_found);
    let connected = connected_components_with_given_node_count(&graph, node_count);
    
    // Parallel
    let output: Vec<(Vec<ResidueMatch>, f32, Vec<ResidueMatch>, f32)>  = connected.par_iter().map(|component| {
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
        // Find mapping between query residues and retrieved residues
        let (query_indices, retrieved_indices) = map_query_and_retrieved_residues(
            &subgraph, query_map, node_count, &query_symmetry_map,
        );

        let mut query_indices_scanned: Vec<usize> = Vec::new();
        let mut retrieved_indices_scanned: Vec<usize> = Vec::new();
        // Sort component to match retrieved indices
        let mut res_vec: Vec<ResidueMatch> = Vec::new();
        let mut res_vec_from_hash: Vec<ResidueMatch> = Vec::new();
        let mut count_map: HashMap<usize, usize> = HashMap::new();
        all_query_indices.iter().for_each(|&i| {
            // If i is in query_indices, get the corresponding retrieved index
            count_map.clear();
            if query_indices.contains(&i) {
                let index = query_indices.iter().position(|&x| x == i).unwrap();
                let (chain, res_ind) = get_chain_and_res_ind(&compact, retrieved_indices[index]);
                res_vec_from_hash.push(Some((chain, res_ind)));
                if !retrieved_indices_scanned.contains(&retrieved_indices[index]) {
                    res_vec.push(Some((chain, res_ind)));
                    query_indices_scanned.push(i);
                    retrieved_indices_scanned.push(retrieved_indices[index]);
                } else {
                    // Substitute res_vec
                    let prev_index = retrieved_indices_scanned.iter().position(|&x| x == retrieved_indices[index]).unwrap();
                    res_vec[prev_index] = None;
                    res_vec.push(Some((chain, res_ind)));
                    // Delete previous indices in query_indices_scanned and retrieved_indices_scanned
                    query_indices_scanned.remove(prev_index);
                    retrieved_indices_scanned.remove(prev_index);
                    query_indices_scanned.push(i);
                    retrieved_indices_scanned.push(retrieved_indices[index]);
                }
            } else {
                res_vec_from_hash.push(None);
                if candidate_pair_map.contains_key(&i) {
                    let pairs = candidate_pair_map.get(&i).unwrap().clone();
                    // pairs.sort_by(|a, b| a.0.cmp(&b.0));
                    // pairs.dedup();
                    for (j, k) in pairs {
                        // If retrieved_indices contains k, add j to mapping
                        if retrieved_indices.contains(&k) {
                            if !count_map.contains_key(&j) {
                                count_map.insert(j, 1);
                            } else {
                                let count = count_map.get_mut(&j).unwrap();
                                *count += 1;
                            }
                        }
                    }
                }
                if !count_map.is_empty() {
                    let max = count_map.iter().max_by(|a, b| a.1.cmp(b.1)).unwrap();
                    if *max.1 > 1 && !retrieved_indices_scanned.contains(max.0) {
                        let (chain, res_ind) = get_chain_and_res_ind(&compact, *max.0);
                        res_vec.push(Some((chain, res_ind)));
                        query_indices_scanned.push(i);
                        retrieved_indices_scanned.push(*max.0);
                    } else {
                        res_vec.push(None);
                    }
                } else {
                    res_vec.push(None);
                }
            }
        });

        let rmsd_from_hash = rmsd_for_matched(
            query_structure, &compact, &query_indices, &retrieved_indices
        );
        
        let rmsd = if res_vec == res_vec_from_hash {
            rmsd_from_hash
        } else {
            rmsd_for_matched(
                query_structure, &compact, &query_indices_scanned, &retrieved_indices_scanned
            )
        };
        
        (res_vec_from_hash, rmsd_from_hash, res_vec, rmsd)
    }).collect();
    // Split 
    let (result_from_hash, result): (Vec<(Vec<ResidueMatch>, f32)>, Vec<(Vec<ResidueMatch>, f32)>) = output.into_iter().map(|(a, b, c, d)| {
        ((a, b), (c, d))
    }).unzip();
    // In result, find the maximum matching node count and minimum RMSD with max match
    let mut max_matching_node_count = 0;
    let mut min_rmsd_with_max_match = 0.0;
    result.iter().for_each(|(res_vec, rmsd)| {
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
    let mut query_symmetry_map = HashMap::with_capacity(query_set.len());
    query_set.iter().for_each(|hash| {
        let symmetry = hash.is_symmetric();
        query_symmetry_map.insert(hash.clone(), symmetry);
    });
    query_symmetry_map
}

pub fn prefilter_amino_acid(query_set: &HashSet<GeometricHash>, _hash_type: HashType, compact: &CompactStructure) -> (BTreeSet<usize>, BTreeSet<usize>) {
    let mut observed_aa1: HashSet<u8> = HashSet::with_capacity(20);
    let mut observed_aa2: HashSet<u8> = HashSet::with_capacity(20);
    let mut index_vec1 = BTreeSet::new();
    let mut index_vec2 = BTreeSet::new();
    let mut feature_holder = vec![0.0; 9];
    query_set.iter().for_each(|hash| {
        hash.reverse_hash_default(&mut feature_holder);
        let aa1 = feature_holder[_hash_type.amino_acid_index().unwrap()[0]] as u8;
        let aa2 = feature_holder[_hash_type.amino_acid_index().unwrap()[1]] as u8;
        if !observed_aa1.contains(&aa1) {
            observed_aa1.insert(aa1);
            compact.residue_name.iter().enumerate().filter_map(|(i, &res)| {
                if res == map_u8_to_aa(aa1).as_bytes() {
                    Some(i)
                } else {
                    None
                }
            }).for_each(|i| {
                index_vec1.insert(i);
            });
        }
        if !observed_aa2.contains(&aa2) {
            observed_aa2.insert(aa2);
            compact.residue_name.iter().enumerate().filter_map(|(i, &res)| {
                if res == map_u8_to_aa(aa2).as_bytes() {
                    Some(i)
                } else {
                    None
                }
            }).for_each(|i| {
                index_vec2.insert(i);
            });
        }
    });
    (index_vec1, index_vec2)
}

pub fn map_query_and_retrieved_residues(
    retrieved: &Graph<usize, GeometricHash>, 
    query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    node_count: usize,
    query_symmetry_map: &HashMap<GeometricHash, bool>,
) -> (Vec<usize>, Vec<usize>) {
    let mut query_indices: Vec<usize> = Vec::new();
    let mut retrieved_indices: Vec<usize> = Vec::new();
    // Iterate over edges in the graph
    retrieved.edge_indices().for_each(|edge| {
        let (i, j) = retrieved.edge_endpoints(edge).unwrap();
        let hash = retrieved[edge];
        let query = query_map.get(&hash);
        if query.is_some() {
            let (query_i, query_j) = query.unwrap().0;
            // If symmetry is true, push both query_i and query_j in ascending order
            if *query_symmetry_map.get(&hash).unwrap() {
                if !query_indices.contains(&query_i) && !retrieved_indices.contains(&retrieved[i]) &&
                   !query_indices.contains(&query_j) && !retrieved_indices.contains(&retrieved[j]) {
                    // Remap by ordering
                    if query_i < query_j {
                        if retrieved[i] < retrieved[j] {
                            query_indices.push(query_i);
                            retrieved_indices.push(retrieved[i]);
                            query_indices.push(query_j);
                            retrieved_indices.push(retrieved[j]);
                        } else {
                            query_indices.push(query_i);
                            retrieved_indices.push(retrieved[j]);
                            query_indices.push(query_j);
                            retrieved_indices.push(retrieved[i]);
                        }
                    } else {
                        if retrieved[i] < retrieved[j] {
                            query_indices.push(query_j);
                            retrieved_indices.push(retrieved[i]);
                            query_indices.push(query_i);
                            retrieved_indices.push(retrieved[j]);
                        } else {
                            query_indices.push(query_j);
                            retrieved_indices.push(retrieved[j]);
                            query_indices.push(query_i);
                            retrieved_indices.push(retrieved[i]);
                        }
                    }
                }
            }
            if !query_indices.contains(&query_i) && !retrieved_indices.contains(&retrieved[i]) {
                query_indices.push(query_i);
                retrieved_indices.push(retrieved[i]);
            }
            if !query_indices.contains(&query_j) && !retrieved_indices.contains(&retrieved[j]){
                query_indices.push(query_j);
                retrieved_indices.push(retrieved[j]);
            }
        }
        if retrieved_indices.len() == node_count {
            return;
        }
    });
    (query_indices, retrieved_indices)
}

pub fn rmsd_for_matched(
    compact1: &CompactStructure, compact2: &CompactStructure, 
    index1: &Vec<usize>, index2: &Vec<usize>
) -> f32 {
    let mut qcp = QCPSuperimposer::new();
    
    let coord_vec1: Vec<Coordinate> = index1.iter().map(
        |&i| (compact1.ca_vector.get_coord(i).unwrap(), compact1.cb_vector.get_coord(i).unwrap())
    ).flat_map(|(a, b)| vec![a, b]).collect();
    
    let coord_vec2: Vec<Coordinate> = index2.iter().map(
        |&i| (compact2.ca_vector.get_coord(i).unwrap(), compact2.cb_vector.get_coord(i).unwrap())
    ).flat_map(|(a, b)| vec![a, b]).collect();

    qcp.set_atoms(&coord_vec1, &coord_vec2);
    qcp.run();
    qcp.get_rms()
}

#[cfg(test)]
mod tests {
    use crate::controller::query::make_query_map;

    use super::*;

    #[test]
    fn test_retrieval_wrapper() {
        let path = String::from("data/serine_peptidases_filtered/4cha.pdb");
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
            &dist_thresholds, &angle_thresholds, &aa_substitutions, dist_cutoff, false
        );
        let queries: Vec<GeometricHash> = query_map.keys().cloned().collect();
        let compact = read_structure_from_path(&path).expect("Error reading structure from path");
        let compact = compact.to_compact();
        let new_path = String::from("data/serine_peptidases_filtered/4cha.pdb");
        let output = measure_time!(retrieval_wrapper(
            &new_path, query_residues.len(), &queries, hash_type, nbin_dist, nbin_angle, &None,
            dist_cutoff, &query_map, &compact, &query_indices, &aa_dist_map, 1.5,
        ));
        println!("{:?}", output);
    }

}