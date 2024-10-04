// 
use std::{collections::{BTreeSet, HashMap, HashSet}, fs::File};
use petgraph::Graph;
use rayon::iter::{IntoParallelRefIterator, ParallelDrainFull, ParallelIterator};

use crate::utils::convert::{map_aa_to_u8, map_u8_to_aa}; 
use crate::prelude::*; 
use crate::structure::{coordinate::Coordinate, core::CompactStructure, qcp::QCPSuperimposer}; 
use crate::utils::combination::{CombinationIterator, CombinationVecIterator};
use crate::controller::graph::{connected_components_with_given_node_count, create_index_graph};
use crate::controller::feature::get_single_feature;

#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;

use super::feature;

pub fn hash_vec_to_aa_pairs(hash_vec: &Vec<GeometricHash>) -> HashSet<(u32, u32)> {
    let mut output: HashSet<(u32, u32)> = HashSet::new();
    let mut feature = vec![0.0; 9];
    for hash in hash_vec {
        hash.reverse_hash_default(&mut feature);
        output.insert((feature[0] as u32, feature[1] as u32));
    }
    output
}

pub fn retrieve_residue_with_hash(
    hash_set: &HashSet<GeometricHash>, _aa_filter: &HashSet<(u32, u32)>, path: &str, hash_type: HashType, nbin_dist: usize, nbin_angle: usize, dist_cutoff: f32,
) -> Option<Vec<((u8, u8), (u64, u64))>> {

    let file = File::open(path).expect("File not found");
    let pdb = PDBReader::new(file);
    let compact = pdb.read_structure().expect("Error reading structure");
    let compact = compact.to_compact();
    let comb = CombinationIterator::new(compact.num_residues);
    let mut output: Vec<((u8, u8), (u64, u64))> = Vec::new();
    let mut feature = vec![0.0; 9];
    comb.for_each(|(i, j)| {
        // let aa_pair = (map_aa_to_u8(&compact.residue_name[i]) as u32, map_aa_to_u8(&compact.residue_name[j]) as u32);
        // if !aa_filter.contains(&aa_pair) {
        //     return;
        // }
        let is_feature = get_single_feature(i, j, &compact, hash_type, dist_cutoff, &mut feature);
        if is_feature {
            let curr_hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(&feature, hash_type)
            } else {
                GeometricHash::perfect_hash(&feature, hash_type, nbin_dist, nbin_angle)
            };
            if hash_set.get(&curr_hash).is_some() {
                output.push((
                    (compact.chain_per_residue[i], compact.chain_per_residue[j]),
                    (compact.residue_serial[i], compact.residue_serial[j])
                ));
            }
        }
    });

    if output.is_empty() {
        None
    } else {
        Some(output)
    }
}

pub fn connected(res_ind_vec: &Vec<((u8, u8), (u64, u64))>, _len: usize) -> usize {
    let mut res_set: HashMap<String, usize> = HashMap::new();
    for (i, j) in res_ind_vec {
        let key1 = format!("{}{}", i.0 as char, j.0);
        let key2 = format!("{}{}", i.1 as char, j.1);
        if !res_set.contains_key(&key1) {
            res_set.insert(key1, 1);
        } else {
            let count = res_set.get_mut(&key1).unwrap();
            *count += 1;
        }
        if !res_set.contains_key(&key2) {
            res_set.insert(key2, 1);
        } else {
            let count = res_set.get_mut(&key2).unwrap();
            *count += 1;
        }
    }
    // Count elements that are greater than 1
    let max = res_set.iter().filter(|(_, &v)| v > 1).count();
    max
}

pub fn cycle(_res_ind_vec: &Vec<((u8, u8), (u64, u64))>) -> usize {
    // Use res_ind_vec as a directed graph and count cycles
    todo!();
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
    compact: &CompactStructure,
    // hash: &GeometricHash, 
    hash_set: &HashSet<GeometricHash>,
    // prefilter: &Vec<(usize, usize)>, 
    prefilter: CombinationVecIterator,
    nbin_dist: usize, nbin_angle: usize, dist_cutoff: f32,
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
                            if (curr_dist - dist).abs() < 1.0 {
                                candidate_pairs.push((*qi, (i, j)));
                            }
                        }
                    }
                }
                let curr_hash = if nbin_dist == 0 || nbin_angle == 0 {
                    GeometricHash::perfect_hash_default(&feature, hash.hash_type())
                } else {
                    GeometricHash::perfect_hash(&feature, hash.hash_type(), nbin_dist, nbin_angle)
                };
                if hash_set.contains(&curr_hash) {
                    output.push((i, j, curr_hash));
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
                            if (curr_dist - dist).abs() < 1.5 {
                                candidate_pairs.push((*qi, (i, j)));
                            }
                        }
                    }
                }

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
    (output, candidate_pairs)
}


pub fn retrieve_whole_structure(
    compact: &CompactStructure,
    hash_set: &HashSet<GeometricHash>,
    nbin_dist: usize, nbin_angle: usize, dist_cutoff: f32,
    query_aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
) -> (Vec<(usize, usize, GeometricHash)>, Vec<(usize, (usize, usize))>) {
    let mut output: Vec<(usize, usize, GeometricHash)> = Vec::new();
    let mut candidate_pairs: Vec<(usize, (usize, usize))> = Vec::new();
    let hash = hash_set.iter().next().cloned().unwrap();
    let mut feature = vec![0.0; 9];
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
                        if (curr_dist - dist).abs() < 1.0 {
                            candidate_pairs.push((*qi, (i, j)));
                        }
                    }
                }
            }
            let curr_hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(&feature, hash.hash_type())
            } else {
                GeometricHash::perfect_hash(&feature, hash.hash_type(), nbin_dist, nbin_angle)
            };
            if hash_set.contains(&curr_hash) {
                output.push((i, j, curr_hash));
            }
        }
    });
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
    _hash_type: HashType, _nbin_dist: usize, _nbin_angle: usize, dist_cutoff: f32,
    query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    query_structure: &CompactStructure, all_query_indices: &Vec<usize>,
    aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
    foldcomp_db_reader: &FoldcompDbReader,
) -> (Vec<(String, f32)>, Vec<(String, f32)>) {
    
    let compact = foldcomp_db_reader.read_single_structure(path).expect("Error reading structure from foldcomp db");
    let compact = compact.to_compact();

    // let mut indices_found: Vec<Vec<(usize, usize)>> = Vec::new();
    // Iterate over query vector and retrieve indices
    // Parallel
    let query_set: HashSet<GeometricHash> = HashSet::from_iter(query_vector.clone());
    let (index_set1, index_set2) = prefilter_amino_acid(&query_set, _hash_type, &compact);
    let aa_filter = CombinationVecIterator::new_from_btreesets(&index_set1, &index_set2);
    let (indices_found , candidate_pairs) = retrieve_with_prefilter(
        &compact, &query_set, aa_filter, _nbin_dist, _nbin_angle, dist_cutoff, aa_dist_map
    );
    
    // Convert candidate_pairs to hashmap
    // let candidate_pair_map: HashMap<usize, Vec<(usize, usize)>> = candidate_pairs.into_iter().fold(
    //     HashMap::new(), |mut map, (qi, pair)| {
    //         if !map.contains_key(&qi) {
    //             map.insert(qi, vec![pair]);
    //         } else {
    //             map.get_mut(&qi).unwrap().push(pair);
    //         }
    //         map
    //     }
    // );
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
    // let indices_found = measure_time!(query_vector.par_iter().map(|hash| {
    //     let prefiltered = prefilter_aa_pair(&compact, hash);
    //     retrieve_with_prefilter(&compact, hash, prefiltered, _nbin_dist, _nbin_angle)
    // }).collect::<Vec<Vec<(usize, usize)>>>());
    
    // Make a graph and find connected components with the same node count
    // NOTE: Naive implementation to find matching components. Need to be improved to handle partial matches
    let graph = create_index_graph(&indices_found);
    let connected = connected_components_with_given_node_count(&graph, node_count);
    
    // Parallel
    let output: Vec<(String, f32, String, f32)>  = connected.par_iter().map(|component| {
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
            &subgraph, query_map, node_count,
        );

        let mut query_indices_scanned: Vec<usize> = Vec::new();
        let mut retrieved_indices_scanned: Vec<usize> = Vec::new();
        // Sort component to match retrieved indices
        let mut res_vec: Vec<String> = Vec::new();
        let mut res_vec_from_hash: Vec<String> = Vec::new();
        let mut count_map: HashMap<usize, usize> = HashMap::new();
        all_query_indices.iter().for_each(|&i| {
            // If i is in query_indices, get the corresponding retrieved index
            count_map.clear();
            if query_indices.contains(&i) {
                let index = query_indices.iter().position(|&x| x == i).unwrap();
                let (chain, res_ind) = get_chain_and_res_ind(&compact, retrieved_indices[index]);
                res_vec_from_hash.push(res_index_to_char(chain, res_ind));
                if !retrieved_indices_scanned.contains(&retrieved_indices[index]) {
                    res_vec.push(res_index_to_char(chain, res_ind));
                    query_indices_scanned.push(i);
                    retrieved_indices_scanned.push(retrieved_indices[index]);
                } else {
                    // Substitute res_vec
                    let prev_index = retrieved_indices_scanned.iter().position(|&x| x == retrieved_indices[index]).unwrap();
                    res_vec[prev_index] = "_".to_string();
                    res_vec.push(res_index_to_char(chain, res_ind));
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
                res_vec_from_hash.push("_".to_string());
                if candidate_pair_map.contains_key(&i) {
                    let mut pairs = candidate_pair_map.get(&i).unwrap().clone();
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
                        res_vec.push(res_index_to_char(chain, res_ind));
                        query_indices_scanned.push(i);
                        retrieved_indices_scanned.push(*max.0);
                    } else {
                        res_vec.push("_".to_string());
                    }
                } else {
                    res_vec.push("_".to_string());
                }
            }
        });
        // retrieved_indices.iter().enumerate().for_each(|(i, &node)| {
        //     let (chain, res_ind) = get_chain_and_res_ind(&compact, node);
        //     res_vec.push(res_index_to_char(chain, res_ind));
        // });

        let res_string_from_hash = res_vec_from_hash.join(",");
        let rmsd_from_hash = rmsd_for_matched(
            query_structure, &compact, &query_indices, &retrieved_indices
        );
        
        let res_string = res_vec.join(",");
        let rmsd = if res_string == res_string_from_hash {
            rmsd_from_hash
        } else {
            rmsd_for_matched(
                query_structure, &compact, &query_indices_scanned, &retrieved_indices_scanned
            )
        };
        
        (res_string_from_hash, rmsd_from_hash, res_string, rmsd)
    }).collect();
    // Split Vec<(String, f32, String, f32)> into Vec<(String, f32)> and Vec<(String, f32)>
    let (result_from_hash, result): (Vec<(String, f32)>, Vec<(String, f32)>) = output.into_iter().map(|(a, b, c, d)| {
        ((a, b), (c, d))
    }).unzip();
    (result_from_hash, result)
}




// Returns a vector of 1) chain+residue index as String and 2) RMSD value as f32
pub fn retrieval_wrapper(
    path: &str, node_count: usize, query_vector: &Vec<GeometricHash>,
    _hash_type: HashType, _nbin_dist: usize, _nbin_angle: usize, dist_cutoff: f32,
    query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    query_structure: &CompactStructure, all_query_indices: &Vec<usize>,
    aa_dist_map: &HashMap<(u8, u8), Vec<(f32, usize)>>,
) -> (Vec<(String, f32)>, Vec<(String, f32)>) {
    // Load structure to retrieve motif
    let pdb_loaded = PDBReader::new(File::open(&path).expect("File not found"));
    let compact = pdb_loaded.read_structure().expect("Error reading structure");
    let compact = compact.to_compact();
    // let mut indices_found: Vec<Vec<(usize, usize)>> = Vec::new();
    // Iterate over query vector and retrieve indices
    // Parallel
    let query_set: HashSet<GeometricHash> = HashSet::from_iter(query_vector.clone());

    let (index_set1, index_set2) = prefilter_amino_acid(&query_set, _hash_type, &compact);
    let aa_filter = CombinationVecIterator::new_from_btreesets(&index_set1, &index_set2);
    let (indices_found , candidate_pairs) = retrieve_with_prefilter(
        &compact, &query_set, aa_filter, _nbin_dist, _nbin_angle, dist_cutoff, aa_dist_map
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
    let output: Vec<(String, f32, String, f32)>  = connected.par_iter().map(|component| {
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
            &subgraph, query_map, node_count,
        );

        let mut query_indices_scanned: Vec<usize> = Vec::new();
        let mut retrieved_indices_scanned: Vec<usize> = Vec::new();
        // Sort component to match retrieved indices
        let mut res_vec: Vec<String> = Vec::new();
        let mut res_vec_from_hash: Vec<String> = Vec::new();
        let mut count_map: HashMap<usize, usize> = HashMap::new();
        all_query_indices.iter().for_each(|&i| {
            // If i is in query_indices, get the corresponding retrieved index
            count_map.clear();
            if query_indices.contains(&i) {
                let index = query_indices.iter().position(|&x| x == i).unwrap();
                let (chain, res_ind) = get_chain_and_res_ind(&compact, retrieved_indices[index]);
                res_vec_from_hash.push(res_index_to_char(chain, res_ind));
                if !retrieved_indices_scanned.contains(&retrieved_indices[index]) {
                    res_vec.push(res_index_to_char(chain, res_ind));
                    query_indices_scanned.push(i);
                    retrieved_indices_scanned.push(retrieved_indices[index]);
                } else {
                    // Substitute res_vec
                    let prev_index = retrieved_indices_scanned.iter().position(|&x| x == retrieved_indices[index]).unwrap();
                    res_vec[prev_index] = "_".to_string();
                    res_vec.push(res_index_to_char(chain, res_ind));
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
                res_vec_from_hash.push("_".to_string());
                if candidate_pair_map.contains_key(&i) {
                    let mut pairs = candidate_pair_map.get(&i).unwrap().clone();
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
                        res_vec.push(res_index_to_char(chain, res_ind));
                        query_indices_scanned.push(i);
                        retrieved_indices_scanned.push(*max.0);
                    } else {
                        res_vec.push("_".to_string());
                    }
                } else {
                    res_vec.push("_".to_string());
                }
            }
        });
        // retrieved_indices.iter().enumerate().for_each(|(i, &node)| {
        //     let (chain, res_ind) = get_chain_and_res_ind(&compact, node);
        //     res_vec.push(res_index_to_char(chain, res_ind));
        // });

        let res_string_from_hash = res_vec_from_hash.join(",");
        let rmsd_from_hash = rmsd_for_matched(
            query_structure, &compact, &query_indices, &retrieved_indices
        );
        
        let res_string = res_vec.join(",");
        let rmsd = if res_string == res_string_from_hash {
            rmsd_from_hash
        } else {
            rmsd_for_matched(
                query_structure, &compact, &query_indices_scanned, &retrieved_indices_scanned
            )
        };
        
        (res_string_from_hash, rmsd_from_hash, res_string, rmsd)
    }).collect();
    // Split Vec<(String, f32, String, f32)> into Vec<(String, f32)> and Vec<(String, f32)>
    let (result_from_hash, result): (Vec<(String, f32)>, Vec<(String, f32)>) = output.into_iter().map(|(a, b, c, d)| {
        ((a, b), (c, d))
    }).unzip();
    (result_from_hash, result)
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
    node_count: usize
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
            // 
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
    // Sort both indices 
    // let mut zipped: Vec<(usize, usize)> = query_indices.iter().zip(retrieved_indices.iter()).map(|(&a, &b)| (a, b)).collect();
    // // Sort by query indices and next by retrieved indices
    // zipped.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    // query_indices = zipped.iter().map(|(a, _)| *a).collect(); ㄴ
    // retrieved_indices = zipped.iter().map(|(_, b)| *b).collect();   
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
        let exact_match = false;
        let dist_thresholds: Vec<f32> = vec![0.5,1.0];
        let angle_thresholds: Vec<f32> = vec![5.0,10.0];
        let dist_cutoff = 20.0;
        let queries = make_query(
            &path, &query_residues, hash_type, nbin_dist, nbin_angle, exact_match, dist_thresholds.clone(), angle_thresholds.clone(), &aa_substitutions, dist_cutoff
        );
        let (query_map, query_indices, aa_dist_map ) = make_query_map(
            &path, &query_residues, hash_type, nbin_dist, nbin_angle, &dist_thresholds, &angle_thresholds, &aa_substitutions, dist_cutoff
        );
        let pdb_loaded = PDBReader::new(File::open(&path).expect("File not found"));
        let compact = pdb_loaded.read_structure().expect("Error reading structure");
        let compact = compact.to_compact();
        let new_path = String::from("data/serine_peptidases_filtered/1azw.pdb");
        let output = measure_time!(retrieval_wrapper(
            &new_path, query_residues.len(), &queries, hash_type, nbin_dist, nbin_angle, dist_cutoff, &query_map, &compact, &query_indices, &aa_dist_map
        ));
        println!("{:?}", output);
    }
    
    
    #[test]
    fn test_retrieve_residue_with_hash() {
        let path = String::from("data/serine_peptidases_filtered/4cha.pdb");
        let query_string = "B57,B102,C195";
        // let path = String::from("analysis/1g91.pdb");
        // let query_string = "A30,A32,A35";
        let (query_residues, aa_substitions) = parse_query_string(query_string, b'A');
        let hash_type = HashType::PDBMotifSinCos;
        let nbin_dist = 16;
        let nbin_angle = 4;
        let exact_match = false;
        let dist_thresholds: Vec<f32> = vec![0.5,1.0];
        let angle_thresholds: Vec<f32> = vec![5.0, 10.0];
        let dist_cutoff = 20.0;
        let queries = make_query(
            &path, &query_residues, hash_type, nbin_dist, nbin_angle, exact_match, dist_thresholds, angle_thresholds, &aa_substitions, dist_cutoff
        );
        let hash_set: HashSet<GeometricHash> = queries.iter().cloned().collect();
    
        let pdb_loaded = PDBReader::new(File::open(&path).expect("File not found"));
        let compact = pdb_loaded.read_structure().expect("Error reading structure");
        let compact = compact.to_compact();
        // print num residues
        println!("Num residues: {}", compact.num_residues);
        let aa_filter = hash_vec_to_aa_pairs(&queries);
        println!("{:?}", aa_filter);
        let retrieved = measure_time!(retrieve_residue_with_hash(
            &hash_set, &aa_filter, &path, hash_type, nbin_dist, nbin_angle, dist_cutoff
        ));
        println!("RETRIEVED: {:?}", retrieved);
        let connected = connected(&retrieved.unwrap(), query_residues.len());
        println!("CONNECTED_COMPONENTS: {:?}", connected);
    }
}