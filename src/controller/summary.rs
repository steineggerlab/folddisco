// File: summary.rs
// Created: 2025-12-27 00:54:13
// Author: Hyunbin Kim (khb7840@gmail.com)
// Enrichment analysis idea from Alex Bott
// Copyright © 2025 Hyunbin Kim, All rights reserved

use crate::controller::DEFAULT_DIST_CUTOFF;
use crate::controller::Folddisco;
use crate::controller::feature::get_single_feature;
use crate::controller::io::read_structure_from_path;
use crate::geometry::core::HashType;
use crate::geometry::core::GeometricHash;
use crate::index::indextable::FolddiscoIndex;
use crate::utils::combination::CombinationIterator;
use crate::utils::convert::map_u8_to_aa;
use crate::utils::log::{log_msg, print_log_msg, INFO, FAIL};
use crate::utils::loader::load_path;
use crate::measure_time;

use std::fs::File;
use std::io::Write;
use std::sync::atomic::AtomicUsize;

use dashmap::DashMap;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator as _};
use rayon::slice::ParallelSliceMut;

pub struct AAPairCounts {
    pub counts: Vec<Vec<AtomicUsize>>,
}
unsafe impl Send for AAPairCounts {}
unsafe impl Sync for AAPairCounts {}

pub struct EncodingStat {
    pub hash_type: HashType,
    pub dist_bin_given: usize,
    pub angle_bin_given: usize,
    pub total_encodings: usize,
    pub total_possible_encodings: usize,
    pub total_empty_encodings: usize,
    pub total_nonempty_encodings: usize,
    pub encoding_density: f32,
    pub aa_pair_counts: AAPairCounts,
    pub dist_bin_counts: Vec<usize>,
    pub angle_bin_counts: Vec<usize>,
    pub hash_count_vec: Vec<(u32, usize)>,
}

impl EncodingStat {
    pub fn new(
        hash_type: HashType, dist_bin_given: usize, angle_bin_given: usize,
        index: &FolddiscoIndex,
    ) -> Self {
        let total_encodings = index.total_hashes;
        let aa_pair_counts = AAPairCounts {
            counts: (0..20)
                .map(|_| (0..20).map(|_| AtomicUsize::new(0)).collect())
                .collect(),
        };
        let total_dist_bins = match hash_type.dist_bins(dist_bin_given) {
            Some(b) => b,
            None => 1,
        };
        let total_angle_bins = match hash_type.angle_bins(angle_bin_given) {
            Some(b) => b,
            None => 1,
        };
        let total_possible = total_possible_encodings(
            hash_type,
            dist_bin_given,
            angle_bin_given,
        );
        let empty_encodings = total_possible - total_encodings;
        let density = (total_encodings as f64 / total_possible as f64 * 100.0) as f32;

        let dist_bin_counts = vec![0usize; total_dist_bins];
        let angle_bin_counts = vec![0usize; total_angle_bins];
        
        let hash_count_vec = get_hash_count_vec(index);
        
        Self {
            hash_type,
            dist_bin_given,
            angle_bin_given,
            total_encodings,
            total_possible_encodings: total_possible,
            total_empty_encodings: empty_encodings,
            total_nonempty_encodings: total_encodings,
            encoding_density: density,
            aa_pair_counts,
            dist_bin_counts,
            angle_bin_counts,
            hash_count_vec,
        }
    }

    pub fn empty() -> Self {
        let aa_pair_counts = AAPairCounts {
            counts: (0..20)
                .map(|_| (0..20).map(|_| AtomicUsize::new(0)).collect())
                .collect(),
        };
        Self {
            hash_type: HashType::PDBTrRosetta,
            dist_bin_given: 0,
            angle_bin_given: 0,
            total_encodings: 0,
            total_possible_encodings: 0,
            total_empty_encodings: 0,
            total_nonempty_encodings: 0,
            encoding_density: 0.0,
            aa_pair_counts: aa_pair_counts,
            dist_bin_counts: vec![],
            angle_bin_counts: vec![],
            hash_count_vec: vec![],
        }
    }
}

pub fn total_possible_encodings(hash_type: HashType, dist_bin_given: usize, angle_bin_given: usize) -> usize {
    hash_type.total_bins(dist_bin_given, angle_bin_given)
}

pub fn count_encodings(
    index: &FolddiscoIndex, hash_type: HashType, nbin_dist: usize, nbin_angle: usize, verbose: bool
) -> EncodingStat {
    if verbose {
        print_log_msg(INFO, "Counting encodings in the index");
    }
    let nbin_dist = match nbin_dist {
        0 => hash_type.default_dist_bin(),
        _ => nbin_dist,
    };
    let nbin_angle = match nbin_angle {
        0 => hash_type.default_angle_bin(),
        _ => nbin_angle,
    };
    // todo!("Calculate total, empty, non-empty encodings, density, etc.");
    let mut stats = EncodingStat::new(
        hash_type,
        nbin_dist,
        nbin_angle,
        index,
    );

    stats.hash_count_vec.par_sort_by(|a, b| b.1.cmp(&a.1)); // Descending order by count
    
    if verbose {
        print_log_msg(INFO, &format!("Hash count vector loaded: {} entries", stats.hash_count_vec.len()));
    }
    
    // todo!("Count features: AA pairs, counts per bin (distance, angle)");
    // DONE: naively count AA pairs from hash_count_vec
    // TODO: Count dist_bin_counts and angle_bin_counts
    stats.hash_count_vec.par_iter().for_each(|(hash, count)| {
        let mut feature_container = vec![0.0f32; 7];
        let geometric_hash = GeometricHash::from_u32(*hash, hash_type);
        // Decode features
        GeometricHash::reverse_hash(
            &geometric_hash, nbin_dist, nbin_angle, &mut feature_container
        );
        let aa1 = feature_container[0] as usize;
        let aa2 = feature_container[1] as usize;
        // Update AA pair counts
        stats.aa_pair_counts.counts[aa1][aa2].fetch_add(
            *count, std::sync::atomic::Ordering::Relaxed
        );
    });
    stats
}

pub fn save_summary(
    encoding_stat: &EncodingStat, output_prefix: &str, 
    top_n: usize, verbose: bool
) -> std::io::Result<()> {
    use std::fs::File;
    use std::io::Write;

    if verbose {
        print_log_msg(INFO, &format!("Saving summary to files with prefix: {}", output_prefix));
    }

    // 1. TSV of overall statistics
    let stats_file = format!("{}_stats.tsv", output_prefix);
    let mut f = File::create(&stats_file)?;
    writeln!(f, "metric\tvalue")?;
    writeln!(f, "total\t{}", encoding_stat.total_encodings)?;
    writeln!(f, "possible\t{}", encoding_stat.total_possible_encodings)?;
    writeln!(f, "empty\t{}", encoding_stat.total_empty_encodings)?;
    writeln!(f, "nonempty\t{}", encoding_stat.total_nonempty_encodings)?;
    writeln!(f, "density\t{:.4}", encoding_stat.encoding_density)?;
    if verbose {
        print_log_msg(INFO, &format!("Saved overall statistics to {}", stats_file));
    }

    // 2. TSV of top N encodings with decoded features
    let top_n_file = format!("{}_top{}.tsv", output_prefix, top_n);
    let mut f = File::create(&top_n_file)?;
    writeln!(f, "rank\thash\tcount\taa1\taa2\tca_dist\tcb_dist\tca_cb_angle\tphi1\tphi2")?;
    
    let nbin_dist = encoding_stat.dist_bin_given;
    let nbin_angle = encoding_stat.angle_bin_given;
    let hash_type = encoding_stat.hash_type;
    let mut feature_container = vec![0.0f32; 7];
    
    for (i, (hash, count)) in encoding_stat.hash_count_vec.iter().take(top_n).enumerate() {
        let geometric_hash = GeometricHash::from_u32(*hash, hash_type);
        GeometricHash::reverse_hash(
            &geometric_hash, nbin_dist, nbin_angle, &mut feature_container
        );
        writeln!(
            f, "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            i + 1, hash, count,
            map_u8_to_aa(feature_container[0] as u8), map_u8_to_aa(feature_container[1] as u8), 
            feature_container[2], feature_container[3], feature_container[4], 
            feature_container[5], feature_container[6]
        )?;
    }
    if verbose {
        print_log_msg(INFO, &format!("Saved top {} encodings to {}", top_n, top_n_file));
    }

    // 3. CSV of amino acid pair counts
    let aa_pairs_file = format!("{}_aa_pairs.csv", output_prefix);
    let mut f = File::create(&aa_pairs_file)?;
    write!(f, "aa1_aa2")?;
    for aa2 in 0..20 {
        write!(f, ",{}", map_u8_to_aa(aa2))?;
    }
    writeln!(f)?;
    for aa1 in 0..20 {
        write!(f, "{}", map_u8_to_aa(aa1))?;
        for aa2 in 0..20 {
            let count = encoding_stat.aa_pair_counts.counts[aa1 as usize][aa2 as usize].load(
                std::sync::atomic::Ordering::Relaxed
            );
            write!(f, ",{}", count)?;
        }
        writeln!(f)?;
    }
    if verbose {
        print_log_msg(INFO, &format!("Saved amino acid pair counts to {}", aa_pairs_file));
    }

    // 4. TSV of count distribution with borders
    let count_dist_file = format!("{}_count_distribution.tsv", output_prefix);
    let mut f = File::create(&count_dist_file)?;
    writeln!(f, "frequency\tcount")?;
    
    let counts = get_counts_from_hash_count_vec(&encoding_stat.hash_count_vec, None);
    for (border, count) in counts.iter() {
        writeln!(f, "{}\t{}", border, count)?;
    }
    if verbose {
        print_log_msg(INFO, &format!("Saved count distribution to {}", count_dist_file));
    }

    Ok(())
}

pub fn analyze_enrichment(
    index: &FolddiscoIndex, pdb_container: &str,
    hash_type: HashType, nbin_dist: usize, nbin_angle: usize, 
    threads: usize, p_value_threshold: f64, output_prefix: &str, 
    count_cutoff: usize, max_residue_count: usize, verbose: bool
) -> std::io::Result<()> {
    if verbose {
        print_log_msg(INFO, &format!(
            "Analyzing enrichment of encodings in PDB container {}",
            pdb_container
        ));
    }
    // Test if the default hashing schemes are working
    let pdb_paths = load_path(pdb_container, false);
    let mut query_handler = Folddisco::new(
        pdb_paths, hash_type, threads,
        nbin_dist, nbin_angle, String::new(),
        DEFAULT_DIST_CUTOFF, None, false,
    );
    measure_time!(query_handler.collect_hash_vec(), verbose);
    
    let mut hashes = query_handler.hash_id_vec.clone();
    hashes.par_sort_unstable();
    // Count ids per hash. Hashes is sorted.
    let mut hash_counts = hashes.iter().fold(
        Vec::new(),
        |mut acc: Vec<(u32, usize)>, &hash| {
            if let Some(last) = acc.last_mut() {
                if last.0 == hash.0 {
                    last.1 += 1;
                } else {
                    acc.push((hash.0, 1));
                }
            } else {
                acc.push((hash.0, 1));
            }
            acc
        }
    );
    // Sort by hash.
    hash_counts.par_sort_by(|a, b| a.0.cmp(&b.0));

    let hash_id_pos_map = measure_time!(query_handler.collect_hash_id_pos(), verbose);
    
    // Load index and get hash_count_vec
    let mut bg_hash_count_vec = get_hash_count_vec(index);
    // Sort by hash.
    bg_hash_count_vec.par_sort_by(|a, b| a.0.cmp(&b.0));
    
    let enriched_hashes = get_enriched_hashes(
        &hash_counts, &bg_hash_count_vec, p_value_threshold
    );

    if verbose {
        print_log_msg(INFO, &format!(
            "Found {} enriched encodings with p-value < {}",
            enriched_hashes.len(), p_value_threshold
        ));
    }

    // Determine actual bin sizes
    let nbin_dist = match nbin_dist {
        0 => hash_type.default_dist_bin(),
        _ => nbin_dist,
    };
    let nbin_angle = match nbin_angle {
        0 => hash_type.default_angle_bin(),
        _ => nbin_angle,
    };

    // Build auxiliary maps for position and query summaries using DashMap
    let pos_hash_map: DashMap<(usize, String), Vec<u32>> = DashMap::new();
    let pdb_positions: DashMap<usize, Vec<String>> = DashMap::new();

    // 1. Write enriched hashes TSV
    let hash_file = format!("{}_enriched_hashes.tsv", output_prefix);
    let mut f = File::create(&hash_file)?;
    writeln!(f, "hash\tp_value\tfeatures\tpositions")?;
    
    for (hash, p_value) in enriched_hashes.iter() {
        // Decode features
        let mut feature_container = vec![0.0f32; 7];
        let geometric_hash = GeometricHash::from_u32(*hash, hash_type);
        GeometricHash::reverse_hash(
            &geometric_hash, nbin_dist, nbin_angle, &mut feature_container
        );
        
        let features_str = format!(
            "{},{},{:.4},{:.4},{:.4},{:.4},{:.4}",
            map_u8_to_aa(feature_container[0] as u8),
            map_u8_to_aa(feature_container[1] as u8),
            feature_container[2],
            feature_container[3],
            feature_container[4],
            feature_container[5],
            feature_container[6],
        );

        // Collect positions
        let mut positions_vec = Vec::new();
        if let Some(entry) = hash_id_pos_map.get(hash) {
            for (pdb_pos, pos1, pos2) in entry.iter() {
                let pos_str = format!("{}-{}-{}", 
                    query_handler.path_vec[*pdb_pos], pos1, pos2);
                positions_vec.push(pos_str);
                
                // Build position map - split pairs into individual positions
                let key1 = (*pdb_pos, pos1.clone());
                pos_hash_map.entry(key1).or_insert_with(Vec::new).push(*hash);
                let key2 = (*pdb_pos, pos2.clone());
                pos_hash_map.entry(key2).or_insert_with(Vec::new).push(*hash);
                
                // Build PDB position map - split pairs into individual positions
                pdb_positions.entry(*pdb_pos).or_insert_with(Vec::new).push(pos1.clone());
                pdb_positions.entry(*pdb_pos).or_insert_with(Vec::new).push(pos2.clone());
            }
        }
        
        let positions_str = positions_vec.join(",");
        writeln!(f, "{}\t{:.4e}\t{}\t{}", hash, p_value, features_str, positions_str)?;
    }
    
    if verbose {
        print_log_msg(INFO, &format!("Saved enriched hashes to {}", hash_file));
    }

    // 2. Write enriched positions TSV
    let pos_file = format!("{}_enriched_positions.tsv", output_prefix);
    let mut f = File::create(&pos_file)?;
    writeln!(f, "id\tpos\tcount\thash_list")?;
    
    // Convert DashMap to Vec for sorting
    let mut pos_list: Vec<((usize, String), Vec<u32>)> = pos_hash_map.into_iter().collect();
    // Sort by pdb id first, then by count (descending)
    pos_list.par_sort_by(|a, b| {
        a.0.0.cmp(&b.0.0)
            .then_with(|| b.1.len().cmp(&a.1.len()))
    });
    
    for ((pdb_pos, pos), hash_list) in pos_list.iter() {
        let count = hash_list.len();
        // Apply count cutoff
        if count <= count_cutoff {
            continue;
        }
        let id = &query_handler.path_vec[*pdb_pos];
        let hash_list_str = hash_list.iter()
            .map(|h| h.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(f, "{}\t{}\t{}\t{}", id, pos, count, hash_list_str)?;
    }
    
    if verbose {
        print_log_msg(INFO, &format!("Saved enriched positions to {}", pos_file));
    }

    // 3. Write query summary TSV
    let query_file = format!("{}_query_summary.tsv", output_prefix);
    let mut f = File::create(&query_file)?;
    writeln!(f, "pdb_path\tpositions")?;
    
    // Convert DashMap to Vec for sorting
    let mut pdb_list: Vec<(usize, Vec<String>)> = pdb_positions.into_iter().collect();
    pdb_list.par_sort_by_key(|(pdb_pos, _)| *pdb_pos);
    
    // Build position lookup map from pos_list
    let pos_count_map: DashMap<(usize, String), usize> = DashMap::new();
    pos_list.par_iter().for_each(|((pdb_pos, pos), hash_list)| {
        pos_count_map.insert((*pdb_pos, pos.clone()), hash_list.len());
    });
    
    for (pdb_pos, positions) in pdb_list.iter() {
        let pdb_path = &query_handler.path_vec[*pdb_pos];
        // Deduplicate positions and apply count cutoff
        let mut unique_positions = positions.clone();
        unique_positions.sort();
        unique_positions.dedup();
        // Filter positions by count cutoff and collect with their counts
        let mut filtered_with_counts: Vec<(String, usize)> = unique_positions.into_iter()
            .filter_map(|pos| {
                pos_count_map.get(&(*pdb_pos, pos.clone()))
                    .map(|count| (pos, *count))
                    .filter(|(_, count)| *count > count_cutoff)
            })
            .collect();
        
        // Sort by count (descending) and apply max residue count limit
        filtered_with_counts.par_sort_by(|a, b| b.1.cmp(&a.1));
        if filtered_with_counts.len() > max_residue_count {
            filtered_with_counts.truncate(max_residue_count);
        }
        
        // Sort final positions by chain id (alphabetically) and residue index (numerically)
        let mut final_positions: Vec<String> = filtered_with_counts.iter()
            .map(|(pos, _)| pos.clone())
            .collect();
        final_positions.sort_by(|a, b| {
            // Extract chain id (first char) and residue index (remaining chars)
            let a_chain = a.chars().next().unwrap_or(' ');
            let b_chain = b.chars().next().unwrap_or(' ');
            let a_idx: i32 = a.chars().skip(1).collect::<String>().parse().unwrap_or(0);
            let b_idx: i32 = b.chars().skip(1).collect::<String>().parse().unwrap_or(0);
            
            a_chain.cmp(&b_chain).then_with(|| a_idx.cmp(&b_idx))
        });
        
        if !final_positions.is_empty() {
            let positions_str = final_positions.join(",");
            writeln!(f, "{}\t{}", pdb_path, positions_str)?;
        }
    }
    
    if verbose {
        print_log_msg(INFO, &format!("Saved query summary to {}", query_file));
    }

    Ok(())
}

// pub fn save_enrichment_result(
//     // analysis_result: &EnrichmentResult,
//     output_prefix: &str,
//     verbose: bool
// ) {
//     todo!("Save enrichment analysis result to output files");
// }

pub fn get_hash_count_vec(
    index: &FolddiscoIndex,
) -> Vec<(u32, usize)> {
    let mut hash_count_vec: Vec<(u32, usize)> = Vec::with_capacity(index.total_hashes);
    for i in 0..index.total_hashes {
        let hash = index.loaded_hashes[i];
        let count = index.loaded_offsets[i + 1] - index.loaded_offsets[i];
        hash_count_vec.push((hash, count));
    }
    hash_count_vec
}

pub fn get_counts_from_hash_count_vec(
    hash_count_vec: &Vec<(u32, usize)>, borders: Option<&Vec<usize>>,
) -> Vec<(usize, usize)> {
    // hash_count vector is sorted by count in descending order
    let max_count = hash_count_vec.first().unwrap().1;
    let borders = match borders {
        Some(b) => b.clone(),
        None => {
            // Starting from 1, include power of twos until max_count
            let mut b: Vec<usize> = vec![];
            let mut next_base = 1;
            // Generate borders like 10, 100, 1000, ,,,
            while next_base <= max_count {
                b.push(next_base);
                next_base *= 2;
            }
            // Add max_count if not already included
            if *b.last().unwrap() < max_count {
                b.push(max_count);
            }
            b
        }
    };
    let mut counts: Vec<(usize, usize)> = Vec::with_capacity(borders.len());
    for (i, &border) in borders.iter().enumerate() {
        let count = if i == 0 {
            // Using partition_point to speed up
            let upper_bound = hash_count_vec.len();
            let lower_bound = hash_count_vec.partition_point(|&(_, c)| c > border);
            upper_bound - lower_bound
        } else {
            let lower_border = borders[i - 1];
            let upper_idx = hash_count_vec.partition_point(|&(_, c)| c > border);
            let lower_idx = hash_count_vec.partition_point(|&(_, c)| c > lower_border);
            lower_idx - upper_idx
        };
        counts.push((border, count));
    }
    counts
}

fn get_enriched_hashes(
    query_hash_counts: &Vec<(u32, usize)>, bg_hash_counts: &Vec<(u32, usize)>, 
    p_value_cutoff: f64
) -> Vec<(u32, f64)> {
    let total_query: usize = query_hash_counts.iter().map(|(_, c)| c).sum();
    let total_bg: usize = bg_hash_counts.iter().map(|(_, c)| c).sum();
    
    // Parallel processing for enrichment testing
    let mut enriched: Vec<(u32, f64)> = query_hash_counts.par_iter()
        .filter_map(|(hash, query_count)| {
            // Binary search for the hash in background
            let bg_count = match bg_hash_counts.binary_search_by_key(hash, |(h, _)| *h) {
                Ok(idx) => bg_hash_counts[idx].1,
                Err(_) => 0, // Hash not found in background
            };
            
            // Hypergeometric test parameters
            let n = total_query; // Sample size (query)
            let k = bg_count + query_count; // Total successes in population
            let n_total = total_bg + total_query; // Population size
            let x = *query_count; // Observed successes in sample
            
            // Calculate p-value (right-tail test for enrichment)
            let p_value = hypergeometric_test(x, n, k, n_total);
            
            if p_value < p_value_cutoff {
                Some((*hash, p_value))
            } else {
                None
            }
        })
        .collect();
    
    // Sort enriched by p-value ascending
    enriched.par_sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    enriched
}

/// Hypergeometric test p-value (right-tail)
/// P(X >= x) where X ~ Hypergeometric(n_total, k, n)
fn hypergeometric_test(x: usize, n: usize, k: usize, n_total: usize) -> f64 {
    let mut p_value = 0.0;
    let max_x = n.min(k);
    
    for i in x..=max_x {
        p_value += hypergeometric_pmf(i, n, k, n_total);
    }
    
    p_value.min(1.0) // Clamp to [0, 1]
}

/// Hypergeometric probability mass function using log-space to avoid overflow
fn hypergeometric_pmf(x: usize, n: usize, k: usize, n_total: usize) -> f64 {
    // P(X = x) = C(k, x) * C(N-k, n-x) / C(N, n)
    let log_prob = log_binomial(k, x) 
        + log_binomial(n_total - k, n - x) 
        - log_binomial(n_total, n);
    
    log_prob.exp()
}

/// Log of binomial coefficient: log(C(n, k))
fn log_binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    if k == 0 || k == n {
        return 0.0;
    }
    
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

/// Log factorial using Stirling's approximation for large n
fn log_factorial(n: usize) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    if n < 20 {
        // Exact for small n
        (2..=n).map(|i| (i as f64).ln()).sum()
    } else {
        // Stirling's approximation: ln(n!) ≈ n*ln(n) - n + 0.5*ln(2πn)
        let n = n as f64;
        n * n.ln() - n + 0.5 * (2.0 * std::f64::consts::PI * n).ln()
    }
}

impl Folddisco {
    pub fn collect_hash_id_pos(&self) -> DashMap<u32, Vec<(usize, String, String)>> {
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating files");
        // For iterating files, apply multi-threading with num_threads_for_file

        let output_map: DashMap<u32, Vec<(usize, String, String)>> = DashMap::new();
        
        pool.install(|| {
            // Preserve locality for multi-threading
            self.path_vec
                .par_iter()
                .enumerate()
                .for_each(|(pdb_pos, pdb_path)| {
                    // #[cfg(not(feature = "foldcomp"))]
                    let compact = read_structure_from_path(pdb_path).expect(
                        log_msg(FAIL, "Failed to read structure").as_str()
                    );

                    // #[cfg(feature = "foldcomp")]
                    // let compact = if self.is_foldcomp_enabled {
                    //     self.foldcomp_db_reader.read_single_structure_by_id(
                    //         self.numeric_db_key_vec[pdb_pos]
                    //     ).expect(
                    //         log_msg(FAIL, "Failed to read structure").as_str()
                    //     )
                    // } else {
                    //     read_structure_from_path(pdb_path).expect(
                    //         log_msg(FAIL, "Failed to read structure").as_str()
                    //     )
                    // };  

                    let compact = compact.to_compact();
                    let res_bound = CombinationIterator::new(compact.num_residues);
                    let mut feature = vec![0.0; 9];
                    res_bound.for_each(|(i, j)| {
                        if i == j {
                            return;
                        }
                        let has_feature = get_single_feature(
                            i, j, &compact, self.hash_type, self.dist_cutoff, &mut feature
                        );
                        if has_feature {
                            if self.num_bin_dist == 0 || self.num_bin_angle == 0 {
                                let hash = GeometricHash::perfect_hash_default_as_u32(&feature, self.hash_type);
                                let pos1 = format!("{}{}", compact.chain_per_residue[i] as char, compact.residue_serial[i]);
                                let pos2 = format!("{}{}", compact.chain_per_residue[j] as char, compact.residue_serial[j]);
                                output_map.entry(hash).or_default().push((pdb_pos, pos1, pos2));
                            } else {
                                let hash = GeometricHash::perfect_hash_as_u32(
                                    &feature, self.hash_type, self.num_bin_dist, self.num_bin_angle
                                );
                                let pos1 = format!("{}{}", compact.chain_per_residue[i] as char, compact.residue_serial[i]);
                                let pos2 = format!("{}{}", compact.chain_per_residue[j] as char, compact.residue_serial[j]);
                                output_map.entry(hash).or_default().push((pdb_pos, pos1, pos2));
                            }
                        }
                    });
                    // Drop intermediate variables
                    drop(compact);
                })
        });
        drop(pool);
        output_map
    }
}


// TODO: unit tests with a real index