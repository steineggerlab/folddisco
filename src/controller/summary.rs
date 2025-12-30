// File: summary.rs
// Created: 2025-12-27 00:54:13
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2025 Hyunbin Kim, All rights reserved

use crate::geometry::core::HashType;
use crate::index::indextable::FolddiscoIndex;
use crate::geometry::core::GeometricHash;
use crate::utils::convert::map_u8_to_aa;
use crate::utils::log::{print_log_msg, INFO};

use std::sync::atomic::AtomicUsize;

use memmap2::Mmap;
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
    index: &FolddiscoIndex, offset_mmap: &Mmap,
    pdb_container: &str,
    verbose: bool
) {
    if verbose {
        print_log_msg(INFO, &format!(
            "Analyzing enrichment of encodings in PDB container {}",
            pdb_container
        ));
    }
    
    // Load PDBs from container
    todo!("Load PDBs from the given container (directory or compressed file)");
    
    // For each PDB, compute encoding distribution
    todo!("For each PDB structure, compute its encoding distribution");
    
    // Compare with index encoding distribution
    todo!("Compare the encoding distribution of PDBs with that of the index");
    
    // Calculate enrichment statistics
    todo!("Calculate enrichment statistics for encodings present in PDBs vs index");
}

pub fn save_enrichment_result(
    // analysis_result: &EnrichmentResult,
    output_prefix: &str,
    verbose: bool
) {
    todo!("Save enrichment analysis result to output files");
}

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


mod tests {
    use super::*;
    #[test]
    fn test_encoding_stat() {
        // 
        todo!("Test with real index");
        let aa_pair_counts = AAPairCounts {
            counts: (0..20)
                .map(|_| (0..20).map(|_| AtomicUsize::new(0)).collect())
                .collect(),
        };
        let stat = EncodingStat {
            hash_type: HashType::PDBTrRosetta,
            dist_bin_given: 16,
            angle_bin_given: 4,
            total_encodings: 1000,
            total_possible_encodings: 100000,
            total_empty_encodings: 99000,
            total_nonempty_encodings: 1000,
            encoding_density: 1.0,
            aa_pair_counts: aa_pair_counts,
            dist_bin_counts: vec![0usize; 16],
            angle_bin_counts: vec![0usize; 4],
            hash_count_vec: vec![],
        };
        println!("EncodingStat: total_encodings = {}", stat.total_encodings);
    }
}