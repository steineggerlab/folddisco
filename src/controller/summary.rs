// File: summary.rs
// Created: 2025-12-27 00:54:13
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2025 Hyunbin Kim, All rights reserved

use crate::geometry::core::HashType;
use crate::index::indextable::FolddiscoIndex;
use crate::utils::log::{print_log_msg, log_msg};

use memmap2::Mmap;

pub struct EncodingStat {
    pub hash_type: HashType,
    pub dist_bin_given: usize,
    pub angle_bin_given: usize,
    pub total_encodings: usize,
    pub total_possible_encodings: usize,
    pub total_empty_encodings: usize,
    pub total_nonempty_encodings: usize,
    pub encoding_density: f32,
    pub aa_pair_counts: [[usize; 20]; 20],
    pub dist_bin_counts: Vec<usize>,
    pub angle_bin_counts: Vec<usize>,
}

impl EncodingStat {
    pub fn new(
        hash_type: HashType, dist_bin_given: usize, angle_bin_given: usize,
        index: &FolddiscoIndex,
    ) -> Self {
        let total_encodings = index.total_hashes;
        let aa_pair_counts = [[0usize; 20]; 20];
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
        }
    }

    pub fn empty() -> Self {
        Self {
            hash_type: HashType::PDBTrRosetta,
            dist_bin_given: 0,
            angle_bin_given: 0,
            total_encodings: 0,
            total_possible_encodings: 0,
            total_empty_encodings: 0,
            total_nonempty_encodings: 0,
            encoding_density: 0.0,
            aa_pair_counts: [[0usize; 20]; 20],
            dist_bin_counts: vec![],
            angle_bin_counts: vec![],
        }
    }

}

pub fn total_possible_encodings(hash_type: HashType, dist_bin_given: usize, angle_bin_given: usize) -> usize {
    hash_type.total_bins(dist_bin_given, angle_bin_given)
}

pub fn count_encodings(
    index: &FolddiscoIndex, offset_mmap: &Mmap,
    hash_type: HashType, nbin_dist: usize, nbin_angle: usize,
    verbose: bool
) -> EncodingStat {
    if verbose {
        print_log_msg("INFO", "Counting encodings in the index");
    }
    // todo!("Calculate total, empty, non-empty encodings, density, etc.");
    let mut stats = EncodingStat::new(
        hash_type,
        nbin_dist,
        nbin_angle,
        index,
    );

    println!("Total encodings: {}", stats.total_encodings);
    println!("Total possible encodings: {}", stats.total_possible_encodings);
    println!("Total empty encodings: {}", stats.total_empty_encodings);
    println!("Encoding density: {:.4}%", stats.encoding_density);
    
    // todo!("Count features: AA pairs, counts per bin (distance, angle)");
    
    // todo!("If doable, calculate secondary structure distribution of encodings by defining boundaries for secondary structures");
    
    // todo!("Top N frequent encodings. sort by frequency");
    
    stats
}

pub fn save_summary(encoding_stat: &EncodingStat, output_prefix: &str, verbose: bool) {
    todo!("Save encoding statistics to a summary file");
}

pub fn analyze_enrichment(
    index: &FolddiscoIndex, offset_mmap: &Mmap,
    pdb_container: &str,
    verbose: bool
) {
    if verbose {
        print_log_msg("INFO", &format!(
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

mod tests {
    use super::*;
    #[test]
    fn test_encoding_stat() {
        // 
        todo!("Test with real index");
        let stat = EncodingStat {
            hash_type: HashType::PDBTrRosetta,
            dist_bin_given: 16,
            angle_bin_given: 4,
            total_encodings: 1000,
            total_possible_encodings: 100000,
            total_empty_encodings: 99000,
            total_nonempty_encodings: 1000,
            encoding_density: 1.0,
            aa_pair_counts: [[0usize; 20]; 20],
            dist_bin_counts: vec![0usize; 16],
            angle_bin_counts: vec![0usize; 4],
        };
        println!("EncodingStat: total_encodings = {}", stat.total_encodings);
    }
}