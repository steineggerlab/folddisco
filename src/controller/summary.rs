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
    pub num_bin_dist: usize,
    pub num_bin_angle: usize,
    pub total_encodings: usize,
    
}

impl EncodingStat {
    pub fn new(
        hash_type: HashType, num_bin_dist: usize, num_bin_angle: usize,
        index: &FolddiscoIndex,
    ) -> Self {
        let total_encodings = index.total_hashes;
        Self {
            hash_type,
            num_bin_dist,
            num_bin_angle,
            total_encodings,
        }
    }

    pub fn empty() -> Self {
        Self {
            hash_type: HashType::PDBTrRosetta,
            num_bin_dist: 0,
            num_bin_angle: 0,
            total_encodings: 0,
        }
    }
    
    fn total_possible_encodings(&self) -> usize {
        todo!("Recalclate with num_bin_dist and num_bin_angle & hash_type");
        2_usize.pow(self.hash_type.encoding_bits() as u32)
    }

    fn total_empty_encodings(&self) -> usize {
        self.total_possible_encodings() - self.total_encodings
    }
    
    fn encoding_density(&self, percentage: bool) -> f32 {
        let density = self.total_encodings as f32 / self.total_possible_encodings() as f32;
        if percentage {
            density * 100.0
        } else {
            density
        }
    }
}

pub fn count_encodings(
    index: &FolddiscoIndex, offset_mmap: &Mmap,
    hash_type: HashType, nbin_dist: usize, nbin_angle: usize,
    verbose: bool
) -> EncodingStat {
    if verbose {
        print_log_msg("INFO", "Counting encodings in the index");
    }

    // 
    todo!("Calculate total, empty, non-empty encodings, density, etc.");
    
    todo!("Count features: AA pairs, counts per bin (distance, angle)");
    
    todo!("If doable, calculate secondary structure distribution of encodings by defining boundaries for secondary structures");
    
    todo!("Top N frequent encodings. sort by frequency");
    
    
    EncodingStat::new(
        hash_type,
        nbin_dist,
        nbin_angle,
        index,
    )
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