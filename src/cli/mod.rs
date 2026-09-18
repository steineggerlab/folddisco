//! Command line interface for FoldDisco

// File: mod.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved


pub mod workflows;
pub mod config;

/// Parsed arguments of each subcommand.
pub enum AppArgs {
    Global {
        help: bool,
    },
    Index {
        pdb_container: Option<String>,
        hash_type: String,
        index_path: String,
        num_threads: usize,
        num_bin_dist: usize,
        num_bin_angle: usize,
        multiple_bins: Option<String>,
        grid_width: f32,
        max_residue: usize,
        recursive: bool,
        mmap_on_disk: bool,
        id_type: String,
        // Build-time expansion; off unless expand_radius > 0 or aa_subst is set
        expand_radius: usize,
        expand_distance: f32,
        expand_angle: f32,
        aa_subst: Option<String>,
        verbose: bool,
        help: bool,
    },
    Query {
        pdb_path: String,
        query_string: String,
        threads: usize,
        index_path: Option<String>,
        skip_match: bool, // Changed from retrieve to skip_match. Now mathcing is default
        // Tolerances; None means not given, so an expanded index can skip them
        dist_threshold: Option<String>,
        angle_threshold: Option<String>,
        ca_dist_threshold: f32,
        expand_radius: Option<usize>, // feature dimensions allowed to deviate at once
        sensitive: bool,              // preset that raises the expansion radius
        confident: bool,              // preset that keeps only full, low-RMSD matches
        aa_subst: Option<String>,     // substitution scheme applied to every residue
        // Structure-level filters
        total_match_count: usize, 
        covered_node_count: usize,
        covered_node_ratio: f32,
        max_matching_node_count: usize,
        max_matching_node_ratio: f32,
        num_res_cutoff: usize,
        plddt_cutoff: f32,
        // Structure- and match-level filter
        idf_score_cutoff: f32,
        // Match-level filters
        connected_node_count: usize,
        connected_node_ratio: f32,
        rmsd_cutoff: f32,
        tm_score_cutoff: f32,
        gdt_ts_cutoff: f32,
        gdt_ha_cutoff: f32,
        chamfer_distance_cutoff: f32,
        hausdorff_distance_cutoff: f32,
        drmsd_cutoff: f32,
        top_n: usize,
        web_mode: bool,
        // Hash sampling
        sampling_count: Option<usize>,
        sampling_ratio: Option<f32>,
        freq_filter: Option<f32>,
        length_penalty: Option<f32>,
        sort_by: String,
        format_output: Option<String>,
        output_per_structure: bool,
        output_per_match: bool,
        output_with_superpose: bool,
        skip_ca_match: bool,
        partial_fit: bool, // Enable LMS based superposition.
        header: bool,
        serial_query: bool,
        // Always write CHAIN_RESIDUE
        chain_separator: bool,
        output: String,
        // One verdict row per query instead of the hit list
        novelty_mode: bool,
        novelty_coverage_threshold: f32,
        novelty_rmsd_threshold: f32,
        verbose: bool,
        help: bool,
    },
    Benchmark {
        // Required tabular files
        result: Option<String>,
        answer: Option<String>,
        // Optional; neutral hits are not counted as false positives
        neutral: Option<String>,
        index: Option<String>,
        input: Option<String>,
        format: String,
        fp: Option<f64>,
        threads: usize,
        afdb_to_uniprot: bool,
        // Column index per file [0]
        column_result: usize,
        column_answer: usize,
        column_neutral: usize,
        // Header line per file
        header_result: bool,
        header_answer: bool,
        header_neutral: bool,
    },
    Analyze {
        // Required 
        index_path: Option<String>,
        // Optional
        pdb_container: Option<String>,
        output: Option<String>,
        // Summary options
        top_n: usize,
        // enrichment options
        p_value: f64,
        min_support: usize,
        max_pos: usize,
        // other general options
        threads: usize,
        verbose: bool,
        help: bool,
    },
    Test {
        index_path: String,
        verbose: bool,
    },
}

/// Print the ASCII logo to stderr.
pub fn print_logo() {
    let logo = [
        "",
        "\x1b[91m░█▀▀░█▀█░█░░░█▀▄░\x1b[93m█▀▄░▀█▀░█▀▀░█▀▀░█▀█\x1b[0m",
        "\x1b[91m░█▀▀░█░█░█░░░█░█░\x1b[93m█░█░░█░░▀▀█░█░░░█░█\x1b[0m",
        "\x1b[91m░▀░░░▀▀▀░▀▀▀░▀▀░░\x1b[93m▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀\x1b[0m",
        "",
    ];

    for line in &logo {
        eprintln!("{}", line);
    }
}
