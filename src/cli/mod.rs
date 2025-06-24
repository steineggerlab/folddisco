//! Command line interface for FoldDisco

// File: mod.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

// Arguments of CLI app are defined here

pub mod workflows;
pub mod config;

pub enum AppArgs {
    Global {
        help: bool,
    },
    Index {
        pdb_container: Option<String>,
        hash_type: String,
        index_path: String,
        mode: String,
        num_threads: usize,
        num_bin_dist: usize,
        num_bin_angle: usize,
        multiple_bins: Option<String>,
        grid_width: f32,
        chunk_size: usize,
        max_residue: usize,
        recursive: bool,
        mmap_on_disk: bool,
        id_type: String,
        verbose: bool,
        help: bool,
    },
    Query {
        pdb_path: String,
        query_string: String,
        threads: usize,
        index_path: Option<String>,
        skip_match: bool, // Changed from retrieve to skip_match. Now mathcing is default
        // Match thresholds
        dist_threshold: String,
        angle_threshold: String,
        ca_dist_threshold: f32,
        // filtering parameters
        // These are for filtering StructQueryResult only
        total_match_count: usize, 
        covered_node_count: usize,
        covered_node_ratio: f32,
        covered_edge_count: usize,
        covered_edge_ratio: f32,
        max_matching_node_count: usize,
        max_matching_node_ratio: f32,
        num_res_cutoff: usize,
        plddt_cutoff: f32,
        // These are for filtering both StructQueryResult and MatchQueryResult
        idf_score_cutoff: f32,
        // These are for filtering MatchQueryResult only
        connected_node_count: usize,
        connected_node_ratio: f32,
        rmsd_cutoff: f32,
        // top N filtering
        top_n: usize,
        web_mode: bool,
        //.Query sampling
        sampling_count: Option<usize>,
        sampling_ratio: Option<f32>,
        freq_filter: Option<f32>,
        length_penalty: Option<f32>,
        // sorting mode
        sort_by_rmsd: bool,
        sort_by_score: bool,
        // output mode
        output_per_structure: bool,
        output_per_match: bool,
        output_with_superpose: bool,
        skip_ca_match: bool,
        header: bool,
        serial_query: bool,
        output: String,
        verbose: bool,
        help: bool,
    },
    Benchmark {
        // Required tabular files
        result: Option<String>,
        answer: Option<String>,
        // Optional tabular files: Neutral list is not considered as false positive
        neutral: Option<String>,
        index: Option<String>,
        input: Option<String>,
        format: String,
        fp: Option<f64>,
        threads: usize,
        afdb_to_uniprot: bool,
        // Column index for each file. Default is 0
        column_result: usize,
        column_answer: usize,
        column_neutral: usize,
        // Use header for each file. Default is false
        header_result: bool,
        header_answer: bool,
        header_neutral: bool,
    },
    Test {
        index_path: String,
        verbose: bool,
    },
}

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
