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
        pdb_dir: Option<String>,
        hash_type: String,
        index_path: String,
        mode: String,
        num_threads: usize,
        num_bin_dist: usize,
        num_bin_angle: usize,
        grid_width: f32,
        chunk_size: usize,
        max_residue: usize,
        recursive: bool,
        id_type: String,
        verbose: bool,
        help: bool,
    },
    Query {
        pdb_path: String,
        query_string: String,
        threads: usize,
        index_path: Option<String>,
        retrieve: bool,
        // Match thresholds
        amino_acid: u8,
        dist_threshold: Option<String>,
        angle_threshold: Option<String>,
        // Cutoffs
        match_cutoff: Option<String>,
        score_cutoff: f32,
        num_res_cutoff: usize,
        plddt_cutoff: f32,
        verbose: bool,
        help: bool,
    },
    Benchmark {
        result: Option<String>,
        answer: Option<String>,
        index: Option<String>,
        format: String,
    },
    Test {
        index_path: String,
        verbose: bool,
    }
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
