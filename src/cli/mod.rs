// File: mod.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

// Arguments of CLI app are defined here

pub mod workflows;

#[derive(Debug)]
enum Subcommand {
    Index,
    Query,
    // Add subcommands here
}

pub enum AppArgs {
    Global {
        help: bool,
    },
    Index {
        pdb_dir: Option<String>,
        hash_type: String,
        index_path: String,
        num_threads: usize,
        num_bin_dist: usize,
        num_bin_angle: usize,
        chunk_size: usize,
        max_residue: usize,
        recursive: bool,
        verbose: bool,
        help: bool,
    },
    Query {
        pdb_path: String,
        query_string: String,
        threads: usize,
        index_path: Option<String>,
        exact_match: bool,
        retrieve: bool,
        // Match thresholds
        dist_threshold: Option<String>,
        angle_threshold: Option<String>,
        // Cutoffs
        match_cutoff: f32,
        score_cutoff: f32,
        num_res_cutoff: usize,
        plddt_cutoff: f32,
        help: bool,
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
