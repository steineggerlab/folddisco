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
        pdb_path_vec: Vec<String>,
        hash_type: String,
        index_path: String,
        num_threads: usize,
        verbose: bool,
        help: bool,
    },
    Query {
        pdb_path: String,
        query_string: String,
        threads: usize,
        index_path: Option<String>,
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
        println!("{}", line);
    }
}
