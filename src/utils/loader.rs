// File: loader.rs
// Created: 2024-02-29 21:54:20
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

const ALLOWED_EXTENSIONS: [&str; 6] = [".pdb", ".ent", ".cif", ".pdb.gz", ".ent.gz", ".cif.gz"];


/// Structure files (`.pdb`, `.ent`, `.cif`, optionally gzipped) in `dir`.
pub fn load_path(dir: &str, recursive: bool) -> Vec<String> {
    let mut pdb_paths = Vec::new();
    let paths = std::fs::read_dir(dir).expect("Unable to read pdb directory");

    for path in paths {
        let path = path.expect("Unable to read path");
        let path = path.path();
        let path = path.to_str().expect("Unable to convert path to string");
        if recursive {
            if std::path::Path::new(path).is_dir() {
                let mut sub_pdb_paths = load_path(path, recursive);
                pdb_paths.append(&mut sub_pdb_paths);
            } else {
                if ALLOWED_EXTENSIONS.iter().any(|&ext| path.ends_with(ext)) {
                    pdb_paths.push(path.to_string());
                }
            }
        } else {
            if ALLOWED_EXTENSIONS.iter().any(|&ext| path.ends_with(ext)) {
                pdb_paths.push(path.to_string());
            }
        }
    }
    pdb_paths
}

/// Fixed homeobox test set under `data/homeobox`.
pub fn load_homeobox_toy() -> Vec<String> {
    vec![
        "data/homeobox/1akha-.pdb".to_string(),
        "data/homeobox/1b72a-.pdb".to_string(),
        "data/homeobox/1b72b-.pdb".to_string(),
        "data/homeobox/1ba5--.pdb".to_string(),
    ]
}

/// `.pdb` files under `data/yeast`.
pub fn load_yeast_proteome() -> Vec<String> {
    let mut pdb_paths = Vec::new();
    let paths = std::fs::read_dir("data/yeast").expect("Unable to read yeast proteome");
    for path in paths {
        let path = path.expect("Unable to read path");
        let path = path.path();
        let path = path.to_str().expect("Unable to convert path to string");
        if path.ends_with(".pdb") {
            pdb_paths.push(path.to_string());
        }
    }
    pdb_paths
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_path() {
        let pdb_paths = load_path("data/io_test", false);
        assert_eq!(pdb_paths.len(), 5);
        println!("Flat: {:?}", pdb_paths);
        let pdb_paths = load_path("data/io_test", true);
        // 14 structures, plus data/io_test/cif/1G2F_multichain.cif
        assert_eq!(pdb_paths.len(), 15);
        println!("Recursive: {:?}", pdb_paths);
    }
}
