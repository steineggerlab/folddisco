// File: loader.rs
// Created: 2024-02-29 21:54:20
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

const ALLOWED_EXTENSIONS: [&str; 2] = ["pdb", "ent"];


pub fn load_path(dir: &str, recursive: bool) -> Vec<String> {
    // Load all pdbs in given path
    let mut pdb_paths = Vec::new();
    let paths = std::fs::read_dir(dir).expect("Unable to read pdb directory");

    for path in paths {
        let path = path.expect("Unable to read path");
        let path = path.path();
        let path = path.to_str().expect("Unable to convert path to string");
        // 
        if recursive {
            if std::path::Path::new(path).is_dir() {
                let mut sub_pdb_paths = load_path(path, recursive);
                pdb_paths.append(&mut sub_pdb_paths);
            } else {
                // Check if allowed extension
                if ALLOWED_EXTENSIONS.iter().any(|&ext| path.ends_with(ext)) {
                    pdb_paths.push(path.to_string());
                }
            }
        } else {
            // Check if allowed extension
            if ALLOWED_EXTENSIONS.iter().any(|&ext| path.ends_with(ext)) {
                pdb_paths.push(path.to_string());
            }
        }
    }
    pdb_paths
}

pub fn get_all_combination(n: usize, include_same: bool) -> Vec<(usize, usize)> {
    let mut res = Vec::new();
    for i in 0..n {
        for j in 0..n {
            if i == j && !include_same {
                continue;
            }
            res.push((i, j));
        }
    }
    res
}

pub fn load_homeobox_toy() -> Vec<String> {
    vec![
        "data/homeobox/1akha-.pdb".to_string(),
        "data/homeobox/1b72a-.pdb".to_string(),
        "data/homeobox/1b72b-.pdb".to_string(),
        "data/homeobox/1ba5--.pdb".to_string(),
    ]
}

pub fn load_yeast_proteome() -> Vec<String> {
    // Load all pdbs in data/yeast
    let mut pdb_paths = Vec::new();
    let paths = std::fs::read_dir("data/yeast").expect("Unable to read yeast proteome");
    for path in paths {
        let path = path.expect("Unable to read path");
        let path = path.path();
        let path = path.to_str().expect("Unable to convert path to string");
        // If the path is a pdb file, add it to the list
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
        let pdb_paths = load_path("data/homeobox", false);
        assert_eq!(pdb_paths.len(), 4);
        println!("Flat: {:?}", pdb_paths);
        let pdb_paths = load_path("data/homeobox", true);
        assert_eq!(pdb_paths.len(), 6);
        println!("Recursive: {:?}", pdb_paths);
    }
    
    #[test]
    fn test_get_all_combination() {
        let combinations = get_all_combination(3, false);
        assert_eq!(combinations.len(), 6);
        println!("{:?}", combinations);
    }
}
