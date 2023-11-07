// File: e_coli.test.rs
// Created: 2023-11-06 21:01:32
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
use std::collections::HashMap;
// Import unique
use std::collections::HashSet;

use motifsearch::controller::{self, Controller, GeometryHashCollector};
use motifsearch::geometry::trrosetta_subfamily::{HashCollection, HashValue};
// use motifsearch::geometry::trrosetta::{HashCollection, HashValue};
use motifsearch::index::builder::IndexBuilder;
use motifsearch::index::query_multiple_with_neighbors;
use motifsearch::index::{IndexTablePrinter, query_single, query_multiple};
use motifsearch::PDBReader;
use motifsearch::index::index_table;
use motifsearch::utils::loader::{load_path, get_all_combination};
use motifsearch::utils::benchmark::{Metrics, calculate_metrics, compare_target_answer};

use rayon::prelude::*;

pub fn main() {
    rayon::ThreadPoolBuilder::new().num_threads(6).build_global().unwrap();
    // Load all PDB files in the given path
    let pdb_paths: Vec<String> = load_path("analysis/raw_ecoli");

    let start = std::time::Instant::now();
    let mut controller = Controller::new(pdb_paths);
    controller.fill_numeric_id_vec();
    let mut index_table = index_table::IndexTable::new();

    for i in 0.. controller.path_vec.len() {
        let pdb_path = &controller.path_vec[i];
        let pdb_reader = PDBReader::from_file(pdb_path).expect("PDB file not found");
        let compact = pdb_reader.read_structure().expect("Failed to read PDB file");
        let compact = compact.to_compact();
        let mut hash_vec: Vec<u64> = Vec::with_capacity(compact.num_residues.pow(2));
        let mut res_pair: Vec<(u16, u16)> = Vec::with_capacity(compact.num_residues.pow(2));
        
        // Generate a combination of n and m
        let res_bound = get_all_combination(compact.num_residues, false);
        let r: Vec<(u64, (u16, u16))> = (&res_bound).into_par_iter().map(
            |(n, m)| {
                if n == m {
                    return (0u64, (0u16, 0u16))
                }
                let trr = compact.get_trrosetta_feature(*n, *m).unwrap_or([0.0; 6]);
                if trr[0] < 2.0 || trr[0] > 20.0 {
                    return (0u64, (0u16, 0u16))
                }
                let hash_value = HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);
                // hash_vec.push(hash_value.as_u64());
                
                let res1 = compact.residue_serial[*n] as u16;
                let res2 = compact.residue_serial[*m] as u16;
                // res_pair.push((res1, res2));
                (hash_value.as_u64(), (res1, res2))
                // println!("n: {:?} m: {:?} hash: {:?} res1: {:?} res2: {:?}", n, m, hash_value.as_u64(), res1, res2);
            }
        ).collect();
        // Filter out the 0s
        let r = r.into_iter().filter(|(x, _)| *x != 0).collect::<Vec<_>>();
        hash_vec = r.iter().map(|(x, _)| *x).collect::<Vec<_>>();
        res_pair = r.iter().map(|(_, x)| *x).collect::<Vec<_>>();
        index_table.concat_structure(i, hash_vec);
    }
    let lap1 = std::time::Instant::now();
    println!("[NEW] Time elapsed for building index table {:?}", lap1 - start);

    // Measure time
    // let homeobox_queries = vec![
    //     HashValue::from_u64(2203335001092u64),
    //     HashValue::from_u64(3298551727108u64),
    //     HashValue::from_u64(3298551661828u64),
    // ];
    let homeobox_queries = vec![
        HashValue::from_u64(5501903505417u64),
        HashValue::from_u64(7696631925000u64),
        HashValue::from_u64(6597120101128u64),
    ];
    let homeobox_queries_as_u64_ref = homeobox_queries.iter().map(|x| x.as_u64()).collect::<Vec<u64>>();
    // let homeobox_neighboring_queries: Vec<Vec<u64>> = homeobox_queries.iter().map(|x| x.neighbors(true).iter().map(|x| x.as_u64()).collect()).collect();
    // let result = index_table.query_multiple_with_neighbors(&homeobox_neighboring_queries);
    // let mut target2 = Vec::new();
    // // Remove duplicates in result
    // if let Some(result) = result {
    //     let mut set = HashSet::new();
    //     let mut dedup_result = Vec::new();
    //     for i in result {
    //         if !set.contains(&i) {
    //             set.insert(i);
    //             dedup_result.push(i);
    //         }
    //     }
    //     for i in dedup_result {
    //         // println!("[NEW] Queried prot (3) {:?}", controller.path_vec.get(i.get_id()));
    //         let str_id = controller.path_vec.get(i.get_id()).unwrap();
    //         target2.push(str_id.clone() as String);
    //     }
    // }
    let lap3 = std::time::Instant::now();
    // println!("[NEW] Time elapsed for quering neighbors {:?}", lap3 - lap1);
    let result = index_table.query_multiple(&homeobox_queries_as_u64_ref);

    let mut target = Vec::new();
    if let Some(result) = result {
        let mut set = HashSet::new();
        let mut dedup_result = Vec::new();
        for i in result {
            if !set.contains(&i) {
                set.insert(i);
                dedup_result.push(i);
            }
        }

        for i in dedup_result {
            println!("[NEW] Queried prot (3) {:?}", controller.path_vec.get(i.get_id()));
            let str_id = controller.path_vec.get(i.get_id()).unwrap();
            target.push(str_id.clone() as String);
        }
    }
    let end = std::time::Instant::now();
    println!("[NEW] Time elapsed for quering {:?}", end - lap3);

    // // Save
    // let result = index_table.save_to_bin("analysis/test_index_table.bin");
    // if let Err(e) = result {
    //     println!("Error saving index table to binary: {:?}", e);
    // }
    let answer = vec![
        // "analysis/raw_ecoli/AF-P76176-F1-model_v4.pdb", "analysis/raw_ecoli/AF-P0C0V0-F1-model_v4.pdb", "analysis/raw_ecoli/AF-P39099-F1-model_v4.pdb", 
        // "analysis/raw_ecoli/AF-P0AEE3-F1-model_v4.pdb", 
        "analysis/raw_ecoli/AF-Q7BSW5-F1-model_v4.pdb", "analysis/raw_ecoli/AF-Q47692-F1-model_v4.pdb", 
        "analysis/raw_ecoli/AF-O68900-F1-model_v4.pdb", "analysis/raw_ecoli/AF-Q7BS42-F1-model_v4.pdb", "analysis/raw_ecoli/AF-Q84GK0-F1-model_v4.pdb", 
        "analysis/raw_ecoli/AF-Q9EZE7-F1-model_v4.pdb", "analysis/raw_ecoli/AF-Q8FDW4-F1-model_v4.pdb",
    ];
    let answer = answer.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    println!("Target: {:?}", target);
    // Calculate metrics
    let metrics = compare_target_answer::<String>(&target, &answer, &controller.path_vec);
    // Print TP, TN, FP, FN
    println!("TP: {}", metrics.true_pos);
    println!("TN: {}", metrics.true_neg);
    println!("FP: {}", metrics.false_pos);
    println!("FN: {}", metrics.false_neg);
    
    println!("Precision: {}", metrics.precision());
    println!("Recall: {}", metrics.recall());
    println!("Accuracy: {}", metrics.accuracy());
    println!("F1 score: {}", metrics.f1_score());
    // let metrics = compare_target_answer::<String>(&target2, &answer, &controller.path_vec);
    // println!("Target2: {:?}", target2);
    // // Print TP, TN, FP, FN
    // println!("TP: {}", metrics.true_pos);
    // println!("TN: {}", metrics.true_neg);
    // println!("FP: {}", metrics.false_pos);
    // println!("FN: {}", metrics.false_neg);
    
    // println!("Precision: {}", metrics.precision());
    // println!("Recall: {}", metrics.recall());
    // println!("Accuracy: {}", metrics.accuracy());
    // println!("F1 score: {}", metrics.f1_score());
    
}


