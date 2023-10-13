use std::collections::HashMap;
// Import unique
use std::collections::HashSet;

use motifsearch::controller::{self, Controller, GeometryHashCollector};
use motifsearch::geometry::trrosetta::{HashCollection, HashValue};
// use motifsearch::geometry::trrosetta_subfamily::{HashCollection, HashValue};
// use motifsearch::geometry::aa_pair::{HashCollection, HashValue};
use motifsearch::index::builder::IndexBuilder;
use motifsearch::index::query_multiple_with_neighbors;
use motifsearch::index::{IndexTablePrinter, query_single, query_multiple};
use motifsearch::PDBReader;

mod common;
use common::loader;

// #[test]
// fn test_geometry_hash_collector() {
//     let pdb_paths = loader::load_homeobox_toy();
//     for pdb_path in pdb_paths {
//         // Start measure time
//         let start = std::time::Instant::now();
//         let pdb_reader = PDBReader::from_file(&pdb_path).expect("Failed to read PDB file");
//         let structure = pdb_reader
//             .read_structure()
//             .expect("Failed to read structure");
//         let compact = structure.to_compact();

//         let mut hash_collector = GeometryHashCollector::new();
//         for i in 0..compact.num_residues {
//             for j in 0..compact.num_residues {
//                 if i == j {
//                     continue;
//                 }
//                 let trr = compact
//                     .get_trrosetta_feature(i, j)
//                     .expect("Failed to get trrosetta feature");
//                 let hash_value =
//                     HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);
//                 hash_collector.collect_hash(hash_value);
//             }
//         }
//         //
//         let before_dedup = hash_collector.hash_collection.len();
//         hash_collector.remove_redundancy();
//         let end = std::time::Instant::now();
//         println!(
//             "{:?} | {} AAs | {:?} | {} | {} hashes",
//             &pdb_path,
//             compact.num_residues,
//             end - start,
//             before_dedup,
//             hash_collector.hash_collection.len()
//         );
//     }
// } // WORKS with errors in Glycine  2023-03-28 18:23:37

#[test]
fn test_controller() {
    let pdb_paths = loader::load_homeobox_toy();
    let mut controller = Controller::new(pdb_paths);
    controller.fill_numeric_id_vec();
    controller.collect_hash();
    for i in 0..controller.hash_collection_vec.len() {
        println!(
            "{:?} | {:?} | {:?} hashes",
            controller.path_vec.get(i),
            controller.numeric_id_vec.get(i),
            controller
                .hash_collection_vec
                .get(i)
                .unwrap_or(&HashCollection::new())
                .len()
        );
    }
} // WORKS 2023-03-28 20:13:33

#[test]
fn test_index_builder() {
    let pdb_paths = loader::load_homeobox_toy();
    let mut controller = Controller::new(pdb_paths);
    controller.fill_numeric_id_vec();
    controller.collect_hash();

    let index_builder = IndexBuilder::new();
    let index_table =
        index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
    let table_printer = IndexTablePrinter::Debug;
    table_printer.print(&index_table, "data/homeobox_index_table.tsv"); // TODO: Change this to work
    controller.save_raw_feature("data/homeobox_hash_raw.tsv", true);
    println!("{:?}", &controller.path_vec);
}

#[test]
fn test_index_builder_2() {
    let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
    let mut controller = Controller::new(pdb_paths);
    controller.fill_numeric_id_vec();
    controller.collect_hash();

    let index_builder = IndexBuilder::new();
    let index_table =
        index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
    let table_printer = IndexTablePrinter::Debug;
    controller.save_raw_feature("data/serine_hash_raw.tsv", true);
    table_printer.print(&index_table, "data/serine_peptidases_index_table.tsv"); // TODO: Change this to work
    // println!("{:?}", &controller.path_vec);
}

#[test]
fn test_temp() {
    let pdb_paths = loader::load_path("analysis/result/msd_human_s1");
    let mut controller = Controller::new(pdb_paths);
    controller.fill_numeric_id_vec();
    // controller.collect_triad_hash();
    controller.collect_hash();
    controller.save_hash_per_pair("analysis/result/msd_human_s1_hash_per_pair.tsv");
    controller.save_raw_feature("analysis/result/msd_human_s1_hash_raw.tsv", true);
    let index_builder = IndexBuilder::new();
    let index_table =
        index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
    let table_printer = IndexTablePrinter::Debug;
    table_printer.print(&index_table, "analysis/result/msd_human_s1_index_table.tsv");
    // println!("{:?}", &controller.path_vec);
}

#[test]
fn test_querying() {
    let pdb_paths: Vec<String> = loader::load_path("data/serine_peptidases_filtered");
    let mut controller = Controller::new(pdb_paths);
    controller.fill_numeric_id_vec();
    controller.collect_hash();
    let index_builder = IndexBuilder::new();
    let index_table =
        index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
    let query: HashValue = HashValue::from_u64(6597069965577u64);
    let result = query_single(&index_table, &query);
    // Remove duplicates in result
    if let Some(result) = result {
        let mut set = HashSet::new();
        let mut dedup_result = Vec::new();
        for i in result {
            if !set.contains(&i) {
                set.insert(i);
                dedup_result.push(i);
            }
        }
        println!("Queried {:?}", dedup_result);
    }

    // Measure time
    let start = std::time::Instant::now();
    let homeobox_queries = [
        HashValue::from_u64(6597069965577u64),
        HashValue::from_u64(4398046577927u64),
        HashValue::from_u64(8800438323719u64),
    ];

    let homeobox_neighboring_queries: Vec<Vec<HashValue>> = homeobox_queries.iter().map(|x| x.neighbors(true)).collect();
    println!("HASHES (3) {:?}", &homeobox_neighboring_queries);
    
    // let result = query_multiple(&index_table, &homeobox_queries);
    let result = query_multiple_with_neighbors(&index_table, homeobox_neighboring_queries);
    // Remove duplicates in result
    if let Some(result) = result {
        let mut set = HashSet::new();
        let mut dedup_result = Vec::new();
        for i in result {
            if !set.contains(&i) {
                set.insert(i);
                dedup_result.push(i);
            }
        }
        println!("Queried ind (3) {:?}", dedup_result);
        for i in dedup_result {
            println!("Queried prot (3) {:?}", controller.path_vec.get(i));
        }
    }
    let end = std::time::Instant::now();
    println!("Time elapsed for quering neighbors {:?}", end - start);
    let start = std::time::Instant::now();
    let result = query_multiple(&index_table, &homeobox_queries);

    if let Some(result) = result {
        let mut set = HashSet::new();
        let mut dedup_result = Vec::new();
        for i in result {
            if !set.contains(&i) {
                set.insert(i);
                dedup_result.push(i);
            }
        }
        println!("Queried ind (3) {:?}", dedup_result);
        for i in dedup_result {
            println!("Queried prot (3) {:?}", controller.path_vec.get(i));
        }
    }
    let end = std::time::Instant::now();
    println!("Time elapsed for quering {:?}", end - start);
    
    // let printer = IndexTablePrinter::Debug;
    // printer.print(&index_table, "data/yeast_index_table.tsv");
    // println!("{:?}", &controller.path_vec);
}
