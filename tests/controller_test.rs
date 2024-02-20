use std::collections::HashMap;
// Import unique
use std::collections::HashSet;

use motifsearch::controller::FoldDisco;
use motifsearch::geometry::core::HashType;
// use motifsearch::geometry::trrosetta::{HashCollection, HashValue};
// use motifsearch::geometry::trrosetta_subfamily::{HashCollection, HashValue};
// use motifsearch::geometry::aa_pair::{HashCollection, HashValue};
use motifsearch::index::builder::IndexBuilder;
use motifsearch::index::query_multiple_with_neighbors;
use motifsearch::index::{IndexTablePrinter, query_single, query_multiple};
use motifsearch::PDBReader;
use motifsearch::index::index_table;
use motifsearch::measure_time;

use rayon::prelude::*;

mod common;
use common::loader;

#[test]
fn test_folddisco_default() {
    // Test if the default hashing schemes are working
    let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
    let mut fold_disco = FoldDisco::new(pdb_paths);
    measure_time!(fold_disco.collect_hash());
    fold_disco.fill_numeric_id_vec();
    for i in 0..fold_disco.hash_collection.len() {
        println!(
            "{:?} | {:?} | {:?} hashes | {:?}",
            fold_disco.numeric_id_vec.get(i).unwrap(),
            fold_disco.path_vec.get(i).unwrap(),
            fold_disco.hash_collection.get(i).unwrap_or(&Vec::new()).len(),
            fold_disco.hash_collection.get(i).unwrap().get(0).unwrap(), 
        );
    }
    measure_time!(fold_disco.set_index_table());
    measure_time!(fold_disco.fill_index_table());
}

#[test]
fn test_folddisco_pdbmotif() {
    // Test if the default hashing schemes are working
    let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
    let mut fold_disco = FoldDisco::new_with_hash_type(pdb_paths, HashType::PDBMotif);
    fold_disco.collect_hash();
    fold_disco.fill_numeric_id_vec();
    for i in 0..fold_disco.hash_collection.len() {
        println!(
            "{:?} | {:?} | {:?} hashes | {:?}",
            fold_disco.numeric_id_vec.get(i).unwrap(),
            fold_disco.path_vec.get(i).unwrap(),
            fold_disco.hash_collection.get(i).unwrap_or(&Vec::new()).len(),
            fold_disco.hash_collection.get(i).unwrap().get(0).unwrap(), 
        );
    }
    fold_disco.set_index_table();
    fold_disco.fill_index_table();
}


#[test]
fn test_folddisco_pdbmotifsincos() {
    // Test if the default hashing schemes are working
    let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
    let mut fold_disco = FoldDisco::new_with_hash_type(pdb_paths, HashType::PDBMotifSinCos);
    fold_disco.collect_hash();
    fold_disco.fill_numeric_id_vec();
    for i in 0..fold_disco.hash_collection.len() {
        println!(
            "{:?} | {:?} | {:?} hashes | {:?}",
            fold_disco.numeric_id_vec.get(i).unwrap(),
            fold_disco.path_vec.get(i).unwrap(),
            fold_disco.hash_collection.get(i).unwrap_or(&Vec::new()).len(),
            fold_disco.hash_collection.get(i).unwrap().get(0).unwrap(), 
        );
    }
    fold_disco.set_index_table();
    fold_disco.fill_index_table();
}


#[test]
fn test_folddisco_trrosetta() {
    // Test if the default hashing schemes are working
    let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
    let mut fold_disco = FoldDisco::new_with_hash_type(pdb_paths, HashType::TrRosetta);
    fold_disco.collect_hash();
    fold_disco.fill_numeric_id_vec();
    for i in 0..fold_disco.hash_collection.len() {
        println!(
            "{:?} | {:?} | {:?} hashes | {:?}",
            fold_disco.numeric_id_vec.get(i).unwrap(),
            fold_disco.path_vec.get(i).unwrap(),
            fold_disco.hash_collection.get(i).unwrap_or(&Vec::new()).len(),
            fold_disco.hash_collection.get(i).unwrap().get(0).unwrap(), 
        );
    }
    fold_disco.set_index_table();
    fold_disco.fill_index_table();
}




// #[test]
// fn test_controller() {
//     let pdb_paths = loader::load_homeobox_toy();
//     let mut controller = Controller::new(pdb_paths);
//     controller.fill_numeric_id_vec();
//     controller.collect_hash();
//     for i in 0..controller.hash_collection_vec.len() {
//         println!(
//             "{:?} | {:?} | {:?} hashes",
//             controller.path_vec.get(i),
//             controller.numeric_id_vec.get(i),
//             controller
//                 .hash_collection_vec
//                 .get(i)
//                 .unwrap_or(&HashCollection::new())
//                 .len()
//         );
//     }
//     // TODO: make it to assertion
// } // WORKS 2023-03-28 20:13:33

// #[test]
// fn test_index_builder() {
//     let pdb_paths = loader::load_homeobox_toy();
//     let mut controller = Controller::new(pdb_paths);
//     controller.fill_numeric_id_vec();
//     controller.collect_hash();

//     let index_builder = IndexBuilder::new();
//     let index_table =
//         index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
//     let table_printer = IndexTablePrinter::Debug;
//     table_printer.print(&index_table, "data/homeobox_index_table.tsv"); // TODO: Change this to work
//     controller.save_raw_feature("data/homeobox_hash_raw.tsv", true);
//     println!("{:?}", &controller.path_vec);
// }

// #[test]
// fn test_index_builder_2() {
//     let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
//     let mut controller = Controller::new(pdb_paths);
//     controller.fill_numeric_id_vec();
//     controller.collect_hash();

//     let index_builder = IndexBuilder::new();
//     let index_table =
//         index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
//     let table_printer = IndexTablePrinter::Debug;
//     controller.save_raw_feature("data/serine_hash_raw.tsv", true);
//     table_printer.print(&index_table, "data/serine_peptidases_index_table.tsv"); // TODO: Change this to work
//     // println!("{:?}", &controller.path_vec);
// }

// #[test]
// fn test_querying() {
//     // let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
//     let pdb_paths: Vec<String> = loader::load_path("analysis/test");
//     let start = std::time::Instant::now();
//     let mut controller = Controller::new(pdb_paths);
//     controller.fill_numeric_id_vec();
//     controller.collect_hash();
//     let index_builder = IndexBuilder::new();
//     let index_table =
//         index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
//     let lap1 = std::time::Instant::now();
//     println!("[OLD] Time elapsed for building index table {:?}", lap1 - start);

//     let query: HashValue = HashValue::from_u64(6597069965577u64);
//     let result = query_single(&index_table, &query);
//     // Remove duplicates in result
//     if let Some(result) = result {
//         let mut set = HashSet::new();
//         let mut dedup_result = Vec::new();
//         for i in result {
//             if !set.contains(&i) {
//                 set.insert(i);
//                 dedup_result.push(i);
//             }
//         }
//         println!("Queried {:?}", dedup_result);
//     }
//     let lap2 = std::time::Instant::now();
//     println!("[OLD] Time elapsed for single query {:?}", lap2 - lap1);
    
//     // Measure time
//     let start = std::time::Instant::now();
//     let homeobox_queries = [
//         HashValue::from_u64(6597069965577u64),
//         HashValue::from_u64(4398046577927u64),
//         HashValue::from_u64(8800438323719u64),
//     ];

//     let homeobox_neighboring_queries: Vec<Vec<HashValue>> = homeobox_queries.iter().map(|x| x.neighbors(true)).collect();
//     // println!("HASHES (3) {:?}", &homeobox_neighboring_queries);
    
//     // let result = query_multiple(&index_table, &homeobox_queries);
//     let result = query_multiple_with_neighbors(&index_table, homeobox_neighboring_queries);
//     // Remove duplicates in result
//     if let Some(result) = result {
//         let mut set = HashSet::new();
//         let mut dedup_result = Vec::new();
//         for i in result {
//             if !set.contains(&i) {
//                 set.insert(i);
//                 dedup_result.push(i);
//             }
//         }
//         // println!("Queried ind (3) {:?}", dedup_result);
//         for i in dedup_result {
//             println!("Queried prot (3) {:?}", controller.path_vec.get(i));
//         }
//     }
//     let lap3 = std::time::Instant::now();
//     println!("[OLD] Time elapsed for quering neighbors {:?}", lap3 - lap2);

//     let result = query_multiple(&index_table, &homeobox_queries);

//     if let Some(result) = result {
//         let mut set = HashSet::new();
//         let mut dedup_result = Vec::new();
//         for i in result {
//             if !set.contains(&i) {
//                 set.insert(i);
//                 dedup_result.push(i);
//             }
//         }
//         // println!("Queried ind (3) {:?}", dedup_result);
//         for i in dedup_result {
//             println!("[OLD] Queried prot (3) {:?}", controller.path_vec.get(i));
//         }
//     }
//     let end = std::time::Instant::now();
//     println!("[OLD] Time elapsed for quering {:?}", end - lap3);
// }


// #[test]
// fn test_querying_using_new_index_table() {
//     rayon::ThreadPoolBuilder::new().num_threads(6).build_global().unwrap();
//     // let pdb_paths: Vec<String> = loader::load_path("data/serine_peptidases_filtered");
//     // let pdb_paths: Vec<String> = loader::load_path("analysis/test");
//     let pdb_paths: Vec<String> = loader::load_path("analysis/raw_ecoli");

//     let start = std::time::Instant::now();
//     let mut controller = Controller::new(pdb_paths);
//     controller.fill_numeric_id_vec();
//     let mut index_table = index_table::IndexTable::new();
    
//     for i in 0.. controller.path_vec.len() {
//         let pdb_path = &controller.path_vec[i];
//         let pdb_reader = PDBReader::from_file(pdb_path).expect("PDB file not found");
//         let compact = pdb_reader.read_structure().expect("Failed to read PDB file");
//         let compact = compact.to_compact();
        
//         let mut hash_vec: Vec<u64> = Vec::with_capacity(compact.num_residues.pow(2));
//         let mut res_pair: Vec<(u16, u16)> = Vec::with_capacity(compact.num_residues.pow(2));
        
        
//         // for n in 0..compact.num_residues {
//         //     for m in 0..compact.num_residues {
//         // Generate a combination of n and m
//         let res_bound = get_all_combination(compact.num_residues, false);
//         let r: Vec<(u64, (u16, u16))> = (&res_bound).into_par_iter().map(
//             |(n, m)| {
//                 if n == m {
//                     return (0u64, (0u16, 0u16))
//                 }
//                 let trr = compact.get_trrosetta_feature(*n, *m).unwrap_or([0.0; 6]);
//                 if trr[0] < 2.0 || trr[0] > 20.0 {
//                     return (0u64, (0u16, 0u16))
//                 }
//                 let hash_value = HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);
//                 // hash_vec.push(hash_value.as_u64());
                
//                 let res1 = compact.residue_serial[*n] as u16;
//                 let res2 = compact.residue_serial[*m] as u16;
//                 // res_pair.push((res1, res2));
//                 (hash_value.as_u64(), (res1, res2))
//                 // println!("n: {:?} m: {:?} hash: {:?} res1: {:?} res2: {:?}", n, m, hash_value.as_u64(), res1, res2);
//             }
//         ).collect();
//         // Filter out the 0s
//         let r = r.into_iter().filter(|(x, _)| *x != 0).collect::<Vec<_>>();
//         hash_vec = r.iter().map(|(x, _)| *x).collect::<Vec<_>>();
//         res_pair = r.iter().map(|(_, x)| *x).collect::<Vec<_>>();
//         index_table.concat_structure(i, hash_vec);
//     }
//     let lap1 = std::time::Instant::now();
//     println!("[NEW] Time elapsed for building index table {:?}", lap1 - start);
    
//     let result = index_table.query_single(&6597069965577u64);
//     // Remove duplicates in result
//     if let Some(result) = result {
//         let mut set = HashSet::new();
//         let mut dedup_result = Vec::new();
//         for i in result {
//             if !set.contains(&i) {
//                 set.insert(i);
//                 dedup_result.push(i);
//             }
//         }
//         println!("[NEW] Single queried {:?}", dedup_result);
//         let queried_ind = dedup_result.iter().map(|x| controller.path_vec.get(x.get_id())).collect::<Vec<_>>();
//         println!("[NEW] Single queried {:?}", queried_ind);
//     }
//     let lap2 = std::time::Instant::now();
//     println!("[NEW] Time elapsed for single query {:?}", lap2 - lap1);

//     // Measure time
//     let homeobox_queries = vec![
//         HashValue::from_u64(6597069965577u64),
//         HashValue::from_u64(4398046577927u64),
//         HashValue::from_u64(8800438323719u64),
//     ];
//     let homeobox_queries_as_u64_ref = homeobox_queries.iter().map(|x| x.as_u64()).collect::<Vec<u64>>();

//     let homeobox_neighboring_queries: Vec<Vec<u64>> = homeobox_queries.iter().map(|x| x.neighbors(true).iter().map(|x| x.as_u64()).collect()).collect();
//     // println!("HASHES (3) {:?}", &homeobox_neighboring_queries);
    
//     // let result = query_multiple(&index_table, &homeobox_queries);
//     let result = index_table.query_multiple_with_neighbors(&homeobox_neighboring_queries);
//     // Remove duplicates in result
//     if let Some(result) = result {
//         let mut set = HashSet::new();
//         let mut dedup_result = Vec::new();
//         for i in result {
//             if !set.contains(&i) {
//                 set.insert(i);
//                 dedup_result.push(i);
//             }
//         }
//         // println!("Queried ind (3) {:?}", dedup_result);
//         for i in dedup_result {
//             println!("[NEW] Queried prot (3) {:?}", controller.path_vec.get(i.get_id()));
//         }
//     }
//     let lap3 = std::time::Instant::now();
//     println!("[NEW] Time elapsed for quering neighbors {:?}", lap3 - lap2);
//     let result = index_table.query_multiple(&homeobox_queries_as_u64_ref);

//     if let Some(result) = result {
//         let mut set = HashSet::new();
//         let mut dedup_result = Vec::new();
//         for i in result {
//             if !set.contains(&i) {
//                 set.insert(i);
//                 dedup_result.push(i);
//             }
//         }
//         // println!("Queried ind (3) {:?}", dedup_result);
//         for i in dedup_result {
//             println!("[NEW] Queried prot (3) {:?}", controller.path_vec.get(i.get_id()));
//         }
//     }
//     let end = std::time::Instant::now();
//     println!("[NEW] Time elapsed for quering {:?}", end - lap3);

//     // Save
//     let result = index_table.save_to_bin("analysis/test_index_table.bin");
//     if let Err(e) = result {
//         println!("Error saving index table to binary: {:?}", e);
//     }
//     let result = index_table.save_to_json("analysis/test_index_table.json");
//     if let Err(e) = result {
//         println!("Error saving index table to json: {:?}", e);
//     }

//     let new_index_table = index_table::IndexTable::load_from_bin("analysis/test_index_table.bin").expect("Failed to load index table from binary");
//     let result = new_index_table.query_multiple(&homeobox_queries_as_u64_ref);
//     println!("Result from loading binary: {:?}", result);
    
// }
