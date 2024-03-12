use std::collections::HashMap;
// Import unique
use std::collections::HashSet;


use motifsearch::prelude::*;
use rayon::prelude::*;

mod common;
use common::loader;

#[test]
fn test_folddisco_default() {
    // Test if the default hashing schemes are working
    let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
    let mut fold_disco = FoldDisco::new(pdb_paths);
    measure_time!(fold_disco.collect_hash());
    // fold_disco.fill_numeric_id_vec();
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
    // fold_disco.fill_numeric_id_vec();
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
    measure_time!(fold_disco.collect_hash_pairs());
    fold_disco.fill_numeric_id_vec();
    let mut hashes = fold_disco.hash_id_pairs.clone();
    // Get the memory usage of hashes
    let size = std::mem::size_of_val(&hashes);
    let total_size = std::mem::size_of_val(&hashes[0]) * hashes.len();
    // sort hashes
    measure_time!(hashes.par_sort_by(|a, b| a.0.cmp(&b.0)));
    let (offset, values) = measure_time!(convert_sorted_pairs_to_offset_and_values_vec(hashes));
}


#[test]
fn test_folddisco_pdbmotifsincos_dashmap() {
    // Test if the default hashing schemes are working
    let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
    let mut fold_disco = FoldDisco::new_with_hash_type(pdb_paths, HashType::PDBMotifSinCos);
    measure_time!(fold_disco.collect_hash());
    // fold_disco.fill_numeric_id_vec();
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
fn test_folddisco_trrosetta() {
    // Test if the default hashing schemes are working
    let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
    let mut fold_disco = FoldDisco::new_with_hash_type(pdb_paths, HashType::TrRosetta);
    fold_disco.collect_hash();
    // fold_disco.fill_numeric_id_vec();
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

// use std::fs::File;
// TODO: IMPLEMENT WRITING RESIDUE COUNTS TO LOOKUP
// #[test]
// fn test_temp() {
//     let pdb_paths = loader::load_path("analysis/h_sapiens_pdb");
//     let mut output = File::create("analysis/h_sapiens_pdb/num_residues.txt").unwrap();
//     for pdb_path in pdb_paths {
//         let pdb = PDBReader::from_file(pdb_path.clone()).unwrap();
//         let structure = pdb.read_structure().unwrap();
//         // Print path and number of residues to a file
//         use std::io::Write;
//         writeln!(output, "{}\t{}", pdb_path, structure.num_residues).unwrap();
//     }
// }