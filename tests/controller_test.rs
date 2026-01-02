
// Import unique



use folddisco::prelude::*;
use rayon::prelude::*;

mod common;
use common::loader;

#[test]
fn test_folddisco_integration() {
    // Indexing
    // Test if the default hashing schemes are working
    let pdb_paths = loader::load_path("data/serine_peptidases");
    
    let mut fold_disco = Folddisco::create_with_hash_type(pdb_paths, HashType::PDBMotifSinCos);
    measure_time!(fold_disco.collect_hash_vec());
    fold_disco.fill_numeric_id_vec();
    let mut hashes = fold_disco.hash_id_vec.clone();
    // Get the memory usage of hashes
    let _size = std::mem::size_of_val(&hashes);
    let _total_size = std::mem::size_of_val(&hashes[0]) * hashes.len();
    // sort hashes
    measure_time!(hashes.par_sort_by(|a, b| a.0.cmp(&b.0)));


    // Querying
}

