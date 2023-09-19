mod common;

#[test]
fn test_serine_peptidase() {
    let path = "data/serine_peptidases_filtered";
    let pdb_paths = common::loader::load_path(path);
    let pdb_features = common::processor::process_pdbs(&pdb_paths);
    //TODO:
    // - change struct Torsion HashMap -> Vector
    // - implement code to save Feature Data
}
