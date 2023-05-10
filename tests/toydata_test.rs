// use crate::{geometry, structure};
// use crate::utils::calculator::Calculate;
// use crate::geometry::hash;
// use motifsearch::structure::core::CompactStructure;
// use crate::index::IndexTablePrinter;
// use crate::structure::io::pdb;
// use crate::controller;

mod common;

#[test]
fn test_serine_peptidase() {
    let path = "data/serine_peptidases";
    let pdb_paths = common::load_path(path);
    let pdb_features = common::process_pdbs(&pdb_paths);
    //TODO:
    // - change struct Torsion HashMap -> Vector
    // - implement code to save Feature Data
}