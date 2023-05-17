// // use crate::{geometry, structure};
// // use crate::utils::calculator::Calculate;
// // use crate::geometry::hash;
// use motifsearch::structure::core::CompactStructure;
// // use crate::index::IndexTablePrinter;
// // use crate::structure::io::pdb;
// // use crate::controller;
// use motifsearch::test::load_path;

// mod common;

// fn process_pdb(compact : &CompactStructure) {
//     for i in 0..compact.num_residues {
//         for j in i+1..compact.num_residues {
//             let ci = compact.get_ca(i).expect("compact failed to get CA");
//             let cj = compact.get_ca(j).expect("compact failed to get CA");

//             let dist = compact.get_distance(i,j).expect("compact failed to get distance");
//             let angle = compact.get_angle(i,j).expect("compact failed to get angle");
//             let torsion = ci.calc_torsion_angle(&cj, &cj, &cj);

//             let hashvalue = geometry::hash::HashValue::perfect_hash(dist, angle);
//             let reverse = hashvalue.reverse_hash();

//             println!("residue1: {}, residue2: {}, dist: {}=={}, angle: {}=={}", i, j, dist.round(), reverse.0, angle.round(), reverse.1);
//             assert_eq!(dist.round() as u16, reverse.0);
//             assert_eq!(angle.round() as u16, reverse.1);
//         }
//     }
// }

// #[test]
// fn test_serine_peptidase() {
//     let path = "data/serine_peptidases";
//     let pdb_paths = load_path(path);

//     let mut controller = Controller::new(pdb_paths);
//     controller.fill_numeric_id_vec();
//     controller.collect_hash();

// }
