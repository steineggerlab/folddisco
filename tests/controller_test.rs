use std::collections::HashMap;

use motifsearch::controller::{self, Controller, GeometryHashCollector};
use motifsearch::geometry::trrosetta::{HashCollection, HashValue};
use motifsearch::index::builder::IndexBuilder;
use motifsearch::index::IndexTablePrinter;
use motifsearch::PDBReader;

mod common;
use common::loader;

#[test]
fn test_geometry_hash_collector() {
    let pdb_paths = loader::load_homeobox_toy();
    for pdb_path in pdb_paths {
        // Start measure time
        let start = std::time::Instant::now();
        let pdb_reader = PDBReader::from_file(&pdb_path).expect("Failed to read PDB file");
        let structure = pdb_reader
            .read_structure()
            .expect("Failed to read structure");
        let compact = structure.to_compact();

        let mut hash_collector = GeometryHashCollector::new();
        for i in 0..compact.num_residues {
            for j in 0..compact.num_residues {
                if i == j {
                    continue;
                }
                let trr = compact
                    .get_trrosetta_feature(i, j)
                    .expect("Failed to get trrosetta feature");
                let hash_value =
                    HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);
                hash_collector.collect_hash(hash_value);
            }
        }
        //
        let before_dedup = hash_collector.hash_collection.len();
        hash_collector.remove_redundancy();
        let end = std::time::Instant::now();
        println!(
            "{:?} | {} AAs | {:?} | {} | {} hashes",
            &pdb_path,
            compact.num_residues,
            end - start,
            before_dedup,
            hash_collector.hash_collection.len()
        );
    }
} // WORKS with errors in Glycine  2023-03-28 18:23:37

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
    controller.save_raw_feature("data/homeobox_hash_raw. tsv", true);
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
    println!("{:?}", &controller.path_vec);
}

    #[test]
    fn test_temp() {
        let pdb_paths = loader::load_path("data/serine_peptidases_filtered");
        let mut controller = Controller::new(pdb_paths);
        controller.fill_numeric_id_vec();
        // controller.collect_triad_hash();
        controller.collect_hash();

        let mut serine_filter: HashMap<String, Vec<u64>> = HashMap::new();
        serine_filter.insert("1aq2.pdb".to_string(), vec![250, 232, 269]);
        serine_filter.insert("1wab.pdb".to_string(), vec![47, 195, 192]);
        serine_filter.insert("1sc9.pdb".to_string(), vec![80, 235, 207]);
        serine_filter.insert("2o7r.pdb".to_string(), vec![169, 306, 276]);
        serine_filter.insert("1bs9.pdb".to_string(), vec![90, 187, 175]);
        serine_filter.insert("1ju3.pdb".to_string(), vec![117, 287, 259]);
        serine_filter.insert("1uk7.pdb".to_string(), vec![34, 252, 224]);
        serine_filter.insert("1okg.pdb".to_string(), vec![255, 75, 61]);
        serine_filter.insert("1qfm.pdb".to_string(), vec![554, 680, 641]);
        controller.save_hash_per_pair("data/archive/serine_hash_per_pair.tsv");
        controller.save_filtered_hash_pair("data/archive/serine_hash_per_pair_filtered.tsv", &serine_filter);
        let index_builder = IndexBuilder::new();
        let index_table =
            index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
        let table_printer = IndexTablePrinter::Debug;
        table_printer.print(&index_table, "data/archive/serine_index_table.tsv");
        println!("{:?}", &controller.path_vec);
    }



// #[test]
// fn test_querying() {
//     let pdb_paths: Vec<String> = load_yeast_proteome();
//     let mut controller = Controller::new(pdb_paths);
//     controller.fill_numeric_id_vec();
//     controller.collect_hash();
//     let index_builder = IndexBuilder::new();
//     let index_table =
//         index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
//     let query: HashValue = HashValue::from_u16(3158u16);
//     let result = query_single(&index_table, &query);
//     println!("{:?}", result);
//     let homeobox_queries = [
//         HashValue::from_u16(3698u16),
//         HashValue::from_u16(2684u16),
//         HashValue::from_u16(2409u16),
//         HashValue::from_u16(3390u16),
//         HashValue::from_u16(3440u16),
//     ];
//     let result = query_multiple(&index_table, &homeobox_queries);
//     println!("{:?}", result);
//     if let Some(result) = result {
//         for i in result {
//             println!("{:?}", controller.path_vec.get(i));
//         }
//     }
//     let printer = IndexTablePrinter::Debug;
//     printer.print(&index_table, "data/yeast_index_table.tsv");
//     println!("{:?}", &controller.path_vec);
// }
