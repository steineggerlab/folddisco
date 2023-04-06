use std::io::Write;

use crate::geometry::hash::{HashCollection, HashValue};
use crate::index::builder::IndexBuilder;
use crate::index::*;
use crate::structure::core::CompactStructure;
use crate::structure::io::pdb::Reader as PDBReader;

pub struct Controller {
    pub path_vec: Vec<String>,
    pub numeric_id_vec: Vec<usize>,
    pub hash_collection_vec: Vec<HashCollection>,
    // pub res_pair_vec: Vec<Vec<(u64, u64)>>, // WARNING: TEMPORARY
}

impl Controller {
    pub fn new(path_vec: Vec<String>) -> Controller {
        Controller {
            path_vec: path_vec,
            numeric_id_vec: Vec::new(),
            hash_collection_vec: Vec::new(),
            // res_pair_vec: Vec::new(),// WARNING: TEMPORARY
        }
    }

    pub fn collect_hash(&mut self) {
        for i in 0..self.path_vec.len() {
            let pdb_path = &self.path_vec[i];
            let pdb_reader = PDBReader::from_file(pdb_path).expect("pdb file not found");
            let structure = pdb_reader.read_structure().expect("structure read failed");
            let compact = structure.to_compact();
            let mut hash_collector = GeometryHashCollector::new();

            // let mut res_pair_vec = Vec::new(); // WARNING: TEMPORARY

            for n in 0..compact.num_residues {
                for m in n..compact.num_residues {
                    if n == m {
                        continue;
                    }
                    let dist = compact.get_distance(n, m).expect("cannot get distance");
                    let angle = compact.get_angle(n, m).unwrap_or(0.0);
                    let hash_value = HashValue::perfect_hash(dist, angle);
                    hash_collector.collect_hash(hash_value);

                    // WARNING: TEMPORARY
                    // let res1 = compact.residues[n];
                    // let res2 = compact.residues[m];
                    // res_pair_vec.push((res1, res2));
                }
            }
            hash_collector.remove_redundancy();
            self.hash_collection_vec
                .push(hash_collector.hash_collection);
            // self.res_pair_vec.push(res_pair_vec); // For debug
        }
        // TEMPORARY
        println!("Collected {} pdbs", self.hash_collection_vec.len()); // TEMP
    }

    pub fn fill_numeric_id_vec(&mut self) {
        string_vec_to_numeric_id_vec(&self.path_vec, &mut self.numeric_id_vec);
    }
}

fn string_vec_to_numeric_id_vec(string_vec: &Vec<String>, numeric_id_vec: &mut Vec<usize>) {
    for i in 0..string_vec.len() {
        numeric_id_vec.push(i);
    }
}

pub struct GeometryHashCollector {
    pub hash_collection: HashCollection,
    // TODO: FILL IN OTHER FIELDS
}

impl GeometryHashCollector {
    pub fn new() -> GeometryHashCollector {
        GeometryHashCollector {
            hash_collection: HashCollection::new(),
        }
    }

    pub fn collect_hash(&mut self, hash_value: HashValue) {
        self.hash_collection.push(hash_value);
    }

    pub fn remove_redundancy(&mut self) {
        self.hash_collection.sort();
        self.hash_collection.dedup();
    }
}

// Temporary function for testing
fn _write_hash_with_res_pair(
    hash_collection: &Vec<HashValue>,
    res_pair_vec: &Vec<(u64, u64)>,
    path: &str,
) {
    let mut file = std::fs::File::create(path).expect("Unable to create file");
    file.write_all(b"hash\tdist\tangle\tres1\tres2\n")
        .expect("Unable to write header");

    for j in 0..hash_collection.len() {
        let hash_value = hash_collection[j];
        let res1 = res_pair_vec[j].0;
        let res2 = res_pair_vec[j].1;
        file.write_all(format!("{}\t{}\t{}\n", hash_value, res1, res2).as_bytes())
            .expect("Unable to write data");
    }
}

#[cfg(test)]
mod controller_tests {
    use super::*;
    use crate::index::IndexTablePrinter;
    use crate::structure::io::pdb;
    use crate::test::{load_homeobox_toy, load_path, load_yeast_proteome};

    #[test]
    fn test_geometry_hash_collector() {
        let pdb_paths = load_homeobox_toy();
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

                    let dist = compact.get_distance(i, j).expect("Failed to get distance");
                    // If angle is None, then set it to 0.0. TODO: Glycine should be handled.
                    let angle = compact.get_angle(i, j).unwrap_or(0.0);
                    let hash_value = HashValue::perfect_hash(dist, angle);
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
        } // WORKS with errors in Glycine  2023-03-28 18:23:37
    }

    #[test]
    fn test_controller() {
        let pdb_paths = load_homeobox_toy();
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
        let pdb_paths = load_path("data/toy/fusion");
        let mut controller = Controller::new(pdb_paths);
        controller.fill_numeric_id_vec();
        controller.collect_hash();
        // Save controller hash vectors and residue pairs
        // for i in 0..controller.hash_collection_vec.len() {
        //     let hash_collection = controller.hash_collection_vec.get(i).expect("cannot get hash collection");
        //     let res_pair_vec = controller.res_pair_vec.get(i).expect("cannot get residue pair vector");
        //     _write_hash_with_res_pair(
        //         hash_collection, &res_pair_vec,
        //         &format!("data/homeobox_{}.tsv", i)
        //     );
        // }
        let index_builder = IndexBuilder::new();
        let index_table =
            index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
        let table_printer = IndexTablePrinter::Debug;
        table_printer.print(&index_table, "data/checking_index_table.tsv"); // TODO: Change this to work
        println!("{:?}", &controller.path_vec);
    }

    #[test]
    fn test_querying() {
        let pdb_paths: Vec<String> = load_yeast_proteome();
        let mut controller = Controller::new(pdb_paths);
        controller.fill_numeric_id_vec();
        controller.collect_hash();
        let index_builder = IndexBuilder::new();
        let index_table =
            index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
        let query: HashValue = HashValue::from_u16(3158u16);
        let result = query_single(&index_table, &query);
        println!("{:?}", result);
        let homeobox_queries = [
            HashValue::from_u16(4690u16),
            HashValue::from_u16(4747u16),
            HashValue::from_u16(4232u16),
            HashValue::from_u16(3895u16),
            HashValue::from_u16(2638u16),
            HashValue::from_u16(3961u16),
            HashValue::from_u16(3940u16),
            HashValue::from_u16(3905u16),
            HashValue::from_u16(3937u16),
            HashValue::from_u16(3930u16),
            HashValue::from_u16(3929u16),
            HashValue::from_u16(3896u16),
            HashValue::from_u16(3647u16),
            HashValue::from_u16(3652u16),
            HashValue::from_u16(3675u16),
            HashValue::from_u16(3661u16),
            HashValue::from_u16(3672u16),
            HashValue::from_u16(3670u16),
            HashValue::from_u16(3700u16),
            HashValue::from_u16(3679u16),
            HashValue::from_u16(3632u16),
            HashValue::from_u16(3687u16),
            HashValue::from_u16(3435u16),
            HashValue::from_u16(3396u16),
            HashValue::from_u16(3444u16),
            HashValue::from_u16(3450u16),
            HashValue::from_u16(3445u16),
            HashValue::from_u16(3390u16),
            HashValue::from_u16(3158u16),
            HashValue::from_u16(3147u16),
            HashValue::from_u16(3166u16),
            HashValue::from_u16(3140u16),
            HashValue::from_u16(3129u16),
            HashValue::from_u16(3178u16),
            HashValue::from_u16(3114u16),
            HashValue::from_u16(3134u16),
            HashValue::from_u16(2843u16),
            HashValue::from_u16(2927u16),
            HashValue::from_u16(2848u16),
            HashValue::from_u16(2663u16),
            HashValue::from_u16(2675u16),
            HashValue::from_u16(2667u16),
            HashValue::from_u16(2660u16),
            HashValue::from_u16(2664u16),
            HashValue::from_u16(2409u16),
            HashValue::from_u16(2427u16),
            HashValue::from_u16(1651u16),
            HashValue::from_u16(1650u16),
            HashValue::from_u16(1391u16),
            HashValue::from_u16(1332u16),
            HashValue::from_u16(1392u16),
            HashValue::from_u16(1395u16),
            HashValue::from_u16(1338u16),
            HashValue::from_u16(1396u16),
            HashValue::from_u16(1111u16),
            HashValue::from_u16(1108u16),
            HashValue::from_u16(1100u16),
            HashValue::from_u16(1105u16),
            HashValue::from_u16(1110u16),
            HashValue::from_u16(1104u16),
            HashValue::from_u16(1106u16),
        ];
        let result = query_multiple(&index_table, &homeobox_queries);
        println!("{:?}", result);
        if let Some(result) = result {
            for i in result {
                println!("{:?}", controller.path_vec.get(i));
            }
        }
        let printer = IndexTablePrinter::Debug;
        printer.print(&index_table, "data/yeast_index_table.tsv");
        println!("{:?}", &controller.path_vec);
    }
}
