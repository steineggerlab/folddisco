use crate::geometry::hash::{HashCollection, HashValue};
use crate::index::{IndexTable, write_index_table};
use crate::index::builder::IndexBuilder;
use crate::structure::io::pdb::Reader as PDBReader;
use crate::structure::core::CompactStructure;
use crate::test::load_homeobox_toy;

pub struct Controller {
    pub path_vec: Vec<String>,
    pub numeric_id_vec: Vec<usize>,
    pub hash_collection_vec: Vec<HashCollection>,
}

impl Controller {
    pub fn new(path_vec: Vec<String>) -> Controller {
        Controller {
            path_vec: path_vec,
            numeric_id_vec: Vec::new(),
            hash_collection_vec: Vec::new(),
        }
    }

    pub fn collect_hash(&mut self) {
        for i in 0..self.path_vec.len() {
            let pdb_path = &self.path_vec[i];
            println!("reading pdb file: {}", pdb_path);
            let pdb_reader = PDBReader::from_file(pdb_path).expect("pdb file not found");
            let structure = pdb_reader.read_structure().expect("structure read failed");
            let compact = structure.to_compact();
            let mut hash_collector = GeometryHashCollector::new();
            for n in 0..compact.num_residues {
                for m in 0..compact.num_residues {
                    if n == m {
                        continue;
                    }
                    let dist = compact.get_distance(n, m).expect("cannot get distance");
                    let angle = compact.get_angle(n, m).unwrap_or(0.0);
                    let hash_value = HashValue::perfect_hash(dist, angle);
                    hash_collector.collect_hash(hash_value);
                }
            }
            hash_collector.remove_redundancy();
            self.hash_collection_vec.push(hash_collector.hash_collection);
        }
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

#[cfg(test)]
mod controller_tests {
    use super::*;

    #[test]
    fn test_geometry_hash_collector() {
        let pdb_paths = load_homeobox_toy();
        for pdb_path in pdb_paths {
            // Start measure time
            let start = std::time::Instant::now();
            let pdb_reader = PDBReader::from_file(&pdb_path).expect("Failed to read PDB file");
            let structure = pdb_reader.read_structure().expect("Failed to read structure");
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
                &pdb_path, compact.num_residues, end - start, before_dedup,
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
                controller.path_vec.get(i), controller.numeric_id_vec.get(i),
                controller.hash_collection_vec.get(i).unwrap_or(&HashCollection::new()).len()
            );
        }
    } // WORKS 2023-03-28 20:13:33

    #[test]
    fn test_index_builder() {
        let pdb_paths = load_homeobox_toy();
        let mut controller = Controller::new(pdb_paths);
        controller.fill_numeric_id_vec();
        controller.collect_hash();
        let index_builder = IndexBuilder::new();
        let index_table = index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
        write_index_table(&index_table, "data/homeobox_index_table.tsv");
    }

}