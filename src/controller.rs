use std::io::Write;

use crate::geometry::hash::{HashCollection, HashValue};
use crate::index::builder::IndexBuilder;
use crate::index::*;
use crate::structure::core::CompactStructure;
use crate::PDBReader;

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