use std::io::Write;
use std::collections::HashMap;

use crate::geometry::simple_hash::{HashCollection, HashValue};
use crate::index::builder::IndexBuilder;
use crate::index::*;
use crate::structure::core::CompactStructure;
use crate::structure::io::pdb::Reader as PDBReader;

pub struct Controller {
    pub path_vec: Vec<String>,
    pub numeric_id_vec: Vec<usize>,
    pub hash_collection_vec: Vec<HashCollection>,
    pub res_pair_vec: Vec<Vec<(u64, u64, [u8; 3], [u8; 3])>>, // WARNING: TEMPORARY
    pub remove_redundancy: bool,
}

impl Controller {
    pub fn new(path_vec: Vec<String>) -> Controller {
        Controller {
            path_vec: path_vec,
            numeric_id_vec: Vec::new(),
            hash_collection_vec: Vec::new(),
            res_pair_vec: Vec::new(),// WARNING: TEMPORARY
            remove_redundancy: false,
        }
    }

    pub fn collect_hash(&mut self) {
        for i in 0..self.path_vec.len() {
            let pdb_path = &self.path_vec[i];
            let pdb_reader = PDBReader::from_file(pdb_path).expect("pdb file not found");
            let structure = pdb_reader.read_structure().expect("structure read failed");
            let compact = structure.to_compact();
            let mut hash_collector = GeometryHashCollector::new();

            let mut res_pair_vec = Vec::new(); // WARNING: TEMPORARY

            for n in 0..compact.num_residues {
                for m in n..compact.num_residues {
                    // if n == m {
                    //     continue;
                    // }
                    if n.abs_diff(m) < 3 { // WARNING: TEMPORARY FOR CHECKING IF ACCUMULATED TORSION WORKS
                        continue;
                    }

                    let dist = compact.get_distance(n, m).expect("cannot get distance");
                    // let angle = compact.get_angle(n, m).unwrap_or(0.0);
                    let angle = compact.get_accumulated_torsion(n, m).expect("cannot get accumulated torsion angle");
                    let hash_value = HashValue::perfect_hash(dist, angle);
                    hash_collector.collect_hash(hash_value);

                    // WARNING: TEMPORARY
                    let res1 = compact.residue_serial[n];
                    let res2 = compact.residue_serial[m];
                    let res1_str = compact.residue_name[n];
                    let res2_str = compact.residue_name[m];
                    res_pair_vec.push((res1, res2, res1_str, res2_str));
                }
            }
            if self.remove_redundancy {
                hash_collector.remove_redundancy();
            }
            self.hash_collection_vec
                .push(hash_collector.hash_collection);
            self.res_pair_vec.push(res_pair_vec); // For debug
        }
        // TEMPORARY
        println!("Collected {} pdbs", self.hash_collection_vec.len()); // TEMP
    }

    pub fn fill_numeric_id_vec(&mut self) {
        string_vec_to_numeric_id_vec(&self.path_vec, &mut self.numeric_id_vec);
    }

    pub fn save_hash_per_pair(&self, path: &str) {
        for i in 0..self.hash_collection_vec.len() {
            let pdb_path = &self.path_vec[i];
            let new_path = format!("{}_{}.tsv", path, pdb_path.split("/").last().unwrap());
            println!("Saving to {}", new_path);
            let mut file = std::fs::File::create(new_path).expect("Unable to create file");
            file.write_all(b"hash\tdist\tangle\tres1_ind\tres2_ind\tres1\tres2\tpdb\n").expect("Unable to write header");
            let hash_collection = &self.hash_collection_vec[i];
            let res_pair_vec = &self.res_pair_vec[i];
            println!("{}: {} hashes, {} pairs", pdb_path, hash_collection.len(), res_pair_vec.len());
            for j in 0..res_pair_vec.len(){
                let hash_value = hash_collection[j];
                let res1 = res_pair_vec[j].0;
                let res2 = res_pair_vec[j].1;
                let res1_str = res_pair_vec[j].2;
                let res2_str = res_pair_vec[j].3;
                file.write_all(format!("{}\t{}\t{}\t{:?}\t{:?}\t{}\n", hash_value, res1, res2, res1_str, res2_str, pdb_path).as_bytes())
                    .expect("Unable to write data");
            }
        }
    }

    pub fn save_filtered_hash_pair(&self, path: &str, res_pair_filter: &HashMap<String, Vec<u64>>) {
        let mut file = std::fs::File::create(path).expect("Unable to create file");
        file.write_all(b"hash\tdist\tangle\tres1_ind\tres2_ind\tres1\tres2\tpdb\n").expect("Unable to write header");
        for i in 0..self.hash_collection_vec.len() {
            let pdb_path = self.path_vec[i].split("/").last().unwrap();
            if !res_pair_filter.contains_key(pdb_path) {
                continue;
            }
            let hash_collection = &self.hash_collection_vec[i];
            let res_pair_vec = &self.res_pair_vec[i];
            println!("{}: {} hashes, {} pairs", pdb_path, hash_collection.len(), res_pair_vec.len());
            for j in 0..res_pair_vec.len(){
                let hash_value = hash_collection[j];
                let res1 = res_pair_vec[j].0;
                let res2 = res_pair_vec[j].1;
                let res1_str = String::from_utf8(res_pair_vec[j].2.to_ascii_uppercase()).unwrap();
                let res2_str = String::from_utf8(res_pair_vec[j].3.to_ascii_uppercase()).unwrap();
                if res_pair_filter[pdb_path].contains(&res1) && res_pair_filter[pdb_path].contains(&res2) {
                    file.write_all(format!("{}\t{}\t{}\t{:?}\t{:?}\t{}\n", hash_value, res1, res2, res1_str, res2_str, pdb_path).as_bytes())
                    .expect("Unable to write data");
                }
            }
        }
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
        let pdb_paths = load_homeobox_toy();
        let mut controller = Controller::new(pdb_paths);
        controller.fill_numeric_id_vec();
        controller.collect_hash();
        // controller.save_hash_per_pair("data/homeobox_hash_per_pair.tsv");
        let index_builder = IndexBuilder::new();
        let index_table =
            index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
        let table_printer = IndexTablePrinter::Debug;
        table_printer.print(&index_table, "data/homeobox_index_table.tsv"); // TODO: Change this to work
        println!("{:?}", &controller.path_vec);
    }

    #[test]
    fn test_temp() {
        let pdb_paths = load_path("data/serine_peptidases");
        let mut controller = Controller::new(pdb_paths);
        controller.fill_numeric_id_vec();
        controller.collect_hash();

        let mut serine_filter: HashMap<String, Vec<u64>> = HashMap::new();
        serine_filter.insert(
            "1aq2.pdb".to_string(), vec![250, 232, 269],
        );
        serine_filter.insert(
            "1wab.pdb".to_string(), vec![47, 195, 192],
        );
        serine_filter.insert(
            "1sc9.pdb".to_string(), vec![80, 235, 207],
        );
        serine_filter.insert(
            "2o7r.pdb".to_string(), vec![169, 306, 276],
        );
        serine_filter.insert(
            "1bs9.pdb".to_string(), vec![90, 187, 175],
        );
        serine_filter.insert(
            "1ju3.pdb".to_string(), vec![117, 287, 259],
        );
        serine_filter.insert(
            "1uk7.pdb".to_string(), vec![34, 252, 224],
        );
        serine_filter.insert(
            "1okg.pdb".to_string(), vec![255, 75, 61],
        );
        serine_filter.insert(
            "1qfm.pdb".to_string(), vec![554, 680, 641],
        );


        controller.save_filtered_hash_pair("data/serine_hash_per_pair.tsv", &serine_filter);
        let index_builder = IndexBuilder::new();
        let index_table =
            index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
        let table_printer = IndexTablePrinter::Debug;
        table_printer.print(&index_table, "data/serine_index_table.tsv");
        println!("{:?}", &controller.path_vec);
    }

    #[test]
    fn test_querying() {
        let pdb_paths: Vec<String> = load_yeast_proteome();
        let mut controller = Controller::new(pdb_paths);
        controller.remove_redundancy = true;
        controller.fill_numeric_id_vec();
        controller.collect_hash();

        let index_builder = IndexBuilder::new();
        let index_table =
            index_builder.concat(&controller.numeric_id_vec, &controller.hash_collection_vec);
        let query: HashValue = HashValue::from_u16(4255u16);
        let result = query_single(&index_table, &query);
        println!("{:?}", result);
        let homeobox_queries = [
            HashValue::from_u16(4255u16),
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