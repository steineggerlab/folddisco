pub mod query;

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::sync::atomic::AtomicUsize;

// use crate::geometry::simple_hash::{HashCollection, HashValue};
// use crate::geometry::triad_hash::{HashCollection, HashValue};
// use crate::geometry::ppf::{HashCollection, HashValue};
// use crate::geometry::trrosetta::{normalize_angle_degree, HashCollection, HashValue, discretize_value};
use crate::geometry::trrosetta_subfamily::{normalize_angle_degree, HashCollection, HashValue, discretize_value};
// use crate::geometry::aa_pair::{HashCollection, HashValue, discretize_value, map_aa_to_u8};
// use crate::geometry::trrosetta_reduced::*;
use crate::index::builder::IndexBuilder;
use crate::index::*;
use crate::structure::core::CompactStructure;
use crate::structure::io::pdb::Reader as PDBReader;
use rayon::prelude::*;

pub struct Controller {
    pub path_vec: Vec<String>,
    pub numeric_id_vec: Vec<usize>,
    pub hash_collection_vec: Vec<HashCollection>,
    pub hash_new_collection_vec: Vec<crate::geometry::two_float::HashCollection>, // WARNING: TEMPORARY
    pub res_pair_vec: Vec<Vec<(u64, u64, [u8; 3], [u8; 3])>>, // WARNING: TEMPORARY
    pub remove_redundancy: bool,
    pub num_threads_file: usize,
    pub num_threads_hash: usize,
}

impl Controller {
    pub fn new(path_vec: Vec<String>) -> Controller {
        Controller {
            path_vec: path_vec,
            numeric_id_vec: Vec::new(),
            hash_collection_vec: Vec::new(),
            hash_new_collection_vec: Vec::new(),
            res_pair_vec: Vec::new(), // WARNING: TEMPORARY
            remove_redundancy: false,
            num_threads_file: 3,
            num_threads_hash: 2,
        }
    }
    pub fn collect_hash(&mut self) {
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads_file)
            .build()
            .expect("Failed to build thread pool for iterating files");
        // For iterating files, apply multi-threading with num_threads_for_file
        let output = pool.install(|| {
            self.path_vec
                .par_iter()
                .map(|pdb_path| {
                    let hash_collection = par_get_feature_per_structure(pdb_path, self.num_threads_hash);
                    hash_collection
                })
                .collect::<Vec<HashCollection>>()
        }) as Vec<HashCollection>;
        
        println!("Collected {} pdbs", output.len()); // TEMP
        println!("Hash collection: {:?}", output[0][0]); // TEMP
        self.hash_collection_vec = output;
        // Delete pool
        
    }
    
    pub fn get_allocation_size(&self) -> usize {
        let mut allocation_size = 0usize;
        self.path_vec.iter().for_each(|pdb_path| {
            let pdb_reader = PDBReader::from_file(pdb_path).expect("PDB file not found");
            let compact = pdb_reader.read_structure().expect("Failed to read PDB file");
            allocation_size += compact.num_residues * (compact.num_residues - 1);
        });
        allocation_size
    }
    
    pub fn _collect_hash(&mut self) {
        // let vae = load_vae();
        for i in 0..self.path_vec.len() {
            let pdb_path = &self.path_vec[i];
            let pdb_reader = PDBReader::from_file(pdb_path).expect("pdb file not found");
            let structure = pdb_reader.read_structure().expect("structure read failed");
            let compact = structure.to_compact();
            let mut hash_collector = GeometryHashCollector::new();

            let mut res_pairs: Vec<(u64, u64, [u8; 3], [u8; 3])> = Vec::new(); // WARNING: TEMPORARY

            for n in 0..compact.num_residues {
                for m in 0..compact.num_residues {
                    if n == m {
                        continue;
                    }
                    // if n.abs_diff(m) < 3 { // WARNING: TEMPORARY FOR CHECKING IF ACCUMULATED TORSION WORKS
                    //     continue;
                    // }

                    let trr = compact.get_trrosetta_feature(n, m).unwrap_or([0.0; 6]);
                    // let ppf = compact.get_ppf(n, m).expect("cannot get ppf");
                    // let angle = compact.get_accumulated_torsion(n, m).expect("cannot get accumulated torsion angle");
                    if trr[0] < 2.0 || trr[0] > 20.0 {
                        continue;
                    }
                    // let reduced_trr = reduce_with_vae(trr, &vae);
                    // let hash_value = HashValue::perfect_hash(reduced_trr[0], reduced_trr[1]);
                    let hash_value =
                        HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);

                    // WARNING: TEMPORARY
                    let res1 = compact.residue_serial[n];
                    let res2 = compact.residue_serial[m];
                    let res1_str = compact.residue_name[n];
                    let res2_str = compact.residue_name[m];
                    // Convert amino acid three letter code to one letter code
                    // let res1_u8 = map_aa_to_u8(&res1_str);
                    // let res2_u8 = map_aa_to_u8(&res2_str);
                    // let hash_value = HashValue::perfect_hash(res1_u8, res2_u8, n as usize, m as usize, trr[0]);
                    hash_collector.collect_hash(hash_value);
                    
                    res_pairs.push((res1, res2, res1_str, res2_str));
                }
            }
            if self.remove_redundancy {
                hash_collector.remove_redundancy();
            }
            self.hash_collection_vec
                .push(hash_collector.hash_collection);
            self.res_pair_vec.push(res_pairs); // For debug
        }
        // TEMPORARY
        println!("Collected {} pdbs", self.hash_collection_vec.len()); // TEMP
    }

    pub fn save_raw_feature(&mut self, path: &str, discretize: bool) {
        let mut file = std::fs::File::create(path).expect("cannot create file");
        for i in 0..self.path_vec.len() {
            let pdb_path = &self.path_vec[i];
            if i == 100 {
                println!("Processing {}th pdb - {}", i, pdb_path);
            }
            let pdb_reader = PDBReader::from_file(pdb_path).expect("pdb file not found");
            let structure = pdb_reader.read_structure().expect("structure read failed");
            let compact = structure.to_compact();
            for n in 0..compact.num_residues {
                for m in 0..compact.num_residues {
                    if n == m {
                        continue;
                    }
                    let trr = compact.get_trrosetta_feature(n, m).unwrap_or([0.0; 6]);
                    if trr[0] < 2.0 || trr[0] > 20.0 {
                        continue;
                    }
                    let res1 = compact.residue_serial[n];
                    let res2 = compact.residue_serial[m];
                    // [u8;3] to String
                    let res1_aa = String::from_utf8_lossy(&compact.residue_name[n]).to_string();
                    let res2_aa = String::from_utf8_lossy(&compact.residue_name[m]).to_string();
                    // let res1_str = compact.residue_name[n];
                    // let res2_str = compact.residue_name[m];
                    // // Convert amino acid three letter code to one letter code
                    // let res1_u8 = map_aa_to_u8(&res1_str);
                    // let res2_u8 = map_aa_to_u8(&res2_str);

                    // normalize trr[0]
                    // let n_cb_dist = (trr[0] - 2.0) / 18.0;
                    let n_cb_dist = trr[0];
                    let omega = trr[1];
                    let phi1 = trr[2];
                    let phi2 = trr[3];
                    let psi1 = trr[4];
                    let psi2 = trr[5];
                    let logdist = ((n as i32 - m as i32).abs() + 1).ilog2();
                    let hash_value =
                        HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);

                    // let hash_value = HashValue::perfect_hash(res1_u8, res2_u8, n as usize, m as usize, trr[0]);
                    let mut line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        pdb_path, res1, res2, res1_aa, res2_aa, logdist,
                        n_cb_dist, omega, phi1, phi2, psi1, psi2, hash_value
                    );

                    if discretize == true {
                        let n_cb_dist = discretize_value(n_cb_dist, 2.0, 20.0, 12.0);
                        // Torsion angles
                        let omega = discretize_value(omega, -1.0, 1.0, 6.0);
                        let phi1 = discretize_value(phi1, -1.0, 1.0, 6.0);
                        let phi2 = discretize_value(phi2, -1.0, 1.0, 6.0);
                        // Planar angles
                        let psi1 = discretize_value(psi1, 0.0, 180.0, 6.0);
                        let psi2 = discretize_value(psi2, 0.0, 180.0, 6.0);
                        let logdist = ((n as i32 - m as i32).abs() + 1).ilog2();
                        line = format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            pdb_path, res1, res2, res1_aa, res2_aa, logdist,
                            n_cb_dist, omega, phi1, phi2, psi1, psi2, hash_value
                        );
                    }
                    file.write_all(line.as_bytes()).expect("cannot write file");
                }
            }
        }
    }

    pub fn fill_numeric_id_vec(&mut self) {
        string_vec_to_numeric_id_vec(&self.path_vec, &mut self.numeric_id_vec);
    }

    pub fn save_id_vec(&self, path: &str) {
        // Save numeric_id_vec & path_vec as headerless tsv
        let mut file = std::fs::File::create(path).expect("Unable to create file");
        for i in 0..self.numeric_id_vec.len() {
            let numeric_id = self.numeric_id_vec[i];
            let path = &self.path_vec[i];
            file.write_all(format!("{}\t{}\n", numeric_id, path).as_bytes())
                .expect("Unable to write data");
        }
    }
    
    pub fn save_hash_per_pair(&self, path: &str) {
        let mut file = std::fs::File::create(path).expect("Unable to create file");
        file.write_all(b"hash\tval1\tval2\tn1_n2\tres1_ind\tres2_ind\tres1\tres2\tpdb\n")
            .expect("Unable to write header");
        println!("{}", self.hash_new_collection_vec.len());
        for i in 0..self.hash_new_collection_vec.len() {
            let pdb_path = &self.path_vec[i];
            // let new_path = format!("{}_{}.tsv", path, pdb_path.split("/").last().unwrap());
            // println!("Saving to {}", new_path);
            // file.write_all(b"hash\tdist\tn1_nd\tn2_nd\tn1_n2\tres1_ind\tres2_ind\tpdb\n").expect("Unable to write header"); // ppf
            let hash_collection = &self.hash_new_collection_vec[i];
            let res_pair_vec = &self.res_pair_vec[i];
            println!(
                "{}: {} hashes, {} pairs",
                pdb_path,
                hash_collection.len(),
                res_pair_vec.len()
            );
            for j in 0..res_pair_vec.len() {
                let hash_value = hash_collection[j];
                let res1 = res_pair_vec[j].0;
                let res2 = res_pair_vec[j].1;
                let res1_str = res_pair_vec[j].2;
                let res2_str = res_pair_vec[j].3;
                file.write_all(
                    format!(
                        "{}\t{}\t{}\t{:?}\t{:?}\t{}\n",
                        hash_value, res1, res2, res1_str, res2_str, pdb_path
                    )
                    .as_bytes(),
                )
                .expect("Unable to write data");
            }
        }
    }

    pub fn save_filtered_hash_pair(&self, path: &str, res_pair_filter: &HashMap<String, Vec<u64>>) {
        let mut file = std::fs::File::create(path).expect("Unable to create file");
        file.write_all(b"hash\tval1\tval2\tres1_ind\tres2_ind\tres1\tres2\tpdb\n")
            .expect("Unable to write header");
        // file.write_all(b"hash\tdist\tn1_nd\tn2_nd\tn1_n2\tres1_ind\tres2_ind\tpdb\n").expect("Unable to write header");
        for i in 0..self.hash_new_collection_vec.len() {
            let pdb_path = self.path_vec[i].split("/").last().unwrap();
            if !res_pair_filter.contains_key(pdb_path) {
                continue;
            }
            let hash_collection = &self.hash_new_collection_vec[i];
            let res_pair_vec = &self.res_pair_vec[i];
            println!(
                "Saving filtered {}: {} hashes, {} pairs",
                pdb_path,
                hash_collection.len(),
                res_pair_vec.len()
            );
            for j in 0..res_pair_vec.len() {
                let hash_value = hash_collection[j];
                let res1 = res_pair_vec[j].0;
                let res2 = res_pair_vec[j].1;
                let res1_str = res_pair_vec[j].2.to_ascii_uppercase();
                let res2_str = res_pair_vec[j].3.to_ascii_uppercase();
                if res_pair_filter[pdb_path].contains(&res1)
                    && res_pair_filter[pdb_path].contains(&res2)
                {
                    file.write_all(
                        format!(
                            "{}\t{}\t{}\t{:?}\t{:?}\t{}\n",
                            hash_value, res1, res2, res1_str, res2_str, pdb_path
                        )
                        .as_bytes(),
                    )
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

fn get_all_combination(n: usize, include_same: bool) -> Vec<(usize, usize)> {
    let mut res = Vec::new();
    for i in 0..n {
        for j in 0..n {
            if i == j && !include_same {
                continue;
            }
            res.push((i, j));
        }
    }
    res
}

fn par_get_feature_per_structure(pdb_path: &String, num_threads: usize) -> Vec<HashValue> {
    let pdb_reader = PDBReader::from_file(pdb_path).expect("PDB file not found");
    let compact = pdb_reader.read_structure().expect("Failed to read PDB file");
    let compact = compact.to_compact();
    let res_bound = get_all_combination(compact.num_residues, false);
    // Set number of threads
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Failed to build thread pool for iterating residues");
    let output =  pool.install(|| {
        res_bound
            .into_par_iter()
            .map(|(n, m)| {
                let trr = compact.get_trrosetta_feature(n, m).unwrap_or([0.0; 6]);
                if trr[0] >= 2.0 && trr[0] <= 20.0 {
                    return HashValue::from_u64(0u64);
                }
                let hash_value = HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);
                hash_value

            })
            .filter(|x| x != &HashValue::from_u64(0u64))
            .collect::<Vec<HashValue>>()
    });
    output
}

#[cfg(test)]
mod tests {
    use core::num;

    #[test]
    fn test_get_all_combination() {
        let res = super::get_all_combination(3, false);
        assert_eq!(res, vec![(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]);
        let res = super::get_all_combination(3, true);
        assert_eq!(res, vec![(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]);
    }
    
    #[test]
    fn test_par_get_feature_per_structure() {
        let pdb_path = String::from("data/homeobox/1akha-.pdb");
        let num_threads = 1;

        let start = std::time::Instant::now();
        let hash_collection = super::par_get_feature_per_structure(&pdb_path, num_threads);
        let end = std::time::Instant::now();
        println!("Time elapsed {}: {:?}", num_threads, end - start);
        // println!("Hash collection: {:?}", hash_collection);
        let num_threads = 2;
        let start = std::time::Instant::now();
        let hash_collection = super::par_get_feature_per_structure(&pdb_path, num_threads);
        let end = std::time::Instant::now();
        println!("Time elapsed {}: {:?}", num_threads, end - start);
        // println!("Hash collection: {:?}", hash_collection);
        let num_threads = 4;
        let start = std::time::Instant::now();
        let hash_collection = super::par_get_feature_per_structure(&pdb_path, num_threads);
        let end = std::time::Instant::now();
        println!("Time elapsed {}: {:?}", num_threads, end - start);
        // println!("Hash collection: {:?}", hash_collection);
    }
}