// File: mod.rs
// Created: 2024-01-18 15:47:16
// Description:
//    a new controller implementation that supports multiple hash types
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

pub mod feature;
pub mod filter;
pub mod graph;
pub mod io;
pub mod query;
pub mod result;
pub mod retrieve;
pub mod count_query;
pub mod mode;
pub mod sort;
pub mod summary;

use std::cell::UnsafeCell;
use std::io::Write;
use std::sync::Arc;
use feature::get_geometric_hash_as_u32_from_structure;
use io::read_structure_from_path;
// External imports
use rayon::prelude::*;

use crate::index::indextable::FolddiscoIndex;
// Internal imports
use crate::PDBReader;
use crate::geometry::core::HashType;
use crate::utils::log::{ print_log_msg, log_msg, FAIL, WARN, INFO };

#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;

// Constants
const DEFAULT_NUM_THREADS: usize = 4;
// const DEFAULT_HASH_TYPE: HashType = HashType::PDBTrRosetta;
const DEFAULT_MAX_RESIDUE: usize = 65535;
const DEFAULT_DIST_CUTOFF: f32 = 20.0;

// Module specific types
pub type ResidueMatch = Option<(u8, u64)>;

unsafe impl Send for Folddisco {}
unsafe impl Sync for Folddisco {}

pub struct Folddisco {
    pub path_vec: Vec<String>,
    pub numeric_id_vec: Vec<usize>,
    pub numeric_db_key_vec: Vec<usize>,
    pub nres_vec: Vec<usize>,
    pub plddt_vec: Vec<f32>,
    pub hash_id_vec: Vec<(u32, usize)>,
    pub hash_type: HashType,
    pub num_threads: usize,
    pub num_bin_dist: usize,
    pub num_bin_angle: usize,
    pub multiple_bins: Option<Vec<(usize, usize)>>,
    pub output_path: String,
    pub max_residue: usize,
    pub dist_cutoff: f32,
    pub fold_disco_index: FolddiscoIndex,
    pub foldcomp_db_path: String,
    #[cfg(not(feature = "foldcomp"))]
    pub foldcomp_db_reader: bool,
    #[cfg(feature = "foldcomp")]
    pub foldcomp_db_reader: FoldcompDbReader,
    pub is_foldcomp_enabled: bool,
}

impl Folddisco {
    // Constructors
    pub fn create_with_hash_type(path_vec: Vec<String>, hash_type: HashType) -> Folddisco {
        let length = path_vec.len();
        let total_hashes = 2usize.pow(hash_type.encoding_bits() as u32);
        Folddisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::with_capacity(length),
            numeric_db_key_vec: Vec::new(),
            nres_vec: Vec::with_capacity(length),
            plddt_vec: Vec::with_capacity(length),
            hash_id_vec: Vec::new(),
            hash_type: hash_type,
            num_threads: DEFAULT_NUM_THREADS,
            num_bin_dist: 0,
            num_bin_angle: 0,
            multiple_bins: None,
            output_path: String::new(),
            max_residue: DEFAULT_MAX_RESIDUE,
            dist_cutoff: DEFAULT_DIST_CUTOFF,
            fold_disco_index: FolddiscoIndex::new(total_hashes, String::new(), false),
            foldcomp_db_path: String::new(),
            #[cfg(not(feature = "foldcomp"))]
            foldcomp_db_reader: false,
            #[cfg(feature = "foldcomp")]
            foldcomp_db_reader: FoldcompDbReader::empty(),
            is_foldcomp_enabled: false,
        }
    }

    pub fn new(
        path_vec: Vec<String>, hash_type: HashType, num_threads: usize,
        num_bin_dist: usize, num_bin_angle: usize, output_path: String,
        dist_cutoff: f32, multiple_bins: Option<Vec<(usize, usize)>>,
        mmap_on_disk: bool,
    ) -> Folddisco {
        let length = path_vec.len();
        let total_hashes = 2usize.pow(hash_type.encoding_bits() as u32);
        Folddisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::with_capacity(length),
            numeric_db_key_vec: Vec::new(),
            nres_vec: Vec::with_capacity(length),
            plddt_vec: Vec::with_capacity(length),
            hash_id_vec: Vec::new(),
            hash_type: hash_type,
            num_threads: num_threads,
            num_bin_dist: num_bin_dist,
            num_bin_angle: num_bin_angle,
            multiple_bins: multiple_bins,
            output_path: output_path.clone(),
            max_residue: DEFAULT_MAX_RESIDUE,
            dist_cutoff: dist_cutoff,
            fold_disco_index: FolddiscoIndex::new(total_hashes, output_path.clone(), mmap_on_disk),
            foldcomp_db_path: String::new(),
            #[cfg(not(feature = "foldcomp"))]
            foldcomp_db_reader: false,
            #[cfg(feature = "foldcomp")]
            foldcomp_db_reader: FoldcompDbReader::empty(),
            is_foldcomp_enabled: false,
        }
    }

    #[cfg(feature = "foldcomp")]
    pub fn new_with_foldcomp_db(
        path_vec: Vec<String>, hash_type: HashType, num_threads: usize,
        num_bin_dist: usize, num_bin_angle: usize, output_path: String,
        dist_cutoff: f32, foldcomp_db_path: &'static str,
        multiple_bins: Option<Vec<(usize, usize)>>,
        mmap_on_disk: bool,
    ) -> Folddisco {
        let length = path_vec.len();
        let total_hashes = 2usize.pow(hash_type.encoding_bits() as u32);
        let foldcomp_db_reader = FoldcompDbReader::new(&foldcomp_db_path);
        
        Folddisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::with_capacity(length),
            numeric_db_key_vec: foldcomp_db_reader.get_db_key_vector(),
            nres_vec: Vec::with_capacity(length),
            plddt_vec: Vec::with_capacity(length),
            hash_id_vec: Vec::new(),
            hash_type: hash_type,
            num_threads: num_threads,
            num_bin_dist: num_bin_dist,
            num_bin_angle: num_bin_angle,
            multiple_bins: multiple_bins,
            output_path: output_path.clone(),
            max_residue: DEFAULT_MAX_RESIDUE,
            dist_cutoff: dist_cutoff,
            fold_disco_index: FolddiscoIndex::new(total_hashes, output_path.clone(), mmap_on_disk),
            foldcomp_db_path: foldcomp_db_path.to_string(),
            foldcomp_db_reader: foldcomp_db_reader,
            is_foldcomp_enabled: true,
        }
    }

    // Setters
    pub fn set_path_vec(&mut self, path_vec: Vec<String>) {
        self.path_vec = path_vec;
    }
    pub fn set_hash_type(&mut self, hash_type: HashType) {
        self.hash_type = hash_type;
    }
    pub fn set_num_threads(&mut self, num_threads: usize) {
        self.num_threads = num_threads;
    }
    pub fn set_output_path(&mut self, output_path: String) {
        self.output_path = output_path;
    }
    pub fn set_num_bin_dist(&mut self, num_bin_dist: usize) {
        self.num_bin_dist = num_bin_dist;
    }
    pub fn set_num_bin_angle(&mut self, num_bin_angle: usize) {
        self.num_bin_angle = num_bin_angle;
    }
    pub fn set_max_residue(&mut self, max_residue: usize) {
        self.max_residue = max_residue;
    }
    
    // Main methods
    pub fn fill_numeric_id_vec(&mut self) {
        string_vec_to_numeric_id_vec(&self.path_vec, &mut self.numeric_id_vec);
    }

    pub fn collect_hash_vec(&mut self) { // THISONE
        // Mutex free version
        let shared_data = SharedData::new(self.path_vec.len());
        
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating files");
        // For iterating files, apply multi-threading with num_threads_for_file

        let collected: Vec<(u32, usize)> = pool.install(|| {
            // Preserve locality for multi-threading
            self.path_vec
                .par_iter()
                .enumerate()
                .map(|(pdb_pos, pdb_path)| {
                    #[cfg(not(feature = "foldcomp"))]
                    let compact = read_structure_from_path(pdb_path).expect(
                        log_msg(FAIL, "Failed to read structure").as_str()
                    );

                    #[cfg(feature = "foldcomp")]
                    let compact = if self.is_foldcomp_enabled {
                        self.foldcomp_db_reader.read_single_structure_by_id(
                            self.numeric_db_key_vec[pdb_pos]
                        ).expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        )
                    } else {
                        read_structure_from_path(pdb_path).expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        )
                    };  

                    if compact.num_residues > self.max_residue {
                        print_log_msg(WARN, &format!("{} has too many residues. Skipping", pdb_path));
                        // skip this file
                        // Drop intermediate variables
                        drop(compact);
                        return Vec::new();
                    }
                    let compact = compact.to_compact();
                    // Directly write num_residues and avg_plddt to the vectors
                    let nres = compact.num_residues;
                    let plddt = compact.get_avg_plddt();
                    // Mutex free version
                    unsafe {
                        let nres_vec = shared_data.get_nres_vec();
                        let plddt_vec = shared_data.get_plddt_vec();
                        let nres_vec = &mut *nres_vec.get();
                        let plddt_vec = &mut *plddt_vec.get();
                        nres_vec[pdb_pos] = nres;
                        plddt_vec[pdb_pos] = plddt;
                    }

                    let mut hash_vec = get_geometric_hash_as_u32_from_structure(
                        &compact, self.hash_type,
                        self.num_bin_dist, self.num_bin_angle,
                        self.dist_cutoff, &self.multiple_bins,
                    );
                    // Drop intermediate variables
                    drop(compact);

                    // If remove_redundancy is true, remove duplicates
                    hash_vec.sort_unstable();
                    hash_vec.dedup();
                    hash_vec.iter().map(|x| (*x, pdb_pos)).collect()
                }).flatten().collect()
            });
        self.hash_id_vec = collected;
        self.nres_vec = shared_data.get_nres_vec_clone();
        self.plddt_vec = shared_data.get_plddt_vec_clone();
        drop(pool);
    }
    
    pub fn collect_and_count(&mut self) {
        let shared_data = SharedData::new(self.path_vec.len());
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating hashes");
        // For iterating files, apply multi-threading with num_threads_for_file
        let chunk_size = self.num_threads * 128;
        // Chunk pdb paths
        let chunked_paths = self.path_vec.chunks(chunk_size);
        let total_chunks = chunked_paths.len();
        chunked_paths.enumerate().for_each(|(chunk_index, chunk)| {
            // Print percentage of completion
            print_log_msg(INFO, &format!("Processing chunk {}/{}", chunk_index + 1, total_chunks));
            let collected: Vec<(u32, usize)> = pool.install(|| {
                chunk
                    .par_iter()
                    .enumerate()
                    .map(|(local_pos, pdb_path)| {
                        let pdb_pos = chunk_index * chunk_size + local_pos;
                        // let pdb_pos = self.path_vec.iter().position(|x| x == pdb_path).unwrap();

                        #[cfg(not(feature = "foldcomp"))]
                        let compact = read_structure_from_path(pdb_path).expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        );
                        #[cfg(feature = "foldcomp")]
                        let compact = if self.is_foldcomp_enabled {
                            self.foldcomp_db_reader.read_single_structure_by_id(
                                self.numeric_db_key_vec[pdb_pos]
                            ).expect(
                                log_msg(FAIL, "Failed to read structure").as_str()
                            )
                        } else {
                            read_structure_from_path(pdb_path).expect(
                                log_msg(FAIL, "Failed to read structure").as_str()
                            )
                        };
                        if compact.num_residues > self.max_residue {
                            print_log_msg(WARN, &format!("{} has too many residues. Skipping", pdb_path));
                            // skip this file
                            // Drop intermediate variables
                            drop(compact);
                            return Vec::new();
                        }
                        let compact = compact.to_compact();
                        // Directly write num_residues and avg_plddt to the vectors
                        let nres = compact.num_residues;
                        let plddt = compact.get_avg_plddt();
                        // Mutex free version
                        unsafe {
                            let nres_vec = shared_data.get_nres_vec();
                            let plddt_vec = shared_data.get_plddt_vec();
                            let nres_vec = &mut *nres_vec.get();
                            let plddt_vec = &mut *plddt_vec.get();
                            nres_vec[pdb_pos] = nres;
                            plddt_vec[pdb_pos] = plddt;
                        }

                        let mut hash_vec = get_geometric_hash_as_u32_from_structure(
                            &compact, self.hash_type, 
                            self.num_bin_dist, self.num_bin_angle,
                            self.dist_cutoff, &self.multiple_bins,
                        );
                        // Drop intermediate variables
                        drop(compact);
                        // If remove_redundancy is true, remove duplicates

                        hash_vec.sort_unstable();
                        hash_vec.dedup();
                        hash_vec.iter().map(|x| (x.clone(), pdb_pos)).collect()

                    }).flatten().collect()
                });
            pool.install(|| {
                (0..self.num_threads).into_par_iter().for_each(| tid | {
                    // Thread only saves hashes with same modulos
                    let _ = &collected.iter().for_each(|(hash, pdb_pos)| {
                        if hash % self.num_threads as u32 == tid as u32 {
                            self.fold_disco_index.count_single_entry(*hash, *pdb_pos);
                        }
                    });
                });
            });
            drop(collected);
        });

        self.nres_vec = shared_data.get_nres_vec_clone();
        self.plddt_vec = shared_data.get_plddt_vec_clone();
        drop(pool);
    }
    
    pub fn add_entries(&mut self) {
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating hashes");
        // For iterating files, apply multi-threading with num_threads_for_file
        let chunk_size = self.num_threads * 128;
        // Chunk pdb paths
        let chunked_paths = self.path_vec.chunks(chunk_size);
        let total_chunks = chunked_paths.len();
        chunked_paths.enumerate().for_each(|(chunk_index, chunk)| {
            print_log_msg(INFO, &format!("Processing chunk {}/{}", chunk_index + 1, total_chunks));
            let collected: Vec<(u32, usize)> = pool.install(|| {
                chunk
                    .par_iter()
                    .enumerate()
                    .map(|(local_pos, pdb_path)| {
                        let pdb_pos = chunk_index * chunk_size + local_pos;

                        #[cfg(not(feature = "foldcomp"))]
                        let compact = read_structure_from_path(pdb_path).expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        );

                        #[cfg(feature = "foldcomp")]
                        let compact = if self.is_foldcomp_enabled {
                            self.foldcomp_db_reader.read_single_structure_by_id(
                                self.numeric_db_key_vec[pdb_pos]
                            ).expect(
                                log_msg(FAIL, "Failed to read structure").as_str()
                            )
                        } else {
                            read_structure_from_path(pdb_path).expect(
                                log_msg(FAIL, "Failed to read structure").as_str()
                            )
                        };

                        if compact.num_residues > self.max_residue {
                            print_log_msg(WARN, &format!("{} has too many residues. Skipping", pdb_path));
                            // skip this file
                            // Drop intermediate variables
                            drop(compact);
                            return Vec::new();
                        }
                        let compact = compact.to_compact();
                        // Directly write num_residues and avg_plddt to the vectors
                        let mut hash_vec = get_geometric_hash_as_u32_from_structure(
                            &compact, self.hash_type, self.num_bin_dist, self.num_bin_angle,
                            self.dist_cutoff, &self.multiple_bins,
                        );
                        // Drop intermediate variables
                        drop(compact);
                        // If remove_redundancy is true, remove duplicates
                        hash_vec.sort_unstable();
                        hash_vec.dedup();
                        hash_vec.iter().map(|x| (x.clone(), pdb_pos)).collect()

                    }).flatten().collect()
                });
            pool.install(|| {
                (0..self.num_threads).into_par_iter().for_each(| tid | {
                    // Thread only saves hashes with same modulos
                    let mut bit_containers = Vec::with_capacity(8);
                    let _ = &collected.iter().for_each(|(hash, pdb_pos)| {
                        if hash % self.num_threads as u32 == tid as u32 {
                            self.fold_disco_index.add_single_entry(*hash, *pdb_pos, &mut bit_containers);
                        }
                    });
                });
            });
            drop(collected);
        });
        drop(pool);
    }
    
    pub fn sort_hash_vec(&mut self) {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect(&log_msg(FAIL, "Failed to build thread pool for sorting"));
        pool.install(|| {
            self.hash_id_vec.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
        });
        drop(pool);
    }
    
    pub fn get_allocation_size(&self) -> usize {
        let mut allocation_size = 0usize;
        self.path_vec.iter().for_each(|pdb_path| {
            let pdb_reader = PDBReader::from_file(pdb_path).expect(
            log_msg(FAIL, "PDB file not found").as_str()
            );
            let compact = pdb_reader.read_structure().expect(
                log_msg(FAIL, "Failed to read structure").as_str()
            );
            allocation_size += compact.num_residues * (compact.num_residues - 1);
        });
        allocation_size
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
}

fn string_vec_to_numeric_id_vec(string_vec: &Vec<String>, numeric_id_vec: &mut Vec<usize>) {
    for i in 0..string_vec.len() {
        numeric_id_vec.push(i);
    }
}

#[allow(dead_code)]
fn numeric_id_vec_from_string_vec(string_vec: &Vec<String>) -> Vec<usize> {
    let mut numeric_id_vec = Vec::with_capacity(string_vec.len());
    for i in 0..string_vec.len() {
        numeric_id_vec.push(i);
    }
    numeric_id_vec
}

struct SharedData {
    nres_vec: Arc<UnsafeCell<Vec<usize>>>,
    plddt_vec: Arc<UnsafeCell<Vec<f32>>>,
}

unsafe impl Sync for SharedData {}

impl SharedData {
    fn new(size: usize) -> Self {
        SharedData {
            nres_vec: Arc::new(UnsafeCell::new(vec![0; size])),
            plddt_vec: Arc::new(UnsafeCell::new(vec![0.0; size])),
        }
    }

    fn get_nres_vec(&self) -> &UnsafeCell<Vec<usize>> {
        &self.nres_vec
    }

    fn get_plddt_vec(&self) -> &UnsafeCell<Vec<f32>> {
        &self.plddt_vec
    }
    
    fn get_nres_vec_clone(&self) -> Vec<usize> {
        unsafe {
            (*self.nres_vec.get()).clone()
        }
    }
    
    fn get_plddt_vec_clone(&self) -> Vec<f32> {
        unsafe {
            (*self.plddt_vec.get()).clone()
        }
    }
}