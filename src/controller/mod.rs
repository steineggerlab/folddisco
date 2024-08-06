// File: mod.rs
// Created: 2024-01-18 15:47:16
// Description:
//    a new controller implementation that supports multiple hash types
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

pub mod feature;
pub mod graph;
pub mod io;
pub mod query;
pub mod retrieve;
pub mod rank;
pub mod map;
pub mod mode;

use std::io::Write;
use std::sync::{Arc, Mutex};
use feature::get_geometric_hash_as_u32_from_structure;
use mode::IndexMode;
// External imports
use rayon::{current_thread_index, prelude::*};

use crate::index::indextable::FolddiscoIndex;
// Internal imports
use crate::prelude::print_log_msg;
use crate::structure::grid::DEFAULT_GRID_WIDTH;
use crate::{measure_time, PDBReader};
use crate::geometry::core::{GeometricHash, HashType};
use crate::index::alloc::IndexBuilder;
use crate::utils::log::{ log_msg, FAIL, WARN };
use crate::controller::feature::get_geometric_hash_from_structure;

use self::feature::get_geometric_hash_with_grid;

#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::{FoldcompDbReader, get_id_vector_subset_out_of_lookup, get_name_vector_subset_out_of_lookup};

const DEFAULT_REMOVE_REDUNDANCY: bool = true;
const DEFAULT_NUM_THREADS: usize = 4;
const DEFAULT_HASH_TYPE: HashType = HashType::PDBTrRosetta;
const DEFAULT_MAX_RESIDUE: usize = 50000;

pub struct FoldDisco {
    pub path_vec: Vec<String>,
    pub numeric_id_vec: Vec<usize>,
    pub nres_vec: Vec<usize>,
    pub plddt_vec: Vec<f32>,
    pub hash_collection: Vec<Vec<GeometricHash>>,
    pub hash_id_pairs: Vec<(GeometricHash, usize)>,
    pub hash_id_grids: Vec<(GeometricHash, usize)>,
    pub hash_id_vec: Vec<(u32, usize)>,
    pub index_builder: IndexBuilder<usize, GeometricHash>,
    pub hash_type: HashType,
    pub remove_redundancy: bool,
    pub num_threads: usize,
    pub num_bin_dist: usize,
    pub num_bin_angle: usize,
    pub output_path: String,
    pub max_residue: usize,
    pub grid_width: f32,
    pub index_mode: IndexMode,
    pub fold_disco_index: FolddiscoIndex,
    pub foldcomp_db_path: String,
    #[cfg(not(feature = "foldcomp"))]
    pub foldcomp_db_reader: bool,
    #[cfg(feature = "foldcomp")]
    pub foldcomp_db_reader: FoldcompDbReader,
    pub is_foldcomp_enabled: bool,
}

impl FoldDisco {
    pub fn new(path_vec: Vec<String>) -> FoldDisco {
        let length = path_vec.len();
        FoldDisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::with_capacity(length),
            nres_vec: Vec::with_capacity(length),
            plddt_vec: Vec::with_capacity(length),
            hash_collection: Vec::new(),
            hash_id_pairs: Vec::new(),
            hash_id_grids: Vec::new(),
            hash_id_vec: Vec::new(),
            index_builder: IndexBuilder::empty(),
            hash_type: DEFAULT_HASH_TYPE,
            remove_redundancy: DEFAULT_REMOVE_REDUNDANCY,
            num_threads: DEFAULT_NUM_THREADS,
            num_bin_dist: 0,
            num_bin_angle: 0,
            output_path: String::new(),
            max_residue: DEFAULT_MAX_RESIDUE,
            grid_width: DEFAULT_GRID_WIDTH,
            index_mode: IndexMode::Id,
            fold_disco_index: FolddiscoIndex::new(0usize, String::new()),
            // By default, foldcomp is not enabled. foldcomp reading is enabled within new_with_params
            foldcomp_db_path: String::new(),
            #[cfg(not(feature = "foldcomp"))]
            foldcomp_db_reader: false,
            #[cfg(feature = "foldcomp")]
            foldcomp_db_reader: FoldcompDbReader::empty(),
            is_foldcomp_enabled: false,
        }
    }
    pub fn new_with_hash_type(path_vec: Vec<String>, hash_type: HashType) -> FoldDisco {
        let length = path_vec.len();
        FoldDisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::with_capacity(length),
            nres_vec: Vec::with_capacity(length),
            plddt_vec: Vec::with_capacity(length),
            hash_collection: Vec::new(),
            hash_id_pairs: Vec::new(),
            hash_id_grids: Vec::new(),
            hash_id_vec: Vec::new(),
            index_builder: IndexBuilder::empty(),
            hash_type: hash_type,
            remove_redundancy: DEFAULT_REMOVE_REDUNDANCY,
            num_threads: DEFAULT_NUM_THREADS,
            num_bin_dist: 0,
            num_bin_angle: 0,
            output_path: String::new(),
            max_residue: DEFAULT_MAX_RESIDUE,
            grid_width: DEFAULT_GRID_WIDTH,
            index_mode: IndexMode::Id,
            fold_disco_index: FolddiscoIndex::new(0usize, String::new()),
            foldcomp_db_path: String::new(),
            #[cfg(not(feature = "foldcomp"))]
            foldcomp_db_reader: false,
            #[cfg(feature = "foldcomp")]
            foldcomp_db_reader: FoldcompDbReader::empty(),
            is_foldcomp_enabled: false,
        }
    }

    pub fn new_with_params(
        path_vec: Vec<String>, hash_type: HashType, remove_redundancy: bool,
        num_threads: usize, num_bin_dist: usize, num_bin_angle: usize, output_path: String,
        grid_width: f32, index_mode: IndexMode
    ) -> FoldDisco {
        let length = path_vec.len();
        let total_hashes = match index_mode {
            IndexMode::Id | IndexMode::Grid | IndexMode::Pos => 0,
            IndexMode::Big => 2usize.pow(hash_type.encoding_bits() as u32),
        };
        FoldDisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::with_capacity(length),
            nres_vec: Vec::with_capacity(length),
            plddt_vec: Vec::with_capacity(length),
            hash_collection: Vec::new(),
            hash_id_pairs: Vec::new(),
            hash_id_grids: Vec::new(),
            hash_id_vec: Vec::new(),
            index_builder: IndexBuilder::empty(),
            hash_type: hash_type,
            remove_redundancy: remove_redundancy,
            num_threads: num_threads,
            num_bin_dist: num_bin_dist,
            num_bin_angle: num_bin_angle,
            output_path: output_path.clone(),
            max_residue: DEFAULT_MAX_RESIDUE,
            grid_width: grid_width,
            index_mode: index_mode,
            fold_disco_index: FolddiscoIndex::new(total_hashes, output_path.clone()),
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
        path_vec: Vec<String>, hash_type: HashType, remove_redundancy: bool,
        num_threads: usize, num_bin_dist: usize, num_bin_angle: usize, output_path: String,
        grid_width: f32, index_mode: IndexMode, foldcomp_db_path: &'static str
    ) -> FoldDisco {
        let length = path_vec.len();
        let total_hashes = match index_mode {
            IndexMode::Id | IndexMode::Grid | IndexMode::Pos => 0,
            IndexMode::Big => 2usize.pow(hash_type.encoding_bits() as u32),
        };
        let mut foldcomp_db_reader = FoldcompDbReader::new(&foldcomp_db_path);
        
        FoldDisco {
            path_vec: path_vec,
            numeric_id_vec: Vec::with_capacity(length),
            nres_vec: Vec::with_capacity(length),
            plddt_vec: Vec::with_capacity(length),
            hash_collection: Vec::new(),
            hash_id_pairs: Vec::new(),
            hash_id_grids: Vec::new(),
            hash_id_vec: Vec::new(),
            index_builder: IndexBuilder::empty(),
            hash_type: hash_type,
            remove_redundancy: remove_redundancy,
            num_threads: num_threads,
            num_bin_dist: num_bin_dist,
            num_bin_angle: num_bin_angle,
            output_path: output_path.clone(),
            max_residue: DEFAULT_MAX_RESIDUE,
            grid_width: grid_width,
            index_mode: index_mode,
            fold_disco_index: FolddiscoIndex::new(total_hashes, output_path.clone()),
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
    pub fn set_remove_redundancy(&mut self, remove_redundancy: bool) {
        self.remove_redundancy = remove_redundancy;
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

    // #[deprecated]
    pub fn collect_hash(&mut self) {
        let path_order: Arc<Vec<usize>> = Arc::new(Vec::with_capacity(self.path_vec.len()));
        let nres_vec: Arc<Mutex<Vec<usize>>> = Arc::new(Mutex::new(vec![0; self.path_vec.len()]));
        let plddt_vec: Arc<Mutex<Vec<f32>>> = Arc::new(Mutex::new(vec![0.0; self.path_vec.len()]));
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating files");
        // For iterating files, apply multi-threading with num_threads_for_file
        let (output, positions): (Vec<Vec<GeometricHash>>, Vec<usize>) = pool.install(|| {
            self.path_vec
                .par_iter()
                .map(|pdb_path| {
                    let _path_order = path_order.clone();
                    // let mut path_order = Arc::make_mut(&mut path_order);
                    let pdb_pos = self.path_vec.iter().position(|x| x == pdb_path).unwrap();
                    let pdb_reader = PDBReader::from_file(pdb_path).expect(
                        log_msg(FAIL, "PDB file not found").as_str()
                    );
                    let compact = pdb_reader.read_structure().expect(
                        log_msg(FAIL, "Failed to read structure").as_str()
                    );
                    let compact = compact.to_compact();
                    // Directly write num_residues and avg_plddt to the vectors
                    let nres = compact.num_residues;
                    let plddt = compact.get_avg_plddt();
                    {
                        let mut nres_vec = nres_vec.lock().unwrap();
                        let mut plddt_vec = plddt_vec.lock().unwrap();
                        nres_vec[pdb_pos] = nres;
                        plddt_vec[pdb_pos] = plddt;
                    }
                    // Skip if the number of residues is too large
                    if compact.num_residues > self.max_residue {
                        print_log_msg(WARN, &format!("{} has too many residues. Skipping", pdb_path));
                        // skip this file
                        return (Vec::new(), pdb_pos);
                    }
 
                    let hash_vec = get_geometric_hash_from_structure(
                        &compact, self.hash_type, self.num_bin_dist, self.num_bin_angle
                    );
                    // Drop intermediate variables
                    drop(compact);
                    drop(pdb_reader);
                    // If remove_redundancy is true, remove duplicates
                    if self.remove_redundancy {
                        let mut hash_vec = hash_vec;
                        hash_vec.sort_unstable();
                        hash_vec.dedup();
                        (hash_vec, pdb_pos)
                    } else {
                        (hash_vec, pdb_pos)
                    }
                })
            }).unzip();
        self.hash_collection = output;
        self.numeric_id_vec = positions;
        // Flatten Arc<Mutex<Vec<T>>> to Vec<T>
        self.nres_vec = Arc::try_unwrap(nres_vec).unwrap().into_inner().unwrap();
        self.plddt_vec = Arc::try_unwrap(plddt_vec).unwrap().into_inner().unwrap();
        // Remove  thread pool
        drop(pool);
    }
    
    pub fn collect_hash_pairs(&mut self) {
        let nres_vec: Arc<Mutex<Vec<usize>>> = Arc::new(Mutex::new(vec![0; self.path_vec.len()]));
        let plddt_vec: Arc<Mutex<Vec<f32>>> = Arc::new(Mutex::new(vec![0.0; self.path_vec.len()]));
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating files");
        // For iterating files, apply multi-threading with num_threads_for_file
        let collected: Vec<(GeometricHash, usize)> = pool.install(|| {
            self.path_vec
                .par_iter()
                .map(|pdb_path| {
                    let pdb_pos = self.path_vec.iter().position(|x| x == pdb_path).unwrap();
                    let pdb_reader = PDBReader::from_file(pdb_path).expect(
                        log_msg(FAIL, "PDB file not found").as_str()
                    );
                    let compact = if pdb_path.ends_with(".gz") {
                        pdb_reader.read_structure_from_gz().expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        )
                    } else {
                        pdb_reader.read_structure().expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        )
                    };
                    if compact.num_residues > self.max_residue {
                        print_log_msg(WARN, &format!("{} has too many residues. Skipping", pdb_path));
                        // skip this file
                        // Drop intermediate variables
                        drop(compact);
                        drop(pdb_reader);
                        return Vec::new();
                    }
                    let compact = compact.to_compact();
                    // Directly write num_residues and avg_plddt to the vectors
                    let nres = compact.num_residues;
                    let plddt = compact.get_avg_plddt();
                    {
                        let mut nres_vec = nres_vec.lock().unwrap();
                        let mut plddt_vec = plddt_vec.lock().unwrap();
                        nres_vec[pdb_pos] = nres;
                        plddt_vec[pdb_pos] = plddt;
                    }
                    let mut hash_vec = get_geometric_hash_from_structure(
                        &compact, self.hash_type, self.num_bin_dist, self.num_bin_angle
                    );
                    // Drop intermediate variables
                    drop(compact);
                    drop(pdb_reader);
                    // If remove_redundancy is true, remove duplicates
                    if self.remove_redundancy {
                        hash_vec.sort_unstable();
                        hash_vec.dedup();
                        hash_vec.iter().map(|x| (*x, pdb_pos)).collect()
                    } else {
                        hash_vec.iter().map(|x| (*x, pdb_pos)).collect()
                    }
                }).flatten().collect()
            });
        self.hash_id_pairs = collected;
        // self.hash_id_grids = collected;
        self.nres_vec = Arc::try_unwrap(nres_vec).unwrap().into_inner().unwrap();
        self.plddt_vec = Arc::try_unwrap(plddt_vec).unwrap().into_inner().unwrap();
        drop(pool);
    }

    pub fn collect_hash_grids(&mut self) {
        let nres_vec: Arc<Mutex<Vec<usize>>> = Arc::new(Mutex::new(vec![0; self.path_vec.len()]));
        let plddt_vec: Arc<Mutex<Vec<f32>>> = Arc::new(Mutex::new(vec![0.0; self.path_vec.len()]));
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating files");
        // For iterating files, apply multi-threading with num_threads_for_file
        let collected: Vec<(GeometricHash, usize)> = pool.install(|| {
            self.path_vec
                .par_iter()
                .map(|pdb_path| {
                    let pdb_pos = self.path_vec.iter().position(|x| x == pdb_path).unwrap();
                    let pdb_reader = PDBReader::from_file(pdb_path).expect(
                        log_msg(FAIL, "PDB file not found").as_str()
                    );
                    let compact = if pdb_path.ends_with(".gz") {
                        pdb_reader.read_structure_from_gz().expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        )
                    } else {
                        pdb_reader.read_structure().expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        )
                    };
                    if compact.num_residues > self.max_residue {
                        print_log_msg(WARN, &format!("{} has too many residues. Skipping", pdb_path));
                        // skip this file
                        // Drop intermediate variables
                        drop(compact);
                        drop(pdb_reader);
                        return Vec::new();
                    }
                    let compact = compact.to_compact();
                    // Directly write num_residues and avg_plddt to the vectors
                    let nres = compact.num_residues;
                    let plddt = compact.get_avg_plddt();
                    {
                        let mut nres_vec = nres_vec.lock().unwrap();
                        let mut plddt_vec = plddt_vec.lock().unwrap();
                        nres_vec[pdb_pos] = nres;
                        plddt_vec[pdb_pos] = plddt;
                    }
                    // let mut hash_vec = get_geometric_hash_from_structure(
                    //     &compact, self.hash_type, self.num_bin_dist, self.num_bin_angle
                    // );
                    let mut hash_with_grid = get_geometric_hash_with_grid(
                        &compact, pdb_pos, self.hash_type, 
                        self.num_bin_dist, self.num_bin_angle, self.grid_width
                    );

                    // Drop intermediate variables
                    drop(compact);
                    drop(pdb_reader);
                    // If remove_redundancy is true, remove duplicates
                    if self.remove_redundancy {
                        hash_with_grid.sort_unstable();
                        hash_with_grid.dedup();
                        hash_with_grid.iter().map(|x| *x).collect()
                    } else {
                        hash_with_grid.iter().map(|x| *x).collect()
                    }
                }).flatten().collect()
            });
        // self.hash_id_pairs = collected;
        self.hash_id_grids = collected;
        self.nres_vec = Arc::try_unwrap(nres_vec).unwrap().into_inner().unwrap();
        self.plddt_vec = Arc::try_unwrap(plddt_vec).unwrap().into_inner().unwrap();
        drop(pool);
    }

    pub fn collect_hash_vec(&mut self) {
        let nres_vec: Arc<Mutex<Vec<usize>>> = Arc::new(Mutex::new(vec![0; self.path_vec.len()]));
        let plddt_vec: Arc<Mutex<Vec<f32>>> = Arc::new(Mutex::new(vec![0.0; self.path_vec.len()]));
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating files");
        // For iterating files, apply multi-threading with num_threads_for_file
        let collected: Vec<(u32, usize)> = pool.install(|| {
            self.path_vec
                .par_iter()
                .map(|pdb_path| {
                    let pdb_pos = self.path_vec.iter().position(|x| x == pdb_path).unwrap();
                    #[cfg(not(feature = "foldcomp"))]
                    let pdb_reader = PDBReader::from_file(pdb_path).expect(
                        log_msg(FAIL, "PDB file not found").as_str()
                    );
                    #[cfg(not(feature = "foldcomp"))]
                    let compact = if pdb_path.ends_with(".gz") {
                        pdb_reader.read_structure_from_gz().expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        )
                    } else {
                        pdb_reader.read_structure().expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        )
                    };

                    #[cfg(feature = "foldcomp")]
                    let compact = if self.is_foldcomp_enabled {
                        self.foldcomp_db_reader.read_single_structure(pdb_path).expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        )
                    } else {
                        let pdb_reader = PDBReader::from_file(pdb_path).expect(
                            log_msg(FAIL, "PDB file not found").as_str()
                        );
                        if pdb_path.ends_with(".gz") {
                            pdb_reader.read_structure_from_gz().expect(
                                log_msg(FAIL, "Failed to read structure").as_str()
                            )
                        } else {
                            pdb_reader.read_structure().expect(
                                log_msg(FAIL, "Failed to read structure").as_str()
                            )
                        }
                    };  

                    if compact.num_residues > self.max_residue {
                        print_log_msg(WARN, &format!("{} has too many residues. Skipping", pdb_path));
                        // skip this file
                        // Drop intermediate variables
                        drop(compact);
                        #[cfg(not(feature = "foldcomp"))]
                        drop(pdb_reader);
                        return Vec::new();
                    }
                    let compact = compact.to_compact();
                    // Directly write num_residues and avg_plddt to the vectors
                    let nres = compact.num_residues;
                    let plddt = compact.get_avg_plddt();
                    {
                        let mut nres_vec = nres_vec.lock().unwrap();
                        let mut plddt_vec = plddt_vec.lock().unwrap();
                        nres_vec[pdb_pos] = nres;
                        plddt_vec[pdb_pos] = plddt;
                    }
                    let mut hash_vec = get_geometric_hash_as_u32_from_structure(
                        &compact, self.hash_type, self.num_bin_dist, self.num_bin_angle
                    );
                    // Drop intermediate variables
                    drop(compact);
                    #[cfg(not(feature = "foldcomp"))]
                    drop(pdb_reader);
                    // If remove_redundancy is true, remove duplicates
                    if self.remove_redundancy {
                        hash_vec.sort_unstable();
                        hash_vec.dedup();
                        hash_vec.iter().map(|x| (*x, pdb_pos)).collect()
                    } else {
                        hash_vec.iter().map(|x| (*x, pdb_pos)).collect()
                    }
                }).flatten().collect()
            });
        self.hash_id_vec = collected;
        // self.hash_id_grids = collected;
        self.nres_vec = Arc::try_unwrap(nres_vec).unwrap().into_inner().unwrap();
        self.plddt_vec = Arc::try_unwrap(plddt_vec).unwrap().into_inner().unwrap();
        drop(pool);
    }
    
    pub fn collect_and_count(&mut self) {
        let nres_vec: Arc<Mutex<Vec<usize>>> = Arc::new(Mutex::new(vec![0; self.path_vec.len()]));
        let plddt_vec: Arc<Mutex<Vec<f32>>> = Arc::new(Mutex::new(vec![0.0; self.path_vec.len()]));
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating hashes");
        // For iterating files, apply multi-threading with num_threads_for_file
        let chunk_size = 65536;
        // Chunk pdb paths
        let chunked_paths = self.path_vec.chunks(chunk_size);
        chunked_paths.for_each(|chunk| {
            let collected: Vec<(u32, usize)> = pool.install(|| {
                chunk
                    .par_iter()
                    .map(|pdb_path| {
                        let pdb_pos = self.path_vec.iter().position(|x| x == pdb_path).unwrap();
                        #[cfg(not(feature = "foldcomp"))]
                        let pdb_reader = PDBReader::from_file(pdb_path).expect(
                            log_msg(FAIL, "PDB file not found").as_str()
                        );
                        #[cfg(not(feature = "foldcomp"))]
                        let compact = if pdb_path.ends_with(".gz") {
                            pdb_reader.read_structure_from_gz().expect(
                                log_msg(FAIL, "Failed to read structure").as_str()
                            )
                        } else {
                            pdb_reader.read_structure().expect(
                                log_msg(FAIL, "Failed to read structure").as_str()
                            )
                        };
                        #[cfg(feature = "foldcomp")]
                        let compact = self.foldcomp_db_reader.read_single_structure(pdb_path).expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        );

                        if compact.num_residues > self.max_residue {
                            print_log_msg(WARN, &format!("{} has too many residues. Skipping", pdb_path));
                            // skip this file
                            // Drop intermediate variables
                            drop(compact);
                            #[cfg(not(feature = "foldcomp"))]
                            drop(pdb_reader);
                            return Vec::new();
                        }
                        let compact = compact.to_compact();
                        // Directly write num_residues and avg_plddt to the vectors
                        let nres = compact.num_residues;
                        let plddt = compact.get_avg_plddt();
                        {
                            let mut nres_vec = nres_vec.lock().unwrap();
                            let mut plddt_vec = plddt_vec.lock().unwrap();
                            nres_vec[pdb_pos] = nres;
                            plddt_vec[pdb_pos] = plddt;
                        }
                        let mut hash_vec = get_geometric_hash_from_structure(
                            &compact, self.hash_type, self.num_bin_dist, self.num_bin_angle
                        );
                        // Drop intermediate variables
                        drop(compact);
                        #[cfg(not(feature = "foldcomp"))]
                        drop(pdb_reader);
                        // If remove_redundancy is true, remove duplicates
                        if self.remove_redundancy {
                            hash_vec.sort_unstable();
                            hash_vec.dedup();
                            hash_vec.iter().map(|x| (x.as_u32(), pdb_pos)).collect()
                        } else {
                            hash_vec.iter().map(|x| (x.as_u32(), pdb_pos)).collect()
                        }
                    }).flatten().collect()
                });
            pool.install(|| {
                (0..self.num_threads).into_par_iter().for_each(| tid | {
                    // Thread only saves hashes with same modulos
                    &collected.iter().for_each(|(hash, pdb_pos)| {
                        if hash % self.num_threads as u32 == tid as u32 {
                            self.fold_disco_index.count_single_entry(*hash, *pdb_pos);
                        }
                    });
                });
            });
            drop(collected);
        });
        // self.hash_id_grids = collected;
        self.nres_vec = Arc::try_unwrap(nres_vec).unwrap().into_inner().unwrap();
        self.plddt_vec = Arc::try_unwrap(plddt_vec).unwrap().into_inner().unwrap();
        drop(pool);
    }
    
    pub fn add_entries(&mut self) {
        // Set file threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect("Failed to build thread pool for iterating hashes");
        // For iterating files, apply multi-threading with num_threads_for_file
        let chunk_size = 65536;
        // Chunk pdb paths
        let chunked_paths = self.path_vec.chunks(chunk_size);
        chunked_paths.for_each(|chunk| {
            let collected: Vec<(u32, usize)> = pool.install(|| {
                chunk
                    .par_iter()
                    .map(|pdb_path| {
                        let pdb_pos = self.path_vec.iter().position(|x| x == pdb_path).unwrap();
                        #[cfg(not(feature = "foldcomp"))]
                        let pdb_reader = PDBReader::from_file(pdb_path).expect(
                            log_msg(FAIL, "PDB file not found").as_str()
                        );
                        #[cfg(not(feature = "foldcomp"))]
                        let compact = if pdb_path.ends_with(".gz") {
                            pdb_reader.read_structure_from_gz().expect(
                                log_msg(FAIL, "Failed to read structure").as_str()
                            )
                        } else {
                            pdb_reader.read_structure().expect(
                                log_msg(FAIL, "Failed to read structure").as_str()
                            )
                        };

                        #[cfg(feature = "foldcomp")]
                        let compact = self.foldcomp_db_reader.read_single_structure(pdb_path).expect(
                            log_msg(FAIL, "Failed to read structure").as_str()
                        );

                        if compact.num_residues > self.max_residue {
                            print_log_msg(WARN, &format!("{} has too many residues. Skipping", pdb_path));
                            // skip this file
                            // Drop intermediate variables
                            drop(compact);
                            #[cfg(not(feature = "foldcomp"))]
                            drop(pdb_reader);
                            return Vec::new();
                        }
                        let compact = compact.to_compact();
                        // Directly write num_residues and avg_plddt to the vectors
                        let mut hash_vec = get_geometric_hash_from_structure(
                            &compact, self.hash_type, self.num_bin_dist, self.num_bin_angle
                        );
                        // Drop intermediate variables
                        drop(compact);
                        #[cfg(not(feature = "foldcomp"))]
                        drop(pdb_reader);
                        // If remove_redundancy is true, remove duplicates
                        if self.remove_redundancy {
                            hash_vec.sort_unstable();
                            hash_vec.dedup();
                            hash_vec.iter().map(|x| (x.as_u32(), pdb_pos)).collect()
                        } else {
                            hash_vec.iter().map(|x| (x.as_u32(), pdb_pos)).collect()
                        }
                    }).flatten().collect()
                });
            pool.install(|| {
                (0..self.num_threads).into_par_iter().for_each(| tid | {
                    // Thread only saves hashes with same modulos
                    &collected.iter().for_each(|(hash, pdb_pos)| {
                        if hash % self.num_threads as u32 == tid as u32 {
                            self.fold_disco_index.add_single_entry(*hash, *pdb_pos);
                        }
                    });
                });
            });
            drop(collected);
        });
        drop(pool);
    }
    
    
    pub fn sort_hash_pairs(&mut self) {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect(&log_msg(FAIL, "Failed to build thread pool for sorting"));
        pool.install(|| {
            self.hash_id_pairs.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
        });
        drop(pool);
    }
    
    pub fn sort_hash_grids(&mut self) {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .expect(&log_msg(FAIL, "Failed to build thread pool for sorting"));
        pool.install(|| {
            self.hash_id_grids.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
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

    pub fn save_raw_feature(&mut self, _path: &str, _discretize: bool) {
        todo!("Implement save_raw_feature method!!!");
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
    
    pub fn return_hash_collection(&mut self) -> Vec<Vec<GeometricHash>> {
        // Move self.hash_collection and clear it
        let hash_collection = std::mem::take(&mut self.hash_collection);
        hash_collection
    }

    pub fn set_index_table(&mut self) {
        let index_builder = IndexBuilder::new(
            self.numeric_id_vec.clone(), self.return_hash_collection(),
            self.num_threads, 1,
            format!("{}.offset", self.output_path), // offset file path
            format!("{}.index", self.output_path), // data file path
        );
        self.index_builder = index_builder;
    }
    
    pub fn fill_index_table(&mut self) {
        // self.index_builder.fill_with_dashmap();
        let index_map = measure_time!(self.index_builder.fill_and_return_dashmap());
        // TODO: IMPORTANT: Move this.
        let (_offset_map, _value_vec) = measure_time!(self.index_builder.convert_hashmap_to_offset_and_values(index_map));
    }
    
    pub fn save_offset_map(&self) {
        todo!("Implement save_offset_map method");
        // self.index_builder.save_offset_map();
    }
    
    pub fn save_index_table(&self) {
        // self.index_builder.save();
        todo!("Implement save_index_table method");
    }
    

}

fn string_vec_to_numeric_id_vec(string_vec: &Vec<String>, numeric_id_vec: &mut Vec<usize>) {
    for i in 0..string_vec.len() {
        numeric_id_vec.push(i);
    }
}

fn numeric_id_vec_from_string_vec(string_vec: &Vec<String>) -> Vec<usize> {
    let mut numeric_id_vec = Vec::with_capacity(string_vec.len());
    for i in 0..string_vec.len() {
        numeric_id_vec.push(i);
    }
    numeric_id_vec
}