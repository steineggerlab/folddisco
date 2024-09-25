
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