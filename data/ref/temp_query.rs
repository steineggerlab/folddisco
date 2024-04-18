
pub fn query_pdb(env: AppArgs) {
    match env {
        AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            retrieve,
            amino_acid, // TODO:" Implement amino acid mode"
            dist_threshold,
            angle_threshold,
            match_cutoff,
            score_cutoff,
            num_res_cutoff,
            plddt_cutoff,
            verbose,
            help,
        } => {
            if help {
                eprintln!("{}", HELP_QUERY);
                std::process::exit(0);
            }
            // Check if arguments are valid
            if index_path.is_none() {
                eprintln!("{}", HELP_QUERY);
                std::process::exit(1);
            }
            // Print query information
            if verbose {
                print_log_msg(INFO, &format!("Querying {}:{} to {}", &pdb_path, &query_string, &index_path.clone().unwrap()));
                // NOTE: If needed, print filter information
            }
            // Get index paths
            let index_paths = check_and_get_indices(index_path.clone(), verbose);
            if verbose {
                print_log_msg(INFO, &format!("Found {} index file(s). Querying with {} threads", index_paths.len(), threads));
            }
            // Set thread pool
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();

            let queries = if query_string.ends_with(".txt") || query_string.ends_with(".tsv") {
                // Read file and get path, query, output by line
                let mut queries: Vec<(String, String, String)> = Vec::new();
                let file = std::fs::File::open(&query_string).expect(
                    &log_msg(FAIL, &format!("Failed to open query file: {}", &query_string))
                );
                let reader = std::io::BufReader::new(file);
                for line in reader.lines() {
                    let line = line.expect("Failed to read line");
                    let mut split = line.split('\t');
                    let pdb_path = split.next().expect("Failed to get pdb path").to_string();
                    let query_string = split.next().unwrap_or("").to_string();
                    let output_path = split.next().unwrap_or("").to_string();
                    queries.push((pdb_path, query_string, output_path));
                }
                queries
            } else {
                vec![(pdb_path.clone(), query_string.clone(), String::new())]
            };
            let dist_thresholds = parse_threshold_string(dist_threshold.clone());
            let angle_thresholds = parse_threshold_string(angle_threshold.clone());
            
            let loaded_index_vec = index_paths.into_par_iter().map(|index_path| {
                let (offset_path, value_path, lookup_path, hash_type_path) = get_offset_value_lookup_type(index_path);
                let config = read_index_config_from_file(&hash_type_path);
                let offset_table = read_offset_map(&offset_path, config.hash_type).expect(
                    &log_msg(FAIL, &format!("Failed to load offset table: {}", &offset_path))
                );
                let lookup = load_lookup_from_file(&lookup_path);
                (offset_table, lookup, config, value_path)
            }).collect::<Vec<()>>();
            // Iterate over queries
            queries.into_par_iter().for_each(|(pdb_path, query_string, output_path)| {
                let pdb_query = read_pdb(&pdb_path).expect(
                    &log_msg(FAIL, &format!("Failed to read pdb file: {}", &pdb_path))
                );
                let query_residues = parse_query_string(&query_string);
                // Get query map for each query in all indices
                let queried_from_indices = loaded_index_vec.into_par_iter().for_each(|(offset_table, lookup, config, value_path)| {
                    let hash_type = config.hash_type;
                    let num_bin_dist = config.num_bin_dist;
                    let num_bin_angle = config.num_bin_angle;
                    let mode = config.mode;
                    let pdb_query_map = make_query_map(
                        &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, &dist_thresholds, &angle_thresholds
                    );
                    let pdb_query = pdb_query_map.keys().cloned().collect::<Vec<_>>();
                    match mode {
                        IndexMode::Id => {
                            let (mmap, value_vec) = read_u16_vector(&value_path).expect(
                                &log_msg(FAIL, &format!("Failed to load value vector: {}", &value_path))
                            );
                            let query_count_map = count_query_idmode(
                                &pdb_query, &pdb_query_map, &offset_table, value_vec, &lookup
                            );
                            let match_count_filter = get_match_count_filter(
                                match_cutoff.clone(), pdb_query.len(), query_residues.len()
                            );
                            let query_count_vec: Vec<(usize, QueryResult)> = query_count_map.into_iter().filter(|(k, v)| {
                                v.total_match_count >= match_count_filter[0] && v.node_count >= match_count_filter[1] && 
                                v.edge_count >= match_count_filter[2] && v.exact_match_count >= match_count_filter[3] &&
                                v.overflow_count <= match_count_filter[4] && v.grid_count <= match_count_filter[5] &&
                                v.idf >= score_cutoff && v.nres <= num_res_cutoff && v.plddt >= plddt_cutoff 
                            });
                            query_count_vec
                        },
                        IndexMode::Grid => {
                            let (mmap, value_vec) = read_u8_vector(&value_path).expect(
                                &log_msg(FAIL, &format!("Failed to load value vector: {}", &value_path))
                            );
                            let query_count_map = count_query_gridmode(
                                &pdb_query, &pdb_query_map, &offset_table, value_vec, &lookup
                            );
                            let match_count_filter = get_match_count_filter(
                                match_cutoff.clone(), pdb_query.len(), query_residues.len()
                            );
                            let query_count_vec: Vec<(usize, QueryResult)> = query_count_map.into_iter().filter(|(k, v)| {
                                v.total_match_count >= match_count_filter[0] && v.node_count >= match_count_filter[1] && 
                                v.edge_count >= match_count_filter[2] && v.exact_match_count >= match_count_filter[3] &&
                                v.overflow_count <= match_count_filter[4] && v.grid_count <= match_count_filter[5] &&
                                v.idf >= score_cutoff && v.nres <= num_res_cutoff && v.plddt >= plddt_cutoff 
                            });
                            query_count_vec
                        },
                        IndexMode::Pos => {
                            todo!()
                        },
                    }; // match mode
                }).reduce(|| Vec::new(), |mut acc, x| { // merge queries from all indices
                    acc.extend(x);
                    acc
                }); // loaded_index_vec
                // Sort query_count_vec by idf
                measure_time!(queried_from_indices.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap()));
                // If output path is not empty, write to file
                if !output_path.is_empty() {
                    // Create file
                    let file = std::fs::File::create(&output_path).expect(
                        &log_msg(FAIL, &format!("Failed to create file: {}", &output_path))
                    );
                    let mut writer = std::io::BufWriter::new(file);
                    for (k, v) in query_count_vec.iter() {
                        writer.write_all(format!("{:?}\t{}\t{}\n", v, query_residues, index_path.clone().unwrap()
                    ).as_bytes()).expect(
                            &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
                        );
                    }
                } else {
                    for (k, v) in query_count_vec.iter() {
                        println!("{:?}\t{}\t{}", v, query_residues, index_path.clone().unwrap());
                    }
                }
            }); // queries
        }, // AppArgs::Query
        _ => {
            eprintln!("{}", HELP_QUERY);
            std::process::exit(1);
        }
    }
}
