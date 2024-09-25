                            IndexMode::Grid => {
                                let (value_mmap, value_vec) = measure_time!(read_u8_vector(&value_path).expect(
                                    &log_msg(FAIL, &format!("Failed to load value vector: {}", &value_path))
                                ));
                                let query_count_map = measure_time!(count_query_gridmode(
                                    &pdb_query, &pdb_query_map, &offset_table, value_vec, &lookup
                                ));
                                let mut match_count_filter = get_match_count_filter(
                                    match_cutoff.clone(), pdb_query.len(), query_residues.len()
                                );
                                if retrieve {
                                    match_count_filter[1] = if node_count > match_count_filter[1] {node_count} else {match_count_filter[1]};
                                    if verbose {
                                        print_log_msg(INFO, &format!("Filtering result with residue >= {}", match_count_filter[1]));
                                    }
                                }
                                let mut query_count_vec: Vec<(usize, QueryResult)> = query_count_map.into_par_iter().filter(|(_k, v)| {
                                    v.total_match_count >= match_count_filter[0] && v.node_count >= match_count_filter[1] && 
                                    v.edge_count >= match_count_filter[2] && v.exact_match_count >= match_count_filter[3] &&
                                    v.overflow_count <= match_count_filter[4] && v.grid_count <= match_count_filter[5] &&
                                    v.idf >= score_cutoff && v.nres <= num_res_cutoff && v.plddt >= plddt_cutoff 
                                }).collect();
                                if verbose {
                                    print_log_msg(INFO, &format!("Found {} structures from inverted index", query_count_vec.len()));
                                }
                                // IF retrieve is true, retrieve matching residues
                                if retrieve {
                                    measure_time!(query_count_vec.par_iter_mut().for_each(|(_, v)| {
                                        #[cfg(not(feature = "foldcomp"))]
                                        let retrieval_result = retrieval_wrapper(
                                            &v.id, node_count, &pdb_query,
                                            hash_type, num_bin_dist, num_bin_angle,
                                            &pdb_query_map, &query_structure, &query_indices,
                                            &aa_dist_map
                                        );
                                        #[cfg(feature = "foldcomp")]
                                        let retrieval_result = if using_foldcomp {
                                            retrieval_wrapper_for_foldcompdb(
                                                &v.id, node_count, &pdb_query,
                                                hash_type, num_bin_dist, num_bin_angle,
                                                &pdb_query_map, &query_structure, &query_indices,
                                                &aa_dist_map, &foldcomp_db_reader
                                            )
                                        } else {
                                            retrieval_wrapper(
                                                &v.id, node_count, &pdb_query,
                                                hash_type, num_bin_dist, num_bin_angle,
                                                &pdb_query_map, &query_structure, &query_indices,
                                                &aa_dist_map
                                            )
                                        };
                                        v.matching_residues = retrieval_result.0;
                                        v.matching_residues_processed = retrieval_result.1;
                                    }));
                                    query_count_vec.retain(|(_, v)| v.matching_residues.len() > 0);
                                    println!("{:?}", query_count_vec.len());
                                    drop(value_mmap);
                                    drop(match_count_filter);
                                    return query_count_vec;
                                }
                                drop(value_mmap);
                                drop(match_count_filter);
                                query_count_vec
                            },
                            IndexMode::Pos => {
                                todo!()
                            },