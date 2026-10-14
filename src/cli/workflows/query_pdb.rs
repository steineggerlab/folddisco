// File: query_pdb.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved
//! `folddisco query`: search motifs against an index and print hits or evidence rows.

use std::io::BufRead;
use std::io::Write;

use rayon::prelude::*;

use crate::cli::config::read_index_config_from_file;
use crate::controller::filter::{MatchFilter, StructureFilter};
use crate::controller::mode::QueryMode;
use crate::controller::ResidueMatch;
use crate::controller::sort::{MatchSortStrategy, StructureSortStrategy};
use crate::cli::*;
use crate::controller::io::{
    read_compact_structure, 
    check_and_get_indices,
    get_lookup_and_type,
    resolve_tid_path_from_index_prefix,
};
use crate::controller::expand::{IndexExpansion, ToleranceConfig};
use crate::controller::query::{index_matching_substitutions, make_query_map, parse_threshold_string, resolve_query_substitutions};
use crate::controller::substitution::SubstitutionScheme;
use crate::controller::count_query::count_query;
use crate::controller::result::{
    convert_structure_query_result_to_match_query_results, format_residue_matches,
    sort_and_print_match_query_result, sort_and_print_structure_query_result, StructureResult
};
use crate::controller::retrieve::retrieval_wrapper;
use crate::index::indextable::load_folddisco_index;
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::*;
use crate::structure::chain_id::{format_chain_residue, residue_list_needs_separator, ChainId};

#[cfg(feature = "foldcomp")]
use crate::controller::retrieve::retrieval_wrapper_for_foldcompdb;
#[cfg(feature = "foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;
#[cfg(feature = "foldcomp")]
use crate::structure::io::StructureFileFormat;
#[cfg(feature = "foldcomp")]
use crate::controller::io::get_foldcomp_db_path_with_prefix;

pub const HELP_QUERY: &str = "\
usage: folddisco query -p <i:PDB> -q <QUERY> -i <i:INDEX> [OPTIONS] 

input/output:
 -p, --pdb <PATH>                 Path of PDB file to query
 -q, --query <STR>                Query string that specifies residues or a text file containing query
 -i, --index <PATH>               Path of index table to load [REQUIRED]
 -o, --output <PATH>              Output file path [stdout]
 
search parameters:
 -t, --threads <INT>              Number of threads [1]
 -d, --distance <FLOAT>           Distance tolerance in Angstroms; widest value wins [0.5]
 -a, --angle <FLOAT>              Angle tolerance in degrees; widest value wins [5.0]
 --ca-distance <FLOAT>            C-alpha distance threshold in matching residues [1.0]
 --sampling-count <INT>           Number of sampled hashes to search [all]
 --sampling-ratio <FLOAT>         Sampling ratio for hashes used in searching. For long queries, smaller ratio is recommended [1.0]
 --freq-filter <FLOAT>            Skip hashes found in more than this fraction of structures [off]
 --length-penalty <FLOAT>         Length penalty for searching. Zero means no penalty and higher value gives more penalty to longer structures [0.5]
 --skip-match                     Skip matching residues
 --serial-index                   Handle residue indices serially

query expansion:
 --expand-radius <INT>            Feature dimensions of a residue pair allowed in a neighbouring bin at once [1]
 --sensitive                      Wider, slower search for deformed motifs: --expand-radius 2. Pair with --max-node <n_residues>
 --confident                      Keep only confident, full matches: at least 80% of the query residues
                                  matched (all of them for a 3-4 residue motif) within 1.0 A RMSD.
                                  Past 12 residues only the coverage is required; filters given explicitly
                                  win, and with --skip-match only the coverage applies
 --aa-subst <MODE>                Substitute every query residue without an explicit :ALT. Also sets the scheme for :*
                                  - blosum62: positive BLOSUM62 score (default for :*)
                                  - group: same class (RHK, DE, NQST, FWY, AVLIMC, GP)
                                  - size: same IMGT side-chain volume class (GAS, CDPNT, QEHV, MILKR, FWY)
                                  On an index built with expansion, what the index covers is off unless given.
                                  Hits match query-side expansion within the tolerance; bin-rounding margins differ

filtering options:
 --total-match <INT>              Filter out structures with less than total match count [0]
 --covered-node <INT>             Filter out structures not covered by given number of nodes with hashes [0]
 --covered-node-ratio <FLOAT>     Filter out structures not covered by given ratio of nodes with hashes [0.0]
 --max-node <INT>                 Filter out structures of maximum matching node size smaller than given value [0]
 --max-node-ratio <FLOAT>         Filter out structures of maximum matching node size smaller than given ratio [0.0]
 --score <FLOAT>                  IDF score cutoff [0.0]
 --connected-node <INT>           Filter out structures/matches with connected node count smaller than given value [0]
 --connected-node-ratio <FLOAT>   Filter out structures/matches with connected node count smaller than given ratio [0.0]
 --num-residue <INT>              Number of residues cutoff [50000]
 --plddt <FLOAT>                  pLDDT cutoff [0.0]
 --rmsd <FLOAT>                   Maximum RMSD cutoff [no limit]
 --tm-score <FLOAT>               Minimum TM-score cutoff [0.0]
 --gdt-ts <FLOAT>                 Minimum GDT-TS cutoff. Thresholds: 1.0Å, 2.0Å, 4.0Å, 8.0Å [0.0]
 --gdt-ha <FLOAT>                 Minimum GDT-HA cutoff. Thresholds: 0.5Å, 1.0Å, 2.0Å, 4.0Å [0.0]
 --chamfer <FLOAT>                Maximum Chamfer distance cutoff (mean nearest-neighbour distance) [no limit]
 --hausdorff <FLOAT>              Maximum Hausdorff distance cutoff (max nearest-neighbour distance) [no limit]
 --drmsd <FLOAT>                  Maximum dRMSD cutoff. Superposition-free, so tolerant of hinged motifs [no limit]
 --top <INT>                      Limit output to top N structures based on IDF score [all]

display options:
 --header                         Print header in output
 --web                            Print output for web
 --per-structure                  Print output per structure
 --per-match                      Print output per match. Not working with --skip-match
 --format-output <KEYS>           Comma-separated column names to output
                                  - Per-match: qid, tid, nid, db_key, node_count, idf, coverage_idf, match_score, rmsd,
                                    matching_residues, u_matrix, t_vector, matching_coordinates, query_residues,
                                    tm_score, gdt_ts, gdt_ha, chamfer_distance, hausdorff_distance, drmsd,
                                    max_dist_deviation
                                  - Per-structure: qid, tid, nid, db_key, total_match_count, node_count, edge_count, idf, nres,
                                    plddt, max_node_cov, min_rmsd, min_drmsd, structure_score, matching_residues, query_residues
                                  - Example: --format-output tid,idf,rmsd,tm_score
 --sort-by <KEYS>                 Comma-separated sort keys with optional :asc or :desc
                                  [per-match match_score:desc,rmsd:asc; per-structure structure_score:desc,min_rmsd:asc]
                                  - Per-match: match_score (idf x coverage^2 x tm_score), node_count, idf,
                                    coverage_idf (idf x matched fraction), rmsd, tm_score,
                                    gdt_ts, gdt_ha, chamfer_distance, hausdorff_distance,
                                    drmsd, max_dist_deviation
                                  - Per-structure: structure_score (matched^2 x sqrt(idf) / (1 + rmsd)), max_node_count,
                                    node_count, idf, min_rmsd, min_drmsd, total_match_count, edge_count, nres, plddt
                                  - Example: --sort-by tm_score,rmsd or --sort-by idf:desc
 --chain-sep                      Always write residues as CHAIN_RESIDUE (A_21). Otherwise `_` is used only
                                  for multi-character or numeric chains. -q accepts both
 --skip-ca-match                  Print matching residues before C-alpha distance check
 --partial-fit                    Superposition will find the best aligning substructure using LMS (Least Median of Squares)
 --superpose                      Print U, T, CA of matching residues

novelty screening:
 --novelty-mode                   One verdict row per query instead of the hit list. Columns: query_id,
                                  verdict, candidates, hits, index_coverage, best_hit, best_coverage,
                                  best_rmsd, best_residues, query_residues
                                  - verdict: KNOWN (best hit covers enough, close enough), PARTIAL
                                    (covered but under either threshold), NOVEL (nothing covered),
                                    NO_HASHES (query not searchable)
                                  - index_coverage: best --covered-node ratio before filters (hash level)
                                  - best_*: best hit after filters, its matched residues included;
                                    RMSD and residues are NA with --skip-match
                                  Filters such as a high --max-node remove partial matches; screen with low ones
 --novelty-coverage <FLOAT>       Residue coverage of the best hit needed to call a motif KNOWN [0.8]
 --novelty-rmsd <FLOAT>           Best-hit RMSD, in Angstroms, still allowed for KNOWN [2.0]. Covered but
                                  over the limit is PARTIAL. With --skip-match coverage alone decides

general options:
 -v, --verbose                    Print verbose messages
 -h, --help                       Print this help menu

examples:
# Search with default settings (ranked by match_score, then RMSD)
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6

# Print custom columns (tid, idf, RMSD, and TM-score only)
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --format-output tid,idf,rmsd,tm_score

# Print matches sorted by node count and TM-score
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --sort-by node_count,tm_score

# Print per-structure results sorted by max node count and IDF score
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --per-structure --sort-by max_node_count,idf

# Query a multi-character or numeric chain ID; CHAIN_RESIDUE is required here
folddisco query -p query/9A1O.cif -q 10_250,AA_312,AA_318 -i index/pdb_folddisco -t 6

# Query file given as separate text file
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 -d 0.5 -a 5

# Substitutions and ranges: residues 1-10, and residue 11 as any amino acid
folddisco query -p query/4CHA.pdb -q 1-10,11:X -i index/h_sapiens_folddisco -t 6 --serial-index

# BLOSUM62 substitutions for Asp102 only (:*), or for every residue (--aa-subst)
folddisco query -p query/4CHA.pdb -q B57,B102:*,C195 -i index/h_sapiens_folddisco -t 6
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --aa-subst group

# Filtering by connected node and RMSD
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --connected-node 0.75 --rmsd 1.0

# Coverage based filtering & top N filtering without residue matching
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --covered-node 3 --top 1000 --per-structure --skip-match

# Sensitive search for a deformed motif, ranked by dRMSD
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --sensitive \\
  --sort-by node_count,drmsd --format-output tid,node_count,idf,rmsd,drmsd,matching_residues

# Only confident, full matches (>= 80% of the query residues within 1 A)
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --confident

# Novelty evidence for designed motifs, one row per design
folddisco query -q designs.txt -i index/pdb_folddisco -t 6 --novelty-mode --header
";

pub const MIN_CONNECTED_COMPONENT_SIZE: usize = 2;
pub const MAX_NUM_LINES_FOR_WEB: usize = 1000;

/// Expansion radius `--sensitive` raises the search to; radius 3 measured no better
/// (see docs/feature_evaluation.md).
const SENSITIVE_EXPAND_RADIUS: usize = 2;
const DEFAULT_EXPAND_RADIUS: usize = 1;
const DEFAULT_DIST_THRESHOLD: &str = "0.5";
const DEFAULT_ANGLE_THRESHOLD: &str = "5.0";

/// `--confident`: fraction of query residues a match must cover. 0.8 keeps every residue of a
/// 3-4 residue motif and allows one missing residue from five residues up.
const CONFIDENT_COVERAGE: f32 = 0.8;
/// `--confident`: RMSD cap on the matched residues, in Angstroms.
const CONFIDENT_RMSD: f32 = 1.0;
/// Longest query `--confident` caps the RMSD of. Above it, matches are assembled from several
/// parts and run to several Angstroms while coverage alone is already precise
/// (docs/feature_evaluation.md §15).
const CONFIDENT_RMSD_MAX_RESIDUES: usize = 12;

/// Column names of a `--novelty-mode` row.
const NOVELTY_HEADER: &str = "query_id\tverdict\tcandidates\thits\tindex_coverage\tbest_hit\tbest_coverage\tbest_rmsd\tbest_residues\tquery_residues";
/// `--novelty-coverage` default: residue coverage of the best hit needed to call a motif KNOWN.
pub const NOVELTY_COVERAGE: f32 = 0.8;
/// `--novelty-rmsd` default, in Angstroms: best-hit RMSD still allowed for KNOWN.
pub const NOVELTY_RMSD: f32 = 2.0;

/// Run `folddisco query` for parsed `AppArgs::Query`.
pub fn query_pdb(env: AppArgs) {
    match env {
        AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            skip_match,
            dist_threshold,
            angle_threshold,
            ca_dist_threshold,
            expand_radius,
            sensitive,
            confident,
            aa_subst,
            total_match_count,
            covered_node_count,
            covered_node_ratio,
            max_matching_node_count,
            max_matching_node_ratio,
            idf_score_cutoff,
            connected_node_count,
            connected_node_ratio,
            num_res_cutoff,
            plddt_cutoff,
            rmsd_cutoff,
            tm_score_cutoff,
            gdt_ts_cutoff,
            gdt_ha_cutoff,
            chamfer_distance_cutoff,
            hausdorff_distance_cutoff,
            drmsd_cutoff,
            top_n,
            web_mode,
            sampling_count,
            sampling_ratio,
            freq_filter,
            length_penalty,
            sort_by,
            format_output,
            output_per_structure,
            output_per_match,
            output_with_superpose,
            skip_ca_match,
            partial_fit,
            header,
            serial_query,
            chain_separator,
            output,
            novelty_mode,
            novelty_coverage_threshold,
            novelty_rmsd_threshold,
            verbose,
            help: _,
        } => {
            if verbose { print_logo(); }
            // --help is handled in main.rs
            if index_path.is_none() {
                eprintln!("{}", HELP_QUERY);
                std::process::exit(1);
            }
            
            // Query mode decides which sort strategy applies
            let query_mode = QueryMode::from_flags(
                skip_match, web_mode, output_per_structure, output_per_match
            );
            
            if query_mode == QueryMode::ContradictoryPrintError {
                print_log_msg(FAIL, 
                    "Cannot print output per structure and per match at the same time. Use either --per-structure or --per-match"
                );
                std::process::exit(1);
            }
            
            
            let use_structure_sort = matches!(
                query_mode,
                QueryMode::PerStructure | QueryMode::SkipMatch
            );
            
            let match_sort_strategy = if !use_structure_sort {
                if sort_by.is_empty() {
                    MatchSortStrategy::default()
                } else {
                    MatchSortStrategy::from_str(&sort_by)
                        .unwrap_or_else(|e| {
                            print_log_msg(FAIL, &format!("Error parsing --sort-by: {}", e));
                            std::process::exit(1);
                        })
                }
            } else {
                MatchSortStrategy::default()
            };
            
            let structure_sort_strategy = if use_structure_sort {
                if sort_by.is_empty() {
                    StructureSortStrategy::default()
                } else {
                    StructureSortStrategy::from_str(&sort_by)
                        .unwrap_or_else(|e| {
                            print_log_msg(FAIL, &format!("Error parsing --sort-by: {}", e));
                            std::process::exit(1);
                        })
                }
            } else {
                StructureSortStrategy::default()
            };
            
            let parsed_columns: Option<Vec<String>> = format_output.map(|cols| {
                cols.split(',').map(|s| s.trim().to_string()).collect()
            });
            let column_refs: Option<Vec<&str>> = parsed_columns.as_ref().map(|cols| {
                cols.iter().map(|s| s.as_str()).collect()
            });

            if verbose  {
                if use_structure_sort {
                    print_log_msg(INFO, &format!("Printing results {} sorting with {}", query_mode, structure_sort_strategy));
                } else {
                    print_log_msg(INFO, &format!("Printing results {} sorting with {}", query_mode, match_sort_strategy));
                }
            }
            
            if verbose {
                if pdb_path.is_empty() {
                    print_log_msg(INFO, &format!("Querying {} to {}", &query_string, &index_path.clone().unwrap()));
                } else {
                    print_log_msg(INFO, &format!("Querying {}:{} to {}", &pdb_path, &query_string, &index_path.clone().unwrap()));
                }
            }
            
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            
            let index_paths = check_and_get_indices(index_path.clone(), verbose);
            if verbose {
                print_log_msg(INFO, &format!("Found {} index file(s). Querying with {} threads", index_paths.len(), threads));
            }
            
            let index_prefix = index_paths[0].clone();
            let (index, offset_mmap) = measure_time!(load_folddisco_index(&index_prefix), verbose);
            
            let (lookup_path, hash_type_path) = get_lookup_and_type(&index_prefix);
            let config = read_index_config_from_file(&hash_type_path);
            let lookup = measure_time!(load_lookup_from_file(&lookup_path), verbose);

            
            let queries = if query_string.ends_with(".txt") || query_string.ends_with(".tsv") {
                // TSV lines: structure path, query, output path
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
                vec![(pdb_path.clone(), query_string.clone(), output.clone())]
            };

            // Novelty rows are appended so queries sharing a file accumulate; empty each
            // distinct path once so a rerun replaces it.
            if novelty_mode {
                let mut truncated: Vec<&str> = Vec::new();
                let mut stdout_used = false;
                for (_, _, output_path) in queries.iter() {
                    if output_path.is_empty() {
                        stdout_used = true;
                        continue;
                    }
                    if truncated.contains(&output_path.as_str()) {
                        continue;
                    }
                    std::fs::File::create(output_path).expect(
                        &log_msg(FAIL, &format!("Failed to create file: {}", output_path))
                    );
                    if header {
                        print_novelty_row(NOVELTY_HEADER, output_path);
                    }
                    truncated.push(output_path);
                }
                if header && stdout_used {
                    println!("{}", NOVELTY_HEADER);
                }
            }

            let scheme = aa_subst.as_ref().map(|mode| SubstitutionScheme::from_str(mode).unwrap_or_else(|| {
                print_log_msg(FAIL, &format!("Unknown --aa-subst '{}'; use blosum62, group or size", mode));
                std::process::exit(1);
            }));

            let index_expansion = config.expansion.clone();
            let tolerance = query_tolerance(
                dist_threshold, angle_threshold, expand_radius, sensitive, index_expansion.as_ref(),
            );
            // Matching must accept whatever the expanded index returned
            let matching_tolerance = index_expansion.as_ref().map(|expansion| expansion.matching_tolerance(&tolerance));
            if verbose {
                print_log_msg(INFO, &format!(
                    "Tolerance: distance {:?} A, angle {:?} deg, expansion radius {}, substitution {}",
                    tolerance.dist_thresholds, tolerance.angle_thresholds,
                    tolerance.radius, scheme.map_or("none".to_string(), |s| s.to_string())
                ));
                if let Some(expansion) = &index_expansion {
                    print_log_msg(INFO, &format!(
                        "Index built with expansion (radius {}, substitution {}); lookup uses the tolerance above",
                        expansion.tolerance.radius, expansion.scheme.map_or("none".to_string(), |s| s.to_string())
                    ));
                }
            }
            if scheme.is_some() && config.hash_type.amino_acid_index().is_none() {
                print_log_msg(WARN, &format!(
                    "{:?} does not encode residue identity; --aa-subst has no effect", config.hash_type
                ));
            }

            // Load foldcomp db 
            #[cfg(feature = "foldcomp")]
            let using_foldcomp = config.foldcomp_db.is_some() && config.input_format == StructureFileFormat::FCZDB;

            #[cfg(feature = "foldcomp")]
            let foldcomp_db_reader = match config.input_format {
                StructureFileFormat::FCZDB => {
                    if !skip_match {
                        let mut foldcomp_db_path = config.foldcomp_db.clone().unwrap();
                        // Fall back to a Foldcomp DB next to the index
                        if !std::path::PathBuf::from(&foldcomp_db_path).is_file() {
                            let local_foldcomp_db_path = get_foldcomp_db_path_with_prefix(&index_prefix);
                            if local_foldcomp_db_path.is_some() {
                                foldcomp_db_path = local_foldcomp_db_path.unwrap();
                            }
                        }
                        measure_time!(FoldcompDbReader::new(foldcomp_db_path.as_str()), verbose)
                    } else {
                        FoldcompDbReader::empty()
                    }
                },
                _ => FoldcompDbReader::empty(),
            };

            // #[cfg(not(feature = "foldcomp"))]
            // let using_foldcomp = false;

            queries.into_par_iter().for_each(|(pdb_path, query_string, output_path)| {
                let (query_structure, _) = read_compact_structure(&pdb_path).expect(
                    &log_msg(FAIL, &format!("Failed to read structure: {}", &pdb_path))
                );
                
                let (query_residues, aa_substitutions) = parse_query_string(&query_string, query_structure.chains[0]);
                
                let residue_count = if query_residues.is_empty() {
                    query_structure.num_residues
                } else {
                    query_residues.len()
                };
                let query_string = if query_residues.is_empty() {
                    query_string
                } else {
                    let query_residues = query_residues.clone();
                    // query_residues.sort();
                    res_chain_to_string(&query_residues, chain_separator)
                };

                let (covered_node_ratio, max_matching_node_ratio, connected_node_ratio, rmsd_cutoff) =
                    confident_filters(
                        confident, skip_match, residue_count, covered_node_ratio,
                        max_matching_node_ratio, connected_node_ratio, rmsd_cutoff,
                    );

                let hash_type = config.hash_type;
                let num_bin_dist = config.num_bin_dist;
                let num_bin_angle = config.num_bin_angle;
                let dist_cutoff = config.grid_width;
                let multiple_bin = &config.multiple_bin;
                let total_structures = lookup.len() as f32;
                        
                let substitutions = resolve_query_substitutions(
                    &query_structure, &query_residues, &aa_substitutions,
                    scheme.unwrap_or(SubstitutionScheme::Blosum62), scheme.is_some(), serial_query,
                );
                let (pdb_query_map, query_indices, aa_dist_map ) = measure_time!(make_query_map(
                    &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, multiple_bin,
                    &tolerance, &substitutions, dist_cutoff, serial_query,
                    &Some(&index), total_structures
                ), verbose);
                let pdb_query = pdb_query_map.keys().cloned().collect::<Vec<_>>();

                // Against an expanded index, matching uses the combined expansion
                let matching_map = matching_tolerance.as_ref().map(|matching_tolerance| {
                    let substitutions = match index_expansion.as_ref().and_then(|e| e.scheme) {
                        Some(index_scheme) => index_matching_substitutions(
                            &query_structure, &query_residues, &substitutions, index_scheme, serial_query,
                        ),
                        None => substitutions.clone(),
                    };
                    let (map, _, aa_dist_map) = make_query_map(
                        &pdb_path, &query_residues, hash_type, num_bin_dist, num_bin_angle, multiple_bin,
                        matching_tolerance, &substitutions, dist_cutoff, serial_query,
                        &Some(&index), total_structures
                    );
                    let hashes = map.keys().cloned().collect::<Vec<_>>();
                    (map, hashes, aa_dist_map)
                });
                let (match_query_map, match_query, match_aa_dist_map) = match &matching_map {
                    Some((map, hashes, aa_dist_map)) => (map, hashes, aa_dist_map),
                    None => (&pdb_query_map, &pdb_query, &aa_dist_map),
                };
                if verbose {
                    print_log_msg(INFO, &format!(
                        "Expanded query into {} lookup and {} matching hashes", pdb_query.len(), match_query.len()
                    ));
                }
                // Residues farther apart than the index cutoff produce no hashes: not searchable
                if pdb_query.is_empty() {
                    let missing = query_residues.len().saturating_sub(query_indices.len());
                    print_log_msg(WARN, &format!(
                        "{}:{} produced no hashes: {}",
                        &pdb_path, &query_string,
                        if missing > 0 {
                            format!("{} of {} residues are not in the structure", missing, query_residues.len())
                        } else {
                            "no residue pair is within the index distance cutoff".to_string()
                        }
                    ));
                    if novelty_mode {
                        print_novelty_row(&novelty_no_hash_row(&pdb_path, &query_string), &output_path);
                        return;
                    }
                }
                let structure_filter = StructureFilter::new(
                    total_match_count, covered_node_count, covered_node_ratio,
                    idf_score_cutoff, num_res_cutoff, plddt_cutoff,
                    max_matching_node_count, max_matching_node_ratio, rmsd_cutoff,
                    drmsd_cutoff, residue_count,
                );

                let query_count_map = measure_time!(count_query(
                    &pdb_query, &pdb_query_map, &index, &lookup,
                    sampling_ratio, sampling_count, freq_filter, length_penalty
                ), verbose);
                // Evidence from the index before any filter ran
                let novelty_candidates = query_count_map.len();
                let novelty_index_coverage = if novelty_mode {
                    query_count_map.iter().map(|(_, v)| v.node_count).max().unwrap_or(0)
                } else {
                    0
                };
                let mut query_count_vec: Vec<(usize, StructureResult)> = query_count_map.into_par_iter().filter(|(_k, v)| {
                    structure_filter.filter_before_matching(v)
                }).collect();

                if verbose {
                    print_log_msg(INFO, &format!("Found {} structures from inverted index", query_count_vec.len()));
                }

                measure_time!(query_count_vec.par_sort_by(|a, b| b.1.idf.partial_cmp(&a.1.idf).unwrap()), verbose);
                if top_n != usize::MAX {
                    if verbose {
                        print_log_msg(INFO, &format!("Limiting result to top {} structures", top_n));
                    }
                    query_count_vec.truncate(top_n);
                }
                        
                // Residue matching unless --skip-match
                if !skip_match {
                    measure_time!(query_count_vec.par_iter_mut().for_each(|(_, v)| {
                        let joined_path;
                        let resolved_tid: &str = if std::path::Path::new(&v.tid).is_file() {
                            &v.tid
                        } else {
                            joined_path = resolve_tid_path_from_index_prefix(&v.tid, &index_prefix);
                            &joined_path
                        };
                        
                        #[cfg(not(feature = "foldcomp"))]
                        let retrieval_result = retrieval_wrapper(
                            resolved_tid, MIN_CONNECTED_COMPONENT_SIZE, match_query,
                            hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                            match_query_map, &query_structure, &query_indices,
                            match_aa_dist_map, ca_dist_threshold, partial_fit
                        );
                        #[cfg(feature = "foldcomp")]
                        let retrieval_result = if using_foldcomp {
                            retrieval_wrapper_for_foldcompdb(
                                v.db_key, MIN_CONNECTED_COMPONENT_SIZE, match_query,
                                hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                match_query_map, &query_structure, &query_indices,
                                match_aa_dist_map, ca_dist_threshold, partial_fit,
                                &foldcomp_db_reader
                            )
                        } else {
                            retrieval_wrapper(
                                resolved_tid, MIN_CONNECTED_COMPONENT_SIZE, match_query,
                                hash_type, num_bin_dist, num_bin_angle, multiple_bin, dist_cutoff,
                                match_query_map, &query_structure, &query_indices,
                                match_aa_dist_map, ca_dist_threshold, partial_fit,
                            )
                        };
                        v.matching_residues = retrieval_result.0;
                        v.matching_residues_processed = retrieval_result.1;
                        v.max_matching_node_count = retrieval_result.2;
                        v.min_rmsd_with_max_match = retrieval_result.3;
                        v.min_drmsd_with_max_match = retrieval_result.4;
                    }), verbose);

                    query_count_vec.retain(|(_, v)| structure_filter.filter_after_matching(v));
                }
                let mut queried_from_indices = query_count_vec;
                drop(query_residues);
                let evalue_cutoff = f64::MAX; // Currently not used
                let match_filter= MatchFilter::new(
                    connected_node_count, connected_node_ratio, idf_score_cutoff, evalue_cutoff,
                    rmsd_cutoff, tm_score_cutoff, gdt_ts_cutoff, gdt_ha_cutoff,
                    chamfer_distance_cutoff, hausdorff_distance_cutoff, drmsd_cutoff,
                    residue_count,
                );

                match query_mode {
                    QueryMode::PerMatch => {
                        let mut match_results = convert_structure_query_result_to_match_query_results(
                            &queried_from_indices, skip_ca_match, 
                            total_structures as usize, residue_count
                        );
                        match_results.retain(|(_, v)| match_filter.filter(v));
                        if novelty_mode {
                            // Highest residue coverage wins, lower RMSD breaks ties
                            let best = match_results.iter().map(|(_, v)| v).max_by(
                                |a, b| a.node_count.cmp(&b.node_count).then_with(
                                    || b.rmsd.partial_cmp(&a.rmsd).unwrap_or(std::cmp::Ordering::Equal)
                                )
                            ).map(|v| (v.tid, v.node_count, Some(v.rmsd),
                                       format_residue_matches(&v.matching_residues, chain_separator)));
                            let hits = match_results.iter().filter(|(_, v)| v.node_count > 0)
                                .map(|(_, v)| v.tid).collect::<std::collections::HashSet<_>>().len();
                            novelty_evidence(
                                &pdb_path, &query_string, residue_count, novelty_candidates, hits,
                                novelty_index_coverage, novelty_coverage_threshold,
                                novelty_rmsd_threshold, best, &output_path,
                            );
                        } else {
                            sort_and_print_match_query_result(
                                &mut match_results, top_n, 
                                &output_path, &pdb_path, &query_string, 
                                column_refs.as_deref(), output_with_superpose, header, verbose,
                                match_sort_strategy.clone(), chain_separator,
                            );
                        }
                    }
                    QueryMode::Web => {
                        let mut match_results = convert_structure_query_result_to_match_query_results(
                            &queried_from_indices, skip_ca_match, 
                            total_structures as usize, residue_count
                        );
                        match_results.retain(|(_, v)| match_filter.filter(v));
                        // Web output always includes superposition
                        sort_and_print_match_query_result(
                            &mut match_results, MAX_NUM_LINES_FOR_WEB,
                            &output_path, &pdb_path, &query_string, 
                            column_refs.as_deref(), true, header, verbose,
                            match_sort_strategy.clone(), chain_separator,
                        );
                    }
                    QueryMode::PerStructure | QueryMode::SkipMatch => {
                        if novelty_mode {
                            // Without matching, IDF breaks coverage ties instead of RMSD
                            let best = queried_from_indices.iter().map(|(_, v)| v).max_by(|a, b| if skip_match {
                                a.node_count.cmp(&b.node_count).then_with(
                                    || a.idf.partial_cmp(&b.idf).unwrap_or(std::cmp::Ordering::Equal)
                                )
                            } else {
                                a.max_matching_node_count.cmp(&b.max_matching_node_count).then_with(
                                    || b.min_rmsd_with_max_match.partial_cmp(&a.min_rmsd_with_max_match)
                                        .unwrap_or(std::cmp::Ordering::Equal)
                                )
                            }).map(|v| if skip_match {
                                (v.tid, v.node_count, None, String::new())
                            } else {
                                (v.tid, v.max_matching_node_count, Some(v.min_rmsd_with_max_match),
                                 best_component_residues(v, chain_separator))
                            });
                            let hits = queried_from_indices.iter().filter(|(_, v)| {
                                if skip_match { v.node_count > 0 } else { v.max_matching_node_count > 0 }
                            }).count();
                            novelty_evidence(
                                &pdb_path, &query_string, residue_count, novelty_candidates,
                                hits, novelty_index_coverage, novelty_coverage_threshold,
                                novelty_rmsd_threshold, best, &output_path,
                            );
                        } else {
                            sort_and_print_structure_query_result(
                                &mut queried_from_indices, &output_path, 
                                &pdb_path, &query_string, column_refs.as_deref(), header, verbose, structure_sort_strategy.clone(), chain_separator,
                            );
                        }
                    }
                    QueryMode::ContradictoryPrintError => {
                        print_log_msg(FAIL, "Invalid query mode");
                        std::process::exit(1);
                    }
                }
                drop(queried_from_indices);
                drop(query_structure);
            }); // queries
            drop(lookup);
            drop(offset_mmap);
            drop(index);
        }, // AppArgs::Query
        _ => {
            eprintln!("{}", HELP_QUERY);
            std::process::exit(1);
        }
    }
}

/// Query-side tolerance. Unset options take their defaults, except on an index that
/// already stores geometric expansion, where they stay off unless any is given.
/// `--sensitive` only raises the radius.
fn query_tolerance(
    dist_threshold: Option<String>, angle_threshold: Option<String>, expand_radius: Option<usize>,
    sensitive: bool, index_expansion: Option<&IndexExpansion>,
) -> ToleranceConfig {
    let index_covers_geometry = index_expansion.map_or(false, |e| e.tolerance.radius > 0);
    let geometry_given = dist_threshold.is_some() || angle_threshold.is_some() || expand_radius.is_some() || sensitive;
    let use_defaults = !index_covers_geometry || geometry_given;
    let default = |value: Option<String>, fallback: &str| value.or(use_defaults.then(|| fallback.to_string()));
    let radius = expand_radius.unwrap_or(if use_defaults { DEFAULT_EXPAND_RADIUS } else { 0 });
    ToleranceConfig::new(
        parse_threshold_string(default(dist_threshold, DEFAULT_DIST_THRESHOLD)),
        parse_threshold_string(default(angle_threshold, DEFAULT_ANGLE_THRESHOLD)),
        if sensitive { radius.max(SENSITIVE_EXPAND_RADIUS) } else { radius },
    )
}

/// Coverage and RMSD filters `--confident` sets for a query of `residue_count` residues,
/// leaving any filter given explicitly alone. With `--skip-match` only hash coverage can be
/// checked, and above `CONFIDENT_RMSD_MAX_RESIDUES` residues only coverage is required.
fn confident_filters(
    confident: bool, skip_match: bool, residue_count: usize, covered_node_ratio: f32,
    max_matching_node_ratio: f32, connected_node_ratio: f32, rmsd_cutoff: f32,
) -> (f32, f32, f32, f32) {
    if !confident {
        return (covered_node_ratio, max_matching_node_ratio, connected_node_ratio, rmsd_cutoff);
    }
    let or_default = |given: f32, preset: f32| if given > 0.0 { given } else { preset };
    let rmsd_preset = if residue_count <= CONFIDENT_RMSD_MAX_RESIDUES { CONFIDENT_RMSD } else { 0.0 };
    (
        if skip_match { or_default(covered_node_ratio, CONFIDENT_COVERAGE) } else { covered_node_ratio },
        or_default(max_matching_node_ratio, CONFIDENT_COVERAGE),
        or_default(connected_node_ratio, CONFIDENT_COVERAGE),
        or_default(rmsd_cutoff, rmsd_preset),
    )
}

/// Fraction of the query's residues that `nodes` covers; 0 for an empty query.
fn residue_coverage(nodes: usize, query_residue_count: usize) -> f32 {
    if query_residue_count > 0 {
        nodes as f32 / query_residue_count as f32
    } else {
        0.0
    }
}

/// Residues of the best connected component of a per-structure hit: most matched residues
/// first, lower RMSD breaking ties - the component `max_matching_node_count` counts.
fn best_component_residues(result: &StructureResult, chain_separator: bool) -> String {
    let components = if result.matching_residues_processed.is_empty() {
        &result.matching_residues
    } else {
        &result.matching_residues_processed
    };
    components.iter()
        .max_by(|a, b| {
            let matched = |r: &Vec<ResidueMatch>| r.iter().filter(|x| x.is_some()).count();
            matched(&a.0).cmp(&matched(&b.0))
                .then_with(|| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal))
        })
        .map_or(String::new(), |(residues, _, _, _, _, _, _)| {
            format_residue_matches(residues, chain_separator)
        })
}

/// Novelty row for a query that produced no hashes and so could not be searched. Kept
/// distinct from NOVEL: an unrepresentable query is not a discovery.
fn novelty_no_hash_row(query_id: &str, query_residues: &str) -> String {
    format!("{}\tNO_HASHES\t0\t0\t0.0000\tNA\tNA\tNA\tNA\t{}", query_id, query_residues)
}

/// Print a novelty row, warning when candidates existed but none survived.
fn novelty_evidence(
    query_id: &str, query_residues: &str, query_residue_count: usize,
    candidates: usize, hits: usize, index_coverage: usize,
    coverage_threshold: f32, rmsd_threshold: f32,
    best: Option<(&str, usize, Option<f32>, String)>, output_path: &str,
) {
    let row = novelty_verdict_row(
        query_id, query_residues, query_residue_count, candidates, hits, index_coverage,
        coverage_threshold, rmsd_threshold, best,
    );
    if row.split('\t').nth(1) == Some("NOVEL") && index_coverage > 0 {
        print_log_msg(WARN, &format!(
            "{}:{} - {} candidates covered up to {} of {} residues, but filters and matching kept none",
            query_id, query_residues, candidates, index_coverage, query_residue_count
        ));
    }
    print_novelty_row(&row, output_path);
}

/// One verdict row per query, with the evidence behind it. Columns follow `NOVELTY_HEADER`.
///
/// KNOWN needs both enough coverage and a close enough best hit; anything covered but
/// failing either is PARTIAL; no covered hit at all is NOVEL. Coverage alone does not make
/// a motif known - the same residues in another arrangement cover everything and are still
/// a different motif - so the RMSD the row prints gates the verdict too. With `--skip-match`
/// there is no RMSD and coverage is all there is.
///
/// `best` is the highest-coverage hit after filters as (target, covered residues, RMSD,
/// matched residues); RMSD and residues are absent when matching was skipped.
fn novelty_verdict_row(
    query_id: &str, query_residues: &str, query_residue_count: usize,
    candidates: usize, hits: usize, index_coverage: usize,
    coverage_threshold: f32, rmsd_threshold: f32,
    best: Option<(&str, usize, Option<f32>, String)>,
) -> String {
    let index_coverage_ratio = residue_coverage(index_coverage, query_residue_count);
    let (verdict, best_hit, best_coverage, best_rmsd, best_residues) = match best {
        Some((tid, covered, rmsd, residues)) if covered > 0 => {
            let coverage = residue_coverage(covered, query_residue_count);
            let close_enough = rmsd.map_or(true, |rmsd| rmsd <= rmsd_threshold);
            let verdict = if coverage >= coverage_threshold && close_enough { "KNOWN" } else { "PARTIAL" };
            (verdict, tid.to_string(), format!("{:.4}", coverage),
             rmsd.map_or("NA".to_string(), |rmsd| format!("{:.4}", rmsd)),
             if residues.is_empty() { "NA".to_string() } else { residues })
        }
        _ => ("NOVEL", "NA".to_string(), "0.0000".to_string(), "NA".to_string(), "NA".to_string()),
    };
    format!(
        "{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}",
        query_id, verdict, candidates, hits, index_coverage_ratio, best_hit, best_coverage,
        best_rmsd, best_residues, query_residues
    )
}

/// Append one row to `output_path`, or print it to stdout. One `write_all` per line
/// keeps parallel queries from interleaving on a local filesystem.
fn print_novelty_row(line: &str, output_path: &str) {
    if output_path.is_empty() {
        println!("{}", line);
        return;
    }
    let mut file = std::fs::OpenOptions::new().create(true).append(true).open(output_path).expect(
        &log_msg(FAIL, &format!("Failed to open file: {}", output_path))
    );
    file.write_all(format!("{}\n", line).as_bytes()).expect(
        &log_msg(FAIL, &format!("Failed to write to file: {}", output_path))
    );
}

/// Render the `query_residues` column: `A21,A23`, or `A_21,...` when a chain needs a
/// separator or `--chain-sep` is set.
pub fn res_chain_to_string(res_chain: &Vec<(ChainId, u64)>, chain_sep: bool) -> String {
    let separator = chain_sep
        || residue_list_needs_separator(res_chain.iter().map(|(chain, _)| chain));
    let mut output = String::new();
    for (i, (chain, res)) in res_chain.iter().enumerate() {
        output.push_str(&format_chain_residue(chain, *res, separator));
        if i < res_chain.len() - 1 {
            output.push(',');
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    #[ignore]
    fn test_query_pdb_workflow() {
        let pdb_path = String::from("data/serine_peptidase/4cha.pdb");
        let query_string = String::from("B57,B102,C195");
        let threads = 1;
        let index_path = Some(String::from("data/serine_peptidases_pdbtr_small"));
        let env = AppArgs::Query {
            pdb_path,
            query_string,
            threads,
            index_path,
            skip_match: false,
            dist_threshold: None,
            angle_threshold: None,
            ca_dist_threshold: 1.0,
            expand_radius: None,
            sensitive: false,
            confident: false,
            aa_subst: None,
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            idf_score_cutoff: 0.0,
            connected_node_count: 0,
            connected_node_ratio: 0.0,
            num_res_cutoff: 3000,
            plddt_cutoff: 0.0,
            rmsd_cutoff: 1.0,
            tm_score_cutoff: 0.0,
            gdt_ts_cutoff: 0.0,
            gdt_ha_cutoff: 0.0,
            chamfer_distance_cutoff: 0.0,
            hausdorff_distance_cutoff: 0.0,
            drmsd_cutoff: 0.0,
            top_n: 1000,
            web_mode: false,
            sampling_count: None,
            sampling_ratio: None,
            freq_filter: None,
            length_penalty: None,
            sort_by: String::from("node_count,rmsd"),
            format_output: None,
            output_per_structure: false,
            output_per_match: true,
            output_with_superpose: false,
            skip_ca_match: false,
            partial_fit: false,
            header: true,
            serial_query: false,
            chain_separator: false,
            output: String::from(""),
            novelty_mode: false,
            novelty_coverage_threshold: NOVELTY_COVERAGE,
            novelty_rmsd_threshold: NOVELTY_RMSD,
            verbose: true,
            help: false,
        };
        query_pdb(env);
    }
    #[test]
    #[ignore]
    fn test_query_with_foldcompdb() {
        #[cfg(feature = "foldcomp")] {
            let pdb_path = String::from("data/foldcomp/example_db:d1asha_");
            let query_string = String::from("1,2,3,4");
            let threads = 1;
            let index_path = Some(String::from("data/example_db_folddisco_db"));
            let env = AppArgs::Query {
                pdb_path,
                query_string,
                threads,
                index_path,
                skip_match: false,
                dist_threshold: None,
                angle_threshold: None,
                ca_dist_threshold: 1.0,
                expand_radius: None,
                sensitive: false,
                confident: false,
                aa_subst: None,
                total_match_count: 0,
                covered_node_count: 0,
                covered_node_ratio: 0.0,
                idf_score_cutoff: 0.0,
                connected_node_count: 0,
                connected_node_ratio: 0.0,
                max_matching_node_count: 0,
                max_matching_node_ratio: 0.0,
                num_res_cutoff: 3000,
                plddt_cutoff: 0.0,
                rmsd_cutoff: 1.0,
                tm_score_cutoff: 0.0,
                gdt_ts_cutoff: 0.0,
                gdt_ha_cutoff: 0.0,
                chamfer_distance_cutoff: 0.0,
                hausdorff_distance_cutoff: 0.0,
                drmsd_cutoff: 0.0,
                top_n: 1000,
                web_mode: false,
                sampling_count: None,
                sampling_ratio: None,
                freq_filter: None,
                length_penalty: None,
                sort_by: String::from("node_count,rmsd"),
                format_output: None,
                output_per_structure: true,
                output_per_match: false,
                output_with_superpose: true,
                skip_ca_match: false,
                header: true,
                serial_query: false,
                chain_separator: false,
                output: String::from(""),
                novelty_mode: false,
                novelty_coverage_threshold: NOVELTY_COVERAGE,
                novelty_rmsd_threshold: NOVELTY_RMSD,
                verbose: true,
                partial_fit: false,
                help: false,
            };
            query_pdb(env);
        }
    }
    #[test]
    #[ignore]
    fn test_query_pdb_with_file() {
        let pdb_path = String::from("");
        let query_string = String::from("data/query.tsv");
        let threads = 4;
        let index_path = Some(String::from("analysis/e_coli/test"));
        let env = AppArgs::Query {
            pdb_path,
            query_string,
            threads, 
            index_path,
            skip_match: true,
            dist_threshold: None,
            angle_threshold: None,
            ca_dist_threshold: 1.0,
            expand_radius: None,
            sensitive: false,
            confident: false,
            aa_subst: None,
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            idf_score_cutoff: 0.0,
            connected_node_count: 0,
            connected_node_ratio: 0.0,
            num_res_cutoff: 3000,
            plddt_cutoff: 0.0,
            rmsd_cutoff: 1.0,
            tm_score_cutoff: 0.0,
            gdt_ts_cutoff: 0.0,
            gdt_ha_cutoff: 0.0,
            chamfer_distance_cutoff: 0.0,
            hausdorff_distance_cutoff: 0.0,
            drmsd_cutoff: 0.0,
            top_n: 1000,
            web_mode: false,
            sampling_count: None,
            sampling_ratio: None,
            freq_filter: None,
            length_penalty: None,
            sort_by: String::from("node_count,rmsd"),
            format_output: None,
            output_per_structure: true,
            output_per_match: false,
            output_with_superpose: true,
            skip_ca_match: false,
            partial_fit: false,
            header: true,
            serial_query: false,
            chain_separator: false,
            output: String::from(""),
            novelty_mode: false,
            novelty_coverage_threshold: NOVELTY_COVERAGE,
            novelty_rmsd_threshold: NOVELTY_RMSD,
            verbose: true,
            help: false,
        };
        query_pdb(env);
    }

    #[test]
    fn confident_preset_leaves_explicit_filters_alone() {
        // Off: every value passes through
        assert_eq!(confident_filters(false, false, 3, 0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0));
        // On: coverage on both filter stages, RMSD cap, no pre-match coverage while matching
        assert_eq!(
            confident_filters(true, false, 3, 0.0, 0.0, 0.0, 0.0),
            (0.0, CONFIDENT_COVERAGE, CONFIDENT_COVERAGE, CONFIDENT_RMSD)
        );
        // --skip-match: hash coverage is all there is
        assert_eq!(
            confident_filters(true, true, 3, 0.0, 0.0, 0.0, 0.0),
            (CONFIDENT_COVERAGE, CONFIDENT_COVERAGE, CONFIDENT_COVERAGE, CONFIDENT_RMSD)
        );
        // A long query keeps the coverage rule and drops the RMSD cap
        assert_eq!(
            confident_filters(true, false, CONFIDENT_RMSD_MAX_RESIDUES + 1, 0.0, 0.0, 0.0, 0.0),
            (0.0, CONFIDENT_COVERAGE, CONFIDENT_COVERAGE, 0.0)
        );
        // Explicit values win
        assert_eq!(
            confident_filters(true, true, 23, 0.5, 1.0, 0.9, 2.0),
            (0.5, 1.0, 0.9, 2.0)
        );
    }

    #[test]
    fn query_tolerance_defaults_depend_on_index_expansion() {
        let some = |x: &str| Some(x.to_string());
        let default = ToleranceConfig::new(vec![0.5], vec![5.0], 1);
        // Plain index: identical to the previous defaults
        assert_eq!(query_tolerance(None, None, None, false, None), default);
        assert_eq!(query_tolerance(None, None, None, true, None).radius, 2);
        assert_eq!(query_tolerance(some("1.0"), None, Some(0), true, None), ToleranceConfig::new(vec![1.0], vec![5.0], 2));

        let geometric = IndexExpansion::new(1, 0.5, 5.0, None).unwrap();
        let residues_only = IndexExpansion::new(0, 0.5, 5.0, Some(SubstitutionScheme::Blosum62)).unwrap();
        // Index covers geometry: nothing extra unless asked for
        assert_eq!(query_tolerance(None, None, None, false, Some(&geometric)), ToleranceConfig::new(vec![], vec![], 0));
        // Any explicit option brings the rest of the defaults with it
        assert_eq!(query_tolerance(None, None, None, true, Some(&geometric)), ToleranceConfig::new(vec![0.5], vec![5.0], 2));
        assert_eq!(query_tolerance(some("1.0"), None, None, false, Some(&geometric)), ToleranceConfig::new(vec![1.0], vec![5.0], 1));
        // Substitution-only index keeps the geometric defaults
        assert_eq!(query_tolerance(None, None, None, false, Some(&residues_only)), default);
    }

    #[test]
    fn novelty_verdict_needs_both_coverage_and_a_close_hit() {
        let residues = "A10,A20,A30";
        let row = |best| novelty_verdict_row("d.pdb", residues, 3, 12, 4, 3, 0.8, 2.0, best);
        let hit = |covered, rmsd, res: &str| Some(("1abc", covered, rmsd, res.to_string()));
        // Full coverage within the RMSD limit is KNOWN, and the matched residues come with it
        assert_eq!(row(hit(3, Some(1.5), "B10,B20,B30")),
            "d.pdb\tKNOWN\t12\t4\t1.0000\t1abc\t1.0000\t1.5000\tB10,B20,B30\tA10,A20,A30");
        // Same residues in another arrangement: covered, but too far to be the same motif
        assert_eq!(row(hit(3, Some(4.2), "B10,B20,B30")),
            "d.pdb\tPARTIAL\t12\t4\t1.0000\t1abc\t1.0000\t4.2000\tB10,B20,B30\tA10,A20,A30");
        // Close, but not enough of the motif
        assert_eq!(row(hit(2, Some(0.02), "B10,_,B30")),
            "d.pdb\tPARTIAL\t12\t4\t1.0000\t1abc\t0.6667\t0.0200\tB10,_,B30\tA10,A20,A30");
        // --skip-match: no RMSD and no residues, so coverage alone decides
        assert_eq!(novelty_verdict_row("d.pdb", residues, 3, 5, 5, 2, 0.8, 2.0,
                                       Some(("1abc", 3, None, String::new()))),
            "d.pdb\tKNOWN\t5\t5\t0.6667\t1abc\t1.0000\tNA\tNA\tA10,A20,A30");
    }

    #[test]
    fn nothing_covered_is_novel_whether_or_not_candidates_existed() {
        let residues = "A57,A102,A195";
        // Candidates existed (2 of 3 residues) but filters kept none: still NOVEL, and the
        // caller warns about the filters
        let expected = "c.pdb\tNOVEL\t7\t0\t0.6667\tNA\t0.0000\tNA\tNA\tA57,A102,A195";
        assert_eq!(novelty_verdict_row("c.pdb", residues, 3, 7, 0, 2, 0.8, 2.0, None), expected);
        assert_eq!(novelty_verdict_row("c.pdb", residues, 3, 7, 1, 2, 0.8, 2.0,
                                       Some(("1ab9", 0, Some(0.0), String::new()))),
            "c.pdb\tNOVEL\t7\t1\t0.6667\tNA\t0.0000\tNA\tNA\tA57,A102,A195");
        // An empty query never divides by zero, and is never KNOWN
        assert_eq!(novelty_verdict_row("c.pdb", "", 0, 1, 1, 1, 0.8, 2.0,
                                       Some(("1abc", 0, None, String::new()))),
            "c.pdb\tNOVEL\t1\t1\t0.0000\tNA\t0.0000\tNA\tNA\t");
    }

    #[test]
    fn no_hash_row_and_header_have_the_same_columns() {
        let row = novelty_no_hash_row("query/1SU6.pdb", "A4,A39,A70");
        assert_eq!(row, "query/1SU6.pdb\tNO_HASHES\t0\t0\t0.0000\tNA\tNA\tNA\tNA\tA4,A39,A70");
        let columns = NOVELTY_HEADER.split('\t').count();
        assert_eq!(row.split('\t').count(), columns);
        assert_eq!(novelty_verdict_row("q", "A1", 1, 0, 0, 0, 0.8, 2.0, None).split('\t').count(), columns);
    }

}
