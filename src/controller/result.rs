// QueryResult struct and its implementation

use std::fmt;
use std::io::Write;
use rayon::slice::ParallelSliceMut;

use crate::measure_time;
use crate::prelude::{log_msg, print_log_msg, FAIL, INFO};
use crate::structure::coordinate::Coordinate;

use super::ResidueMatch;


pub const STRUCTURE_QUERY_RESULT_HEADER: &str = "id\tidf_score\ttotal_match_count\tnode_count\tedge_count\tmax_node_cov\tmin_rmsd\tmax_tm_score\tnres\tplddt\tmatching_residues\tkey\tquery_residues";
pub const MATCH_QUERY_RESULT_HEADER: &str = "id\tnode_count\tidf_score\trmsd\tmatching_residues\tkey\tquery_residues";
pub const MATCH_QUERY_RESULT_SUPERPOSE_HEADER: &str = "id\tnode_count\tidf_score\trmsd\tmatching_residues\tu_vector\tt_vector\ttarget_ca_coords\tkey\tquery_residues";

pub struct StructureResult<'a> {
    pub id: &'a str,
    pub nid: usize,
    pub db_key: usize, // Database key for the structure
    pub total_match_count: usize,
    pub node_count: usize,
    pub edge_count: usize,  // Add edge count field
    pub idf: f32,
    pub nres: usize,
    pub plddt: f32,
    pub matching_residues: Vec<(Vec<ResidueMatch>, f32, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>)>, // Match with connected components (added TM-score)
    pub matching_residues_processed: Vec<(Vec<ResidueMatch>, f32, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>)>, // Match with c-alpha distances (added TM-score)
    pub max_matching_node_count: usize,
    pub min_rmsd_with_max_match: f32,
    pub max_tm_score_with_max_match: f32,
}

impl<'a> StructureResult<'a> {
    pub fn new(
        id: &'a str, nid: usize, total_match_count: usize, node_count: usize, edge_count: usize,
        idf: f32, nres: usize, plddt: f32, db_key: usize
    ) -> Self {
        Self {
            id,
            nid,
            db_key,
            total_match_count,
            node_count,
            edge_count: edge_count,
            idf,
            nres,
            plddt,
            matching_residues: Vec::new(),
            matching_residues_processed: Vec::new(),
            max_matching_node_count: 0,
            min_rmsd_with_max_match: 0.0,
            max_tm_score_with_max_match: 0.0,
        }
    }

    pub fn into_match_query_results(&self, skip_ca_dist: bool) -> Vec<MatchResult> {
        match skip_ca_dist {
            false => self.matching_residues_processed.iter().enumerate().map(|(i, (residues, rmsd, tm_score, u_matrix, t_matrix, matching_coordinates))| {
                // WARNING: NOTE: Thinking of getting the match specific idf by saving the idf for each edge
                MatchResult::new(
                    self.id, i, self.idf, residues.clone(), *rmsd, *tm_score,
                    *u_matrix, *t_matrix, matching_coordinates.clone(), self.db_key
                )
            }).collect(),
            true => self.matching_residues.iter().enumerate().map(|(i, (residues, rmsd, tm_score, u_matrix, t_matrix, matching_coordinates))| {
                // WARNING: NOTE: Thinking of getting the match specific idf by saving the idf for each edge
                MatchResult::new(
                    self.id, i, self.idf, residues.clone(), *rmsd, *tm_score,
                    *u_matrix, *t_matrix, matching_coordinates.clone(), self.db_key
                )
            }).collect(),
        }
    }
}

pub fn convert_structure_query_result_to_match_query_results<'a>(
    results: &'a [(usize, StructureResult<'a>)], skip_ca_dist: bool
) -> Vec<(usize, MatchResult<'a>)> {
    results
        .iter()
        .flat_map(|(k, v)| {
            v.into_match_query_results(skip_ca_dist)
             .into_iter()
             .map(|x| (*k, x))
        })
        .collect()
}

impl<'a> fmt::Display for StructureResult<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let matching_residues_processed_with_score = if self.matching_residues_processed.len() == 0 {
            "NA".to_string()
        } else {
            self.matching_residues_processed.iter().map(
                // Only print score with 4 decimal places
                // Join with comma
                |(x, y, z, _, _, _)| format!("{}:{:.4}:{:.4}", x.iter().map(|x| {
                    match x {
                        // Convert u8 to char
                        Some((a, b)) => format!("{}{}", *a as char, b),
                        None => "_".to_string()
                    }
                }).collect::<Vec<String>>().join(","), y, z)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{:.4}\t{}\t{}", 
            self.id, self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.max_matching_node_count, self.min_rmsd_with_max_match, self.max_tm_score_with_max_match,
            self.nres, self.plddt, matching_residues_processed_with_score, self.db_key
        )
    }
}

impl<'a> fmt::Debug for StructureResult<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let matching_residues_processed_with_score = if self.matching_residues_processed.len() == 0 {
            "NA".to_string()
        } else {
            self.matching_residues_processed.iter().map(
                // Only print score with 4 decimal places
                |(x, y, z, _, _, _)| format!("{}:{:.4}:{:.4}", x.iter().map(|x| {
                    match x {
                        Some((a, b)) => format!("{}{}", *a as char, b),
                        None => "_".to_string()
                    }
                }).collect::<Vec<String>>().join(","), y, z)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{:.4}\t{}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.max_matching_node_count, self.min_rmsd_with_max_match, self.max_tm_score_with_max_match,
            self.nres, self.plddt, matching_residues_processed_with_score, self.db_key
        )
    }
}
// write_fmt
impl<'a> StructureResult<'a> {
    pub fn write_fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let matching_residues_processed_with_score = if self.matching_residues_processed.len() == 0 {
            "NA".to_string()
        } else {
            self.matching_residues_processed.iter().map(
                // Only print score with 4 decimal places
                |(x, y, z, _, _, _)| format!("{}:{:.4}:{:.4}", x.iter().map(|x| {
                    match x {
                        Some((a, b)) => format!("{}{}", *a as char, b),
                        None => "_".to_string()
                    }
                }).collect::<Vec<String>>().join(","), y, z)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{:.4}\t{}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.max_matching_node_count, self.min_rmsd_with_max_match, self.max_tm_score_with_max_match,
            self.nres, self.plddt, matching_residues_processed_with_score, self.db_key
        )
    }
}

pub struct MatchResult<'a> {
    pub id: &'a str,
    pub nid: usize,
    pub db_key: usize, // Database key for the structure
    pub node_count: usize,
    pub idf: f32,
    pub matching_residues: Vec<ResidueMatch>,
    pub rmsd: f32,
    pub tm_score: f32,
    pub u_matrix: [[f32; 3]; 3],
    pub t_matrix: [f32; 3],
    pub matching_coordinates: Vec<Coordinate>,
}

impl<'a> MatchResult<'a> {
    pub fn new(
        id: &'a str, nid: usize, avg_idf: f32, matching_residues: Vec<ResidueMatch>, rmsd: f32, tm_score: f32,
        u_matrix: [[f32; 3]; 3], t_matrix: [f32; 3], matching_coordinates: Vec<Coordinate>, db_key: usize
    ) -> Self {
        //
        let node_count = matching_residues.iter().map(|x| {
            match x {
                Some(_) => 1,
                None => 0
            }
        }).sum();
        Self {
            id,
            nid,
            db_key,
            node_count,
            idf: avg_idf,
            matching_residues,
            rmsd,
            tm_score,
            u_matrix,
            t_matrix,
            matching_coordinates,
        }
    }
    
    pub fn to_string(&self, superpose: bool) -> String {
        let matching_residues = self.matching_residues.iter().map(|x| {
            match x {
                Some((a, b)) => format!("{}{}", *a as char, b),
                None => "_".to_string()
            }
        }).collect::<Vec<String>>().join(",");
        if superpose {
            // print u_matrix by flattening it
            let u_string = self.u_matrix.iter().flat_map(
                |x| x.iter()).map(|&val| format!("{:.4}", val)
            ).collect::<Vec<String>>().join(",");
            // print t_matrix by flattening it
            let t_string = self.t_matrix.iter().map(|&val| format!("{:.4}", val)
            ).collect::<Vec<String>>().join(",");
            // c-alpha coordinates. print x, y, z by concatenating with comma
            let matching_coordinates = self.matching_coordinates.iter().map(|x| {
                format!("{:.4},{:.4},{:.4}", x.x, x.y, x.z)
            }).collect::<Vec<String>>().join(",");
            // Return
            format!(
                "{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}",
                self.id, self.node_count, self.idf, self.rmsd,
                matching_residues, u_string, t_string, matching_coordinates, self.db_key
            )
        } else {
            format!(
                "{}\t{}\t{:.4}\t{:.4}\t{}\t{}", 
                self.id, self.node_count, self.idf, self.rmsd,
                matching_residues, self.db_key
            )
        }
    }
}

impl<'a> fmt::Display for MatchResult<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f, "{}\t{}\t{:.4}\t{:.4}\t{}\t{}",
            self.id, self.node_count, self.idf, self.rmsd,
            self.matching_residues.iter().map(|x| {
                match x {
                    Some((a, b)) => format!("{}{}", *a as char, b),
                    None => "_".to_string()
                }
            }).collect::<Vec<String>>().join(","),
            self.db_key
        )
    }
}

impl<'a> fmt::Debug for MatchResult<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f, "{}\t{}\t{:.4}\t{:.4}\t{}\t{}",
            self.id, self.node_count, self.idf, self.rmsd,
            self.matching_residues.iter().map(|x| {
                match x {
                    Some((a, b)) => format!("{}{}", *a as char, b),
                    None => "_".to_string()
                }
            }).collect::<Vec<String>>().join(","),
            self.db_key
        )
    }
}

// write_fmt
impl<'a> MatchResult<'a> {
    pub fn write_fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f, "{}\t{}\t{:.4}\t{:.4}\t{}\t{}",
            self.id, self.node_count, self.idf, self.rmsd,
            self.matching_residues.iter().map(|x| {
                match x {
                    Some((a, b)) => format!("{}{}", *a as char, b),
                    None => "_".to_string()
                }
            }).collect::<Vec<String>>().join(","),
            self.db_key
        )
    }
}

// Simple Z-score based statistics and scoring
#[derive(Debug, Clone)]
pub struct SimpleStats {
    pub idf_mean: f32,
    pub idf_std: f32,
    pub rmsd_mean: f32,
    pub rmsd_std: f32,
    pub tm_mean: f32,
    pub tm_std: f32,
}

impl SimpleStats {
    pub fn new() -> Self {
        Self {
            idf_mean: 0.0,
            idf_std: 1.0,
            rmsd_mean: 0.0,
            rmsd_std: 1.0,
            tm_mean: 0.0,
            tm_std: 1.0,
        }
    }
    
    /// Calculate statistics from structure results
    pub fn calculate_from_structure_results(results: &[(usize, StructureResult)]) -> Self {
        if results.is_empty() {
            return Self::new();
        }
        
        let idfs: Vec<f32> = results.iter().map(|(_, r)| r.idf).collect();
        let rmsds: Vec<f32> = results.iter().map(|(_, r)| r.min_rmsd_with_max_match).collect();
        let tm_scores: Vec<f32> = results.iter().map(|(_, r)| r.max_tm_score_with_max_match).collect();
        
        Self {
            idf_mean: calculate_mean(&idfs),
            idf_std: calculate_std(&idfs),
            rmsd_mean: calculate_mean(&rmsds),
            rmsd_std: calculate_std(&rmsds),
            tm_mean: calculate_mean(&tm_scores),
            tm_std: calculate_std(&tm_scores),
        }
    }
    
    /// Calculate statistics from match results
    pub fn calculate_from_match_results(results: &[(usize, MatchResult)]) -> Self {
        if results.is_empty() {
            return Self::new();
        }
        
        let idfs: Vec<f32> = results.iter().map(|(_, r)| r.idf).collect();
        let rmsds: Vec<f32> = results.iter().map(|(_, r)| r.rmsd).collect();
        let tm_scores: Vec<f32> = results.iter().map(|(_, r)| r.tm_score).collect();
        
        Self {
            idf_mean: calculate_mean(&idfs),
            idf_std: calculate_std(&idfs),
            rmsd_mean: calculate_mean(&rmsds),
            rmsd_std: calculate_std(&rmsds),
            tm_mean: calculate_mean(&tm_scores),
            tm_std: calculate_std(&tm_scores),
        }
    }
}

// Helper functions for statistics
fn calculate_mean(values: &[f32]) -> f32 {
    if values.is_empty() { return 0.0; }
    values.iter().sum::<f32>() / values.len() as f32
}

fn calculate_std(values: &[f32]) -> f32 {
    if values.len() < 2 { return 1.0; }
    let mean = calculate_mean(values);
    let variance = values.iter()
        .map(|x| (x - mean).powi(2))
        .sum::<f32>() / (values.len() - 1) as f32;
    variance.sqrt().max(0.001) // Avoid division by zero
}

// Z-score calculation methods
impl<'a> StructureResult<'a> {
    /// Calculate simple Z-score: 0.5 * (structural_z) + 0.5 * (idf_z)
    /// where structural_z = 0.5 * (tm_z - rmsd_z)
    pub fn calculate_simple_z_score(&self, stats: &SimpleStats) -> f32 {
        let idf_z = (self.idf - stats.idf_mean) / stats.idf_std;
        let tm_z = (self.max_tm_score_with_max_match - stats.tm_mean) / stats.tm_std;
        let rmsd_z = -((self.min_rmsd_with_max_match - stats.rmsd_mean) / stats.rmsd_std); // Negative for higher is better
        
        let structural_z = 0.5 * (tm_z + rmsd_z);
        let composite_z = 0.5 * structural_z + 0.5 * idf_z;
        
        composite_z
    }
}

impl<'a> MatchResult<'a> {
    /// Calculate simple Z-score for individual matches
    pub fn calculate_simple_z_score(&self, stats: &SimpleStats) -> f32 {
        let idf_z = (self.idf - stats.idf_mean) / stats.idf_std;
        let tm_z = (self.tm_score - stats.tm_mean) / stats.tm_std;
        let rmsd_z = -((self.rmsd - stats.rmsd_mean) / stats.rmsd_std); // Negative for higher is better
        
        let structural_z = 0.5 * (tm_z + rmsd_z);
        let composite_z = 0.5 * structural_z + 0.5 * idf_z;
        
        composite_z
    }
}

pub fn sort_and_print_structure_query_result(
    results: &mut Vec<(usize, StructureResult)>, do_sort_by_rmsd: bool, 
    output_path: &str, query_string: &str, header: bool, verbose: bool,
) {
    if do_sort_by_rmsd {
        if verbose {
            measure_time!(results.par_sort_by(|a, b| {
                // Testing solely by max_tm_score_with_max_match
                b.1.max_tm_score_with_max_match.partial_cmp(&a.1.max_tm_score_with_max_match).unwrap()
                // // Primary by max_matching_node_count, secondary by min_rmsd
                // if a.1.max_matching_node_count != b.1.max_matching_node_count {
                //     b.1.max_matching_node_count.partial_cmp(&a.1.max_matching_node_count).unwrap()
                // } else {
                //     // a.1.min_rmsd_with_max_match.partial_cmp(&b.1.min_rmsd_with_max_match).unwrap()
                //     b.1.max_tm_score_with_max_match.partial_cmp(&a.1.max_tm_score_with_max_match).unwrap()
                // }
            }));
        } else {
            results.par_sort_by(|a, b| {
                b.1.max_tm_score_with_max_match.partial_cmp(&a.1.max_tm_score_with_max_match).unwrap()
                // // Primary by max_matching_node_count, secondary by min_rmsd
                // if a.1.max_matching_node_count != b.1.max_matching_node_count {
                //     b.1.max_matching_node_count.partial_cmp(&a.1.max_matching_node_count).unwrap()
                // } else {
                //     // a.1.min_rmsd_with_max_match.partial_cmp(&b.1.min_rmsd_with_max_match).unwrap()
                //     b.1.max_tm_score_with_max_match.partial_cmp(&a.1.max_tm_score_with_max_match).unwrap()
                // }
            });
        }
    }

    // If output path is not empty, write to file
    if !output_path.is_empty() {
        // Create file
        let file = std::fs::File::create(&output_path).expect(
            &log_msg(FAIL, &format!("Failed to create file: {}", &output_path))
        );
        let mut writer = std::io::BufWriter::new(file);
        if header {
            writer.write_all(format!("{}\n", STRUCTURE_QUERY_RESULT_HEADER).as_bytes()).expect(
                &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
            );
        }
        for (_k, v) in results.iter() {
            writer.write_all(format!("{:?}\t{}\n", v, query_string).as_bytes()).expect(
                &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
            );
        }
    } else {
        if header {
            println!("{}", STRUCTURE_QUERY_RESULT_HEADER);
        }
        // let mut id_container = String::new();
        for (_k, v) in results.iter() {
            println!("{:?}\t{}", v, query_string);
        }
    }
}

pub fn sort_and_print_match_query_result(
    results: &mut Vec<(usize, MatchResult)>, top_n: usize, 
    output_path: &str, query_string: &str, superpose: bool, header: bool, verbose: bool,
    do_sort_by_rmsd: bool,
) {
    // Sort query_count_vec by rmsd
    if verbose {
        measure_time!(results.par_sort_by(|a, b| {
            // Primary by max_matching_node_count, secondary by min_rmsd
            // Testing
            b.1.tm_score.partial_cmp(&a.1.tm_score).unwrap()
            // if a.1.node_count != b.1.node_count {
            //     b.1.node_count.partial_cmp(&a.1.node_count).unwrap()
            // } else {
            //     if do_sort_by_rmsd {
            //         a.1.rmsd.partial_cmp(&b.1.rmsd).unwrap()
            //     } else {
            //         // If not sorting by rmsd, sort by idf
            //         b.1.idf.partial_cmp(&a.1.idf).unwrap()
            //     }
            // }
        }));
    } else {
        results.par_sort_by(|a, b| {
            // Primary by max_matching_node_count, secondary by min_rmsd
            // Testing
            b.1.tm_score.partial_cmp(&a.1.tm_score).unwrap()
            // if a.1.node_count != b.1.node_count {
            //     b.1.node_count.partial_cmp(&a.1.node_count).unwrap()
            // } else {
            //     if do_sort_by_rmsd {
            //         a.1.rmsd.partial_cmp(&b.1.rmsd).unwrap()
            //     } else {
            //         // If not sorting by rmsd, sort by idf
            //         b.1.idf.partial_cmp(&a.1.idf).unwrap()
            //     }
            // }
        });
    }
    // Apply top N filter if top_n is not usize::MAX
    if top_n != usize::MAX {
        if verbose {
            print_log_msg(INFO, &format!("Printing top {} results", top_n));
        }
        results.truncate(top_n);
    }

    // If output path is not empty, write to file
    if !output_path.is_empty() {
        // Create file
        let file = std::fs::File::create(&output_path).expect(
            &log_msg(FAIL, &format!("Failed to create file: {}", &output_path))
        );
        let mut writer = std::io::BufWriter::new(file);
        if header {
            // Print header based on superpose
            if superpose {
                writer.write_all(format!("{}\n", MATCH_QUERY_RESULT_SUPERPOSE_HEADER).as_bytes()).expect(
                    &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
                );
            } else {
                writer.write_all(format!("{}\n", MATCH_QUERY_RESULT_HEADER).as_bytes()).expect(
                    &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
                );
            }
        }
        for (_k, v) in results.iter() {
            writer.write_all(format!("{}\t{}\n", v.to_string(superpose), query_string).as_bytes()).expect(
                &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
            );
        }
    } else {
        if header {
            println!("{}", MATCH_QUERY_RESULT_HEADER);
        }
        // let mut id_container = String::new();
        for (_k, v) in results.iter() {
            println!("{}\t{}", v.to_string(superpose), query_string);
        }
    }
}

// Z-score based sorting functions
pub fn sort_and_print_structure_query_result_z_score(
    results: &mut Vec<(usize, StructureResult)>, 
    output_path: &str, 
    query_string: &str, 
    header: bool, 
    verbose: bool,
) {
    // Calculate statistics from the results
    let stats = if verbose {
        print_log_msg(INFO, "Calculating Z-score statistics for ranking");
        measure_time!(SimpleStats::calculate_from_structure_results(results))
    } else {
        SimpleStats::calculate_from_structure_results(results)
    };
    
    if verbose {
        print_log_msg(INFO, &format!(
            "Stats - IDF: {:.2}±{:.2}, RMSD: {:.2}±{:.2}, TM: {:.3}±{:.3}",
            stats.idf_mean, stats.idf_std, stats.rmsd_mean, stats.rmsd_std, stats.tm_mean, stats.tm_std
        ));
    }
    
    // Sort by Z-score (descending order - higher Z-score is better)
    if verbose {
        measure_time!(results.par_sort_by(|a, b| {
            let z_score_a = a.1.calculate_simple_z_score(&stats);
            let z_score_b = b.1.calculate_simple_z_score(&stats);
            z_score_b.partial_cmp(&z_score_a).unwrap()
        }));
    } else {
        results.par_sort_by(|a, b| {
            let z_score_a = a.1.calculate_simple_z_score(&stats);
            let z_score_b = b.1.calculate_simple_z_score(&stats);
            z_score_b.partial_cmp(&z_score_a).unwrap()
        });
    }

    // Output with Z-scores
    if !output_path.is_empty() {
        let file = std::fs::File::create(&output_path).expect(
            &log_msg(FAIL, &format!("Failed to create file: {}", &output_path))
        );
        let mut writer = std::io::BufWriter::new(file);
        if header {
            writer.write_all(format!("{}\tz_score\n", STRUCTURE_QUERY_RESULT_HEADER).as_bytes()).expect(
                &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
            );
        }
        for (_k, v) in results.iter() {
            let z_score = v.calculate_simple_z_score(&stats);
            writer.write_all(format!("{:?}\t{:.4}\t{}\n", v, z_score, query_string).as_bytes()).expect(
                &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
            );
        }
    } else {
        if header {
            println!("{}\tz_score", STRUCTURE_QUERY_RESULT_HEADER);
        }
        for (_k, v) in results.iter() {
            let z_score = v.calculate_simple_z_score(&stats);
            println!("{:?}\t{:.4}\t{}", v, z_score, query_string);
        }
    }
}

pub fn sort_and_print_match_query_result_z_score(
    results: &mut Vec<(usize, MatchResult)>, 
    top_n: usize,
    output_path: &str, 
    query_string: &str, 
    superpose: bool, 
    header: bool, 
    verbose: bool,
) {
    // Calculate statistics from the results
    let stats = if verbose {
        print_log_msg(INFO, "Calculating Z-score statistics for ranking");
        measure_time!(SimpleStats::calculate_from_match_results(results))
    } else {
        SimpleStats::calculate_from_match_results(results)
    };
    
    if verbose {
        print_log_msg(INFO, &format!(
            "Stats - IDF: {:.2}±{:.2}, RMSD: {:.2}±{:.2}, TM: {:.3}±{:.3}",
            stats.idf_mean, stats.idf_std, stats.rmsd_mean, stats.rmsd_std, stats.tm_mean, stats.tm_std
        ));
    }

    // Sort by Z-score (descending order - higher Z-score is better)
    if verbose {
        measure_time!(results.par_sort_by(|a, b| {
            let z_score_a = a.1.calculate_simple_z_score(&stats);
            let z_score_b = b.1.calculate_simple_z_score(&stats);
            z_score_b.partial_cmp(&z_score_a).unwrap()
        }));
    } else {
        results.par_sort_by(|a, b| {
            let z_score_a = a.1.calculate_simple_z_score(&stats);
            let z_score_b = b.1.calculate_simple_z_score(&stats);
            z_score_b.partial_cmp(&z_score_a).unwrap()
        });
    }

    // Apply top N filter
    if top_n != usize::MAX {
        if verbose {
            print_log_msg(INFO, &format!("Limiting result to top {} matches", top_n));
        }
        results.truncate(top_n);
    }

    // Output with Z-scores
    if !output_path.is_empty() {
        let file = std::fs::File::create(&output_path).expect(
            &log_msg(FAIL, &format!("Failed to create file: {}", &output_path))
        );
        let mut writer = std::io::BufWriter::new(file);
        let header_str = if superpose {
            format!("{}\tz_score", MATCH_QUERY_RESULT_SUPERPOSE_HEADER)
        } else {
            format!("{}\tz_score", MATCH_QUERY_RESULT_HEADER)
        };
        if header {
            writer.write_all(format!("{}\n", header_str).as_bytes()).expect(
                &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
            );
        }
        for (_k, v) in results.iter() {
            let z_score = v.calculate_simple_z_score(&stats);
            writer.write_all(format!("{}\t{:.4}\t{}\n", v.to_string(superpose), z_score, query_string).as_bytes()).expect(
                &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
            );
        }
    } else {
        let header_str = if superpose {
            format!("{}\tz_score", MATCH_QUERY_RESULT_SUPERPOSE_HEADER)
        } else {
            format!("{}\tz_score", MATCH_QUERY_RESULT_HEADER)
        };
        if header {
            println!("{}", header_str);
        }
        for (_k, v) in results.iter() {
            let z_score = v.calculate_simple_z_score(&stats);
            println!("{}\t{:.4}\t{}", v.to_string(superpose), z_score, query_string);
        }
    }
}

// TODO: Need testing