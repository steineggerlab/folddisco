// QueryResult struct and its implementation

use std::fmt;
use std::collections::HashSet;
use std::io::Write;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::slice::ParallelSliceMut;

use crate::measure_time;
use crate::prelude::{log_msg, print_log_msg, FAIL, INFO};

use super::ResidueMatch;


// TODO: need to think about the column name
pub const STRUCTURE_QUERY_RESULT_HEADER: &str = "id\tidf_score\ttotal_match_count\tnode_count\tedge_count\tmax_node_cov\tmin_rmsd\tnres\tplddt\tmatching_residues\tresidues_rescued\tquery_residues";
pub const MATCH_QUERY_RESULT_HEADER: &str = "id\tnode_count\tavg_idf\trmsd\tmatching_residues\tquery_residues";

pub struct StructureResult<'a> {
    pub id: &'a str,
    pub nid: usize,
    pub total_match_count: usize,
    pub node_count: usize,
    pub edge_count: usize,
    pub idf: f32,
    pub nres: usize,
    pub plddt: f32,
    pub node_set: HashSet<usize>,
    pub edge_set: HashSet<(usize, usize)>,
    pub matching_residues: Vec<(Vec<ResidueMatch>, f32)>, // Match with connected components
    pub matching_residues_processed: Vec<(Vec<ResidueMatch>, f32)>, // Match with c-alpha distances
    pub max_matching_node_count: usize,
    pub min_rmsd_with_max_match: f32,
}

impl<'a> StructureResult<'a> {
    pub fn new(
        id: &'a str, nid: usize, total_match_count: usize, node_count: usize, edge_count: usize,
        idf: f32, nres: usize, plddt: f32, edge: &(usize, usize), 
    ) -> Self {
        let mut node_set = HashSet::new();
        node_set.insert(edge.0);
        node_set.insert(edge.1);
        let mut edge_set = HashSet::new();
        edge_set.insert(*edge);
        Self {
            id,
            nid,
            total_match_count,
            node_count,
            edge_count,
            idf,
            nres,
            plddt,
            node_set: node_set,
            edge_set: edge_set,
            matching_residues: Vec::new(),
            matching_residues_processed: Vec::new(),
            max_matching_node_count: 0,
            min_rmsd_with_max_match: 0.0,
        }
    }
    
    pub fn into_match_query_results(&self) -> Vec<MatchResult> {
        self.matching_residues_processed.iter().enumerate().map(|(i, (residues, rmsd))| {
            // WARNING: NOTE: Thinking of getting the match specific idf by saving the idf for each edge
            let avg_idf = self.idf / self.matching_residues_processed.len() as f32;
            MatchResult::new(
                self.id, i, avg_idf, residues.clone(), *rmsd
            )
        }).collect()
    }
}

pub fn convert_structure_query_result_to_match_query_results<'a>(
    results: &'a [(usize, StructureResult<'a>)]
) -> Vec<(usize, MatchResult<'a>)> {
    results
        .iter()
        .flat_map(|(k, v)| {
            v.into_match_query_results()
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
                |(x, y)| format!("{}:{:.4}", x.iter().map(|x| {
                    match x {
                        // Convert u8 to char
                        Some((a, b)) => format!("{}{}", *a as char, b),
                        None => "_".to_string()
                    }
                }).collect::<Vec<String>>().join(","), y)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.max_matching_node_count, self.min_rmsd_with_max_match,
            self.nres, self.plddt, matching_residues_processed_with_score
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
                |(x, y)| format!("{}:{:.4}", x.iter().map(|x| {
                    match x {
                        Some((a, b)) => format!("{}{}", *a as char, b),
                        None => "_".to_string()
                    }
                }).collect::<Vec<String>>().join(","), y)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.max_matching_node_count, self.min_rmsd_with_max_match,
            self.nres, self.plddt, matching_residues_processed_with_score
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
                |(x, y)| format!("{}:{:.4}", x.iter().map(|x| {
                    match x {
                        Some((a, b)) => format!("{}{}", *a as char, b),
                        None => "_".to_string()
                    }
                }).collect::<Vec<String>>().join(","), y)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}", 
            self.id ,self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.max_matching_node_count, self.min_rmsd_with_max_match,
            self.nres, self.plddt, matching_residues_processed_with_score
        )
    }
}

pub struct MatchResult<'a> {
    pub id: &'a str,
    pub nid: usize,
    pub node_count: usize,
    pub avg_idf: f32,
    pub matching_residues: Vec<ResidueMatch>,
    pub rmsd: f32,
}

impl<'a> MatchResult<'a> {
    pub fn new(
        id: &'a str, nid: usize, avg_idf: f32, matching_residues: Vec<ResidueMatch>, rmsd: f32
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
            node_count,
            avg_idf,
            matching_residues,
            rmsd,
        }
    }
}

impl<'a> fmt::Display for MatchResult<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f, "{}\t{}\t{:.4}\t{:.4}\t{}", 
            self.id, self.node_count, self.avg_idf, self.rmsd,
            self.matching_residues.iter().map(|x| {
                match x {
                    Some((a, b)) => format!("{}{}", *a as char, b),
                    None => "_".to_string()
                }
            }).collect::<Vec<String>>().join(",")
        )
    }
}

impl<'a> fmt::Debug for MatchResult<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f, "{}\t{}\t{:.4}\t{:.4}\t{}", 
            self.id, self.node_count, self.avg_idf, self.rmsd,
            self.matching_residues.iter().map(|x| {
                match x {
                    Some((a, b)) => format!("{}{}", *a as char, b),
                    None => "_".to_string()
                }
            }).collect::<Vec<String>>().join(",")
        )
    }
}

// write_fmt
impl<'a> MatchResult<'a> {
    pub fn write_fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f, "{}\t{}\t{:.4}\t{:.4}\t{}", 
            self.id, self.node_count, self.avg_idf, self.rmsd,
            self.matching_residues.iter().map(|x| {
                match x {
                    Some((a, b)) => format!("{}{}", *a as char, b),
                    None => "_".to_string()
                }
            }).collect::<Vec<String>>().join(",")
        )
    }
}

pub fn sort_and_print_structure_query_result(
    results: &mut Vec<(usize, StructureResult)>, do_sort_by_rmsd: bool, 
    output_path: &str, query_string: &str, header: bool, verbose: bool,
) {
    if do_sort_by_rmsd {
        if verbose {
            measure_time!(results.par_sort_by(|a, b| {
                // Primary by max_matching_node_count, secondary by min_rmsd
                if a.1.max_matching_node_count != b.1.max_matching_node_count {
                    b.1.max_matching_node_count.partial_cmp(&a.1.max_matching_node_count).unwrap()
                } else {
                    a.1.min_rmsd_with_max_match.partial_cmp(&b.1.min_rmsd_with_max_match).unwrap()
                }
            }));
        } else {
            results.par_sort_by(|a, b| {
                // Primary by max_matching_node_count, secondary by min_rmsd
                if a.1.max_matching_node_count != b.1.max_matching_node_count {
                    b.1.max_matching_node_count.partial_cmp(&a.1.max_matching_node_count).unwrap()
                } else {
                    a.1.min_rmsd_with_max_match.partial_cmp(&b.1.min_rmsd_with_max_match).unwrap()
                }
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
    output_path: &str, query_string: &str, header: bool, verbose: bool,
    // do_sort_by_rmsd: bool, WARNING: not implemented yetf
) {
    // Sort query_count_vec by rmsd
    if verbose {
        measure_time!(results.par_sort_by(|a, b| {
            // Primary by max_matching_node_count, secondary by min_rmsd
            if a.1.node_count != b.1.node_count {
                b.1.node_count.partial_cmp(&a.1.node_count).unwrap()
            } else {
                a.1.rmsd.partial_cmp(&b.1.rmsd).unwrap()
            }
        }));
    } else {
        results.par_sort_by(|a, b| {
            // Primary by max_matching_node_count, secondary by min_rmsd
            if a.1.node_count != b.1.node_count {
                b.1.node_count.partial_cmp(&a.1.node_count).unwrap()
            } else {
                a.1.rmsd.partial_cmp(&b.1.rmsd).unwrap()
            }
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
            writer.write_all(format!("{}\n", MATCH_QUERY_RESULT_HEADER).as_bytes()).expect(
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
            println!("{}", MATCH_QUERY_RESULT_HEADER);
        }
        // let mut id_container = String::new();
        for (_k, v) in results.iter() {
            println!("{:?}\t{}", v, query_string);
        }
    }
}

// TODO: Need testing