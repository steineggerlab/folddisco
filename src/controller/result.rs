// QueryResult struct and its implementation

use std::fmt;
use rayon::slice::ParallelSliceMut;

use crate::measure_time;
use crate::prelude::{log_msg, print_log_msg, FAIL, INFO};
use crate::structure::coordinate::Coordinate;
use crate::structure::metrics::StructureSimilarityMetrics;
use crate::utils::formatter::{Column, TsvFormatter, Value, DEFAULT_FLOAT_PRECISION};
use rustc_hash::FxHashMap as HashMap;

use super::ResidueMatch;
use super::sort::{MatchSortStrategy, StructureSortStrategy};

pub struct StructureResult<'a> {
    pub tid: &'a str,
    pub nid: usize,
    pub db_key: usize, // Database key for the structure
    pub total_match_count: usize,
    pub node_count: usize,
    pub edge_count: usize,  // Add edge count field
    pub idf: f32,
    pub nres: usize,
    pub plddt: f32,
    pub matching_residues: Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>, // Match with connected components, with subgraph IDF
    pub matching_residues_processed: Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>, // Match with c-alpha distances, with subgraph IDF
    pub max_matching_node_count: usize,
    pub min_rmsd_with_max_match: f32,
}

impl<'a> StructureResult<'a> {
    pub fn new(
        tid: &'a str, nid: usize, total_match_count: usize, node_count: usize, edge_count: usize,
        idf: f32, nres: usize, plddt: f32, db_key: usize
    ) -> Self {
        Self {
            tid,
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
        }
    }

    pub fn into_match_query_results(&self, skip_ca_dist: bool, index_size: usize, query_length: usize) -> Vec<MatchResult<'_>> {
        match skip_ca_dist {
            false => self.matching_residues_processed.iter().enumerate().map(|(i, (residues, rmsd, u_matrix, t_matrix, matching_coordinates, metrics, subgraph_idf))| {
                MatchResult::new(
                    self.tid, i, *subgraph_idf, residues.clone(), *rmsd,
                    *u_matrix, *t_matrix, matching_coordinates.clone(), self.db_key, index_size, query_length, metrics.clone()
                )
            }).collect(),
            true => self.matching_residues.iter().enumerate().map(|(i, (residues, rmsd, u_matrix, t_matrix, matching_coordinates, metrics, subgraph_idf))| {
                MatchResult::new(
                    self.tid, i, *subgraph_idf, residues.clone(), *rmsd,
                    *u_matrix, *t_matrix, matching_coordinates.clone(), self.db_key, index_size, query_length, metrics.clone()
                )
            }).collect(),
        }
    }
}

pub fn convert_structure_query_result_to_match_query_results<'a>(
    results: &'a [(usize, StructureResult<'a>)], skip_ca_dist: bool, index_size: usize, query_length: usize
) -> Vec<(usize, MatchResult<'a>)> {
    results
        .iter()
        .flat_map(|(k, v)| {
            v.into_match_query_results(skip_ca_dist, index_size, query_length)
             .into_iter()
             .map(|x| (*k, x))
        })
        .collect()
}

impl<'a> fmt::Display for StructureResult<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use the default formatter logic inline for Display
        let matching_residues_str = if self.matching_residues_processed.is_empty() {
            "NA".to_string()
        } else {
            self.matching_residues_processed.iter().map(
                |(x, y, _, _, _, _, _)| format!("{}:{:.4}", x.iter().map(|x| {
                    match x {
                        Some((a, b)) => format!("{}{}", *a as char, b),
                        None => "_".to_string()
                    }
                }).collect::<Vec<String>>().join(","), y)
            ).collect::<Vec<String>>().join(";")
        };
        write!(
            f, "{}\t{:.4}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}\t{}", 
            self.tid, self.idf, self.total_match_count, self.node_count, self.edge_count,
            self.max_matching_node_count, self.min_rmsd_with_max_match,
            self.nres, self.plddt, matching_residues_str, self.db_key
        )
    }
}

impl<'a> fmt::Debug for StructureResult<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

pub struct MatchResult<'a> {
    pub tid: &'a str,
    pub nid: usize,
    pub db_key: usize, // Database key for the structure
    pub node_count: usize,
    pub idf: f32,
    pub matching_residues: Vec<ResidueMatch>,
    pub rmsd: f32,
    pub evalue: f64,
    pub u_matrix: [[f32; 3]; 3],
    pub t_matrix: [f32; 3],
    pub matching_coordinates: Vec<Coordinate>,
    pub metrics: StructureSimilarityMetrics,
}

impl<'a> MatchResult<'a> {
    pub fn new(
        tid: &'a str, nid: usize, avg_idf: f32, matching_residues: Vec<ResidueMatch>, rmsd: f32,
        u_matrix: [[f32; 3]; 3], t_matrix: [f32; 3], matching_coordinates: Vec<Coordinate>, db_key: usize, 
        index_size: usize, query_length: usize,
        metrics: StructureSimilarityMetrics,
    ) -> Self {
        //
        let node_count = matching_residues.iter().map(|x| {
            match x {
                Some(_) => 1,
                None => 0
            }
        }).sum();

        let evalue =  evalue_fitting(avg_idf, index_size as f32,query_length as f32);
        
        Self {
            tid,
            nid,
            db_key,
            node_count,
            idf: avg_idf,
            matching_residues,
            rmsd,
            evalue,
            u_matrix,
            t_matrix,
            matching_coordinates,
            metrics,
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
                "{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}",
                self.tid, self.node_count, self.idf, self.rmsd, self.evalue,
                matching_residues, u_string, t_string, matching_coordinates, self.db_key, self.metrics
            )
        } else {
            format!(
                "{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}",
                self.tid, self.node_count, self.idf, self.rmsd, self.evalue,
                matching_residues, self.db_key, self.metrics
            )
        }
    }
}

impl<'a> fmt::Display for MatchResult<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f, "{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}",
            self.tid, self.node_count, self.idf, self.rmsd, self.evalue,
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
        fmt::Display::fmt(self, f)
    }
}



// Column registries and defaults


/// Build all available columns for StructureResult
fn build_structure_result_columns<'a>(qid: String, query_residues: String) -> HashMap<&'static str, Column<StructureResult<'a>>> {
    vec![
        Column::new("qid", "Query structure ID", move |_r: &StructureResult| qid.clone().into()),
        Column::new("tid", "Target structure ID", |r: &StructureResult| r.tid.into()),
        Column::new("nid", "Numeric ID", |r: &StructureResult| (r.nid as u64).into()),
        Column::new("db_key", "Database key", |r: &StructureResult| (r.db_key as u64).into()),
        Column::new("total_match_count", "Total match count", |r: &StructureResult| (r.total_match_count as u64).into()),
        Column::new("node_count", "Node count", |r: &StructureResult| (r.node_count as u64).into()),
        Column::new("edge_count", "Edge count", |r: &StructureResult| (r.edge_count as u64).into()),
        Column::new("idf", "IDF score", |r: &StructureResult| Value::Float(r.idf, DEFAULT_FLOAT_PRECISION)),
        Column::new("nres", "Number of residues", |r: &StructureResult| (r.nres as u64).into()),
        Column::new("plddt", "pLDDT score", |r: &StructureResult| Value::Float(r.plddt, 2)),
        Column::new("max_node_cov", "Max node coverage", |r: &StructureResult| (r.max_matching_node_count as u64).into()),
        Column::new("min_rmsd", "Min RMSD", |r: &StructureResult| Value::Float(r.min_rmsd_with_max_match, DEFAULT_FLOAT_PRECISION)),
        Column::new("matching_residues", "Matching residues with RMSD", |r: &StructureResult| {
            if r.matching_residues_processed.is_empty() {
                "NA".into()
            } else {
                r.matching_residues_processed.iter().map(
                    |(x, y, _, _, _, _, _)| format!("{}:{:.4}", x.iter().map(|x| {
                        match x {
                            Some((a, b)) => format!("{}{}", *a as char, b),
                            None => "_".to_string()
                        }
                    }).collect::<Vec<String>>().join(","), y)
                ).collect::<Vec<String>>().join(";").into()
            }
        }),
        Column::new("query_residues", "Query residues", move |_r: &StructureResult| {
            query_residues.clone().into()
        }),
    ].into_iter().map(|col| (col.key, col)).collect()
}

/// Build all available columns for MatchResult
fn build_match_result_columns<'a>(qid: String, query_residues: String) -> HashMap<&'static str, Column<MatchResult<'a>>> {
    vec![
        Column::new("qid", "Query structure ID", move |_r: &MatchResult| qid.clone().into()),
        Column::new("tid", "Target structure ID", |r: &MatchResult| r.tid.into()),
        Column::new("nid", "Numeric ID", |r: &MatchResult| (r.nid as i64).into()),
        Column::new("db_key", "Database key", |r: &MatchResult| (r.db_key as u64).into()),
        Column::new("node_count", "Node count", |r: &MatchResult| (r.node_count as u64).into()),
        Column::new("idf", "IDF score", |r: &MatchResult| Value::Float(r.idf, DEFAULT_FLOAT_PRECISION)),
        Column::new("rmsd", "RMSD", |r: &MatchResult| Value::Float(r.rmsd, DEFAULT_FLOAT_PRECISION)),
        Column::new("e_value", "E-value", |r: &MatchResult| Value::ScientificFloat(r.evalue, 4)),
        Column::new("u_matrix", "Rotation matrix", |r: &MatchResult| Value::Float3DMatrix(r.u_matrix, DEFAULT_FLOAT_PRECISION, ",")),
        Column::new("t_vector", "Translation vector", |r: &MatchResult| Value::Float3DVector(r.t_matrix, DEFAULT_FLOAT_PRECISION, ",")),
        Column::new("matching_residues", "Matching residues", |r: &MatchResult| {
            r.matching_residues.iter().map(|x| {
                match x {
                    Some((a, b)) => format!("{}{}", *a as char, b),
                    None => "_".to_string()
                }
            }).collect::<Vec<String>>().join(",").into()
        }),
        Column::new("matching_coordinates", "Matching C-alpha coordinates", |r: &MatchResult| {
            Value::FloatVector(
                r.matching_coordinates.iter().flat_map(|c| vec![c.x, c.y, c.z]).collect(),
                4,
                ","
            )
        }),
        Column::new("query_residues", "Query residues", move |_r: &MatchResult| {
            query_residues.clone().into()
        }),
        Column::new("tm_score", "TM-score", |r: &MatchResult| Value::Float(r.metrics.tm_score, DEFAULT_FLOAT_PRECISION)),
        Column::new("gdt_ts", "GDT-TS", |r: &MatchResult| Value::Float(r.metrics.gdt_ts, DEFAULT_FLOAT_PRECISION)),
        Column::new("gdt_ha", "GDT-HA", |r: &MatchResult| Value::Float(r.metrics.gdt_ha, DEFAULT_FLOAT_PRECISION)),
        Column::new("chamfer_distance", "Chamfer distance", |r: &MatchResult| Value::Float(r.metrics.chamfer_distance, DEFAULT_FLOAT_PRECISION)),
        Column::new("hausdorff_distance", "Hausdorff distance", |r: &MatchResult| Value::Float(r.metrics.hausdorff_distance, DEFAULT_FLOAT_PRECISION)),
    ].into_iter().map(|col| (col.key, col)).collect()
}


/// Default column keys for StructureResult output
pub const STRUCTURE_RESULT_DEFAULT_COLUMNS: &[&str] = &[
    "tid",
    "idf",
    "total_match_count",
    "node_count",
    "edge_count",
    "max_node_cov",
    "min_rmsd",
    "nres",
    "plddt",
    "matching_residues",
    "db_key",
    "query_residues",
];

/// Create a TsvFormatter for StructureResult with specified column keys
pub fn structure_result_formatter<'a>(column_keys: &[&str], qid: &str, query_residues: &str) -> TsvFormatter<StructureResult<'a>> {
    let registry = build_structure_result_columns(qid.to_string(), query_residues.to_string());
    let columns: Vec<Column<StructureResult>> = column_keys.iter()
        .filter_map(|&key| registry.get(key).cloned())
        .collect();
    TsvFormatter::new(columns)
}

/// Create a TsvFormatter for StructureResult with default columns
pub fn structure_result_default_formatter<'a>(qid: &str, query_residues: &str) -> TsvFormatter<StructureResult<'a>> {
    structure_result_formatter(STRUCTURE_RESULT_DEFAULT_COLUMNS, qid, query_residues)
}

/// Default column keys for MatchResult output
pub const MATCH_RESULT_DEFAULT_COLUMNS: &[&str] = &[
    "tid",
    "node_count",
    "idf",
    "rmsd",
    "e_value",
    "matching_residues",
    "query_residues",
];

/// Column keys for MatchResult output with superpose/web mode
pub const MATCH_RESULT_SUPERPOSE_COLUMNS: &[&str] = &[
    "tid",
    "node_count",
    "idf",
    "rmsd",
    "matching_residues",
    "u_matrix",
    "t_vector",
    "matching_coordinates",
    "db_key",
    "query_residues",
];



//Create a evalue fitting function to compute evalues based on IDF score
pub fn evalue_fitting(x: f32, m: f32, l: f32) -> f64 {
    // TODO: Finalize fitting.
    // WARNING: This is not validated yet.
    // x: score, m: index size, l: query residue length 
    let x_d = x as f64;
    let m_d = m as f64;
    let l_d = l as f64;

    let mu = 4.2161 * (l_d * 0.0489).exp() + 3.6661;
    let lam = 0.2894 * (l_d * -0.0762).exp() + 0.0316;
    
    let ref_db_size = 10546.0; 
    let search_space_ref = ref_db_size;    // let search_space_ref = ref_db_size;
    
    let k_val = (lam * mu).exp() / search_space_ref;
    let real_search_space = m_d;
    let e_val_raw = k_val * real_search_space * l_d * (-lam * x_d).exp();

    let e_val = (e_val_raw * real_search_space) / (e_val_raw + real_search_space);

    e_val as f32
}

/// Create a TsvFormatter for MatchResult with specified column keys
pub fn match_result_formatter<'a>(column_keys: &[&str], qid: &str, query_residues: &str) -> TsvFormatter<MatchResult<'a>> {
    let registry = build_match_result_columns(qid.to_string(), query_residues.to_string());
    let columns: Vec<Column<MatchResult>> = column_keys.iter()
        .filter_map(|&key| registry.get(key).cloned())
        .collect();
    TsvFormatter::new(columns)
}

/// Create a TsvFormatter for MatchResult with default columns
pub fn match_result_default_formatter<'a>(qid: &str, query_residues: &str) -> TsvFormatter<MatchResult<'a>> {
    match_result_formatter(MATCH_RESULT_DEFAULT_COLUMNS, qid, query_residues)
}

/// Create a TsvFormatter for MatchResult with superpose columns
pub fn match_result_superpose_formatter<'a>(qid: &str, query_residues: &str) -> TsvFormatter<MatchResult<'a>> {
    match_result_formatter(MATCH_RESULT_SUPERPOSE_COLUMNS, qid, query_residues)
}

pub fn sort_and_print_structure_query_result(
    results: &mut Vec<(usize, StructureResult)>, 
    output_path: &str, qid: &str, query_residues: &str, columns: Option<&[&str]>, header: bool, verbose: bool,
    sort_strategy: StructureSortStrategy,
) {
    // Sort using the strategy
    if verbose {
        measure_time!(results.par_sort_by(|a, b| {
            sort_strategy.compare(&a.1, &b.1)
        }));
    } else {
        results.par_sort_by(|a, b| {
            sort_strategy.compare(&a.1, &b.1)
        });
    }

    let formatter = match columns {
        Some(cols) => structure_result_formatter(cols, qid, query_residues),
        None => structure_result_default_formatter(qid, query_residues),
    };

    // If output path is not empty, write to file
    if !output_path.is_empty() {
        // Create file
        let file = std::fs::File::create(&output_path).expect(
            &log_msg(FAIL, &format!("Failed to create file: {}", &output_path))
        );
        let mut writer = std::io::BufWriter::new(file);
        if header {
            formatter.write_header(&mut writer).expect(
                &log_msg(FAIL, &format!("Failed to write header to file: {}", &output_path))
            );
        }
        for (_k, v) in results.iter() {
            formatter.write_record(&mut writer, v).expect(
                &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
            );
        }
    } else {
        if header {
            let mut stdout = std::io::stdout();
            formatter.write_header(&mut stdout).expect("Failed to write header to stdout");
        }
        for (_k, v) in results.iter() {
            let mut stdout = std::io::stdout();
            formatter.write_record(&mut stdout, v).expect("Failed to write to stdout");
        }
    }
}

pub fn sort_and_print_match_query_result(
    results: &mut Vec<(usize, MatchResult)>, top_n: usize, 
    output_path: &str, qid: &str, query_residues: &str, columns: Option<&[&str]>, superpose: bool, header: bool, verbose: bool,
    sort_strategy: MatchSortStrategy,
) {
    // Sort using the strategy
    if verbose {
        measure_time!(results.par_sort_by(|a, b| {
            sort_strategy.compare(&a.1, &b.1)
        }));
    } else {
        results.par_sort_by(|a, b| {
            sort_strategy.compare(&a.1, &b.1)
        });
    }
    // Apply top N filter if top_n is not usize::MAX
    if top_n != usize::MAX {
        if verbose {
            print_log_msg(INFO, &format!("Printing top {} results", top_n));
        }
        results.truncate(top_n);
    }

    let formatter = match columns {
        Some(cols) => match_result_formatter(cols, qid, query_residues),
        None => {
            if superpose {
                match_result_superpose_formatter(qid, query_residues)
            } else {
                match_result_default_formatter(qid, query_residues)
            }
        }
    };

    // If output path is not empty, write to file
    if !output_path.is_empty() {
        // Create file
        let file = std::fs::File::create(&output_path).expect(
            &log_msg(FAIL, &format!("Failed to create file: {}", &output_path))
        );
        let mut writer = std::io::BufWriter::new(file);
        if header {
            formatter.write_header(&mut writer).expect(
                &log_msg(FAIL, &format!("Failed to write header to file: {}", &output_path))
            );
        }
        for (_k, v) in results.iter() {
            formatter.write_record(&mut writer, v).expect(
                &log_msg(FAIL, &format!("Failed to write to file: {}", &output_path))
            );
        }
    } else {
        if header {
            let mut stdout = std::io::stdout();
            formatter.write_header(&mut stdout).expect("Failed to write header to stdout");
        }
        for (_k, v) in results.iter() {
            let mut stdout = std::io::stdout();
            formatter.write_record(&mut stdout, v).expect("Failed to write to stdout");
        }
    }
}

// TODO: Need testing
