// Query results (per structure and per match) and their TSV output.

use std::fmt;
use rayon::slice::ParallelSliceMut;

use crate::measure_time;
use crate::prelude::{log_msg, print_log_msg, FAIL, INFO};
use crate::structure::chain_id::{format_chain_residue, residue_list_needs_separator};
use crate::structure::coordinate::Coordinate;
use crate::structure::metrics::StructureSimilarityMetrics;
use crate::utils::formatter::{Column, TsvFormatter, Value, DEFAULT_FLOAT_PRECISION};
use rustc_hash::FxHashMap as HashMap;

use super::ResidueMatch;
use super::sort::{MatchSortStrategy, StructureSortStrategy};

/// Placeholder printed for a query residue that the target did not match.
const UNMATCHED: &str = "_";

/// Render one `matching_residues` field: `A21,A23`, or `AA_250,10_250` for the whole
/// field if any chain is ambiguous without `_` or `force_separator` (`--chain-sep`) is set.
pub fn format_residue_matches(residues: &[ResidueMatch], force_separator: bool) -> String {
    let separator = force_separator
        || residue_list_needs_separator(residues.iter().filter_map(|r| r.as_ref()).map(|(c, _)| c));
    residues.iter().map(|residue| {
        match residue {
            Some((chain, res_ind)) => format_chain_residue(chain, *res_ind, separator),
            None => UNMATCHED.to_string(),
        }
    }).collect::<Vec<String>>().join(",")
}

/// Per-structure result: index counts, then residue matches once retrieved.
pub struct StructureResult<'a> {
    pub tid: &'a str,
    pub nid: usize,
    pub db_key: usize, // Database key for the structure
    pub total_match_count: usize,
    pub node_count: usize,
    pub edge_count: usize,
    pub idf: f32,
    pub nres: usize,
    pub plddt: f32,
    pub matching_residues: Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>, // (residues, RMSD, U, T, Ca coords, metrics, subgraph IDF) from connected components
    pub matching_residues_processed: Vec<(Vec<ResidueMatch>, f32, [[f32; 3]; 3], [f32; 3], Vec<Coordinate>, StructureSimilarityMetrics, f32)>, // Same, after the Ca-distance check
    pub max_matching_node_count: usize,
    pub min_rmsd_with_max_match: f32,
    /// Lowest superposition-free deformation among the best matches
    pub min_drmsd_with_max_match: f32,
}

impl<'a> StructureResult<'a> {
    /// Default ranking score: matched residues^2 x sqrt(IDF) / (1 + RMSD). Query length is
    /// constant within a query, so counts rank the same as coverage.
    /// Chosen in docs/feature_evaluation.md §15.
    pub fn structure_score(&self) -> f32 {
        let matched = self.max_matching_node_count as f32;
        let rmsd = if self.min_rmsd_with_max_match.is_finite() { self.min_rmsd_with_max_match } else { 0.0 };
        matched * matched * self.idf.max(0.0).sqrt() / (1.0 + rmsd)
    }

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
            min_drmsd_with_max_match: 0.0,
        }
    }

    /// One `MatchResult` per match; `skip_ca_dist` uses matches before the Ca-distance check.
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

/// Flatten structure results into match results.
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
        let matching_residues_str = if self.matching_residues_processed.is_empty() {
            "NA".to_string()
        } else {
            self.matching_residues_processed.iter().map(
                |(x, y, _, _, _, _, _)| format!("{}:{:.4}", format_residue_matches(x, false), y)
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

/// One residue-level match in a target structure.
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
    /// IDF scaled by the fraction of query residues this match covers.
    pub fn coverage_idf(&self) -> f32 {
        if self.matching_residues.is_empty() {
            return self.idf;
        }
        self.idf * self.node_count as f32 / self.matching_residues.len() as f32
    }

    /// Default ranking score: evidence (IDF) x coverage^2 x geometric quality (TM-score).
    /// Chosen in docs/feature_evaluation.md §15.
    pub fn match_score(&self) -> f32 {
        let coverage = if self.matching_residues.is_empty() {
            1.0
        } else {
            self.node_count as f32 / self.matching_residues.len() as f32
        };
        let tm_score = if self.metrics.tm_score.is_finite() { self.metrics.tm_score } else { 0.0 };
        self.idf * coverage * coverage * tm_score
    }

    pub fn new(
        tid: &'a str, nid: usize, avg_idf: f32, matching_residues: Vec<ResidueMatch>, rmsd: f32,
        u_matrix: [[f32; 3]; 3], t_matrix: [f32; 3], matching_coordinates: Vec<Coordinate>, db_key: usize, 
        index_size: usize, query_length: usize,
        metrics: StructureSimilarityMetrics,
    ) -> Self {
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
    
    /// Legacy TSV line; `superpose` adds U, T and Ca coordinates.
    pub fn to_string(&self, superpose: bool) -> String {
        let matching_residues = format_residue_matches(&self.matching_residues, false);
        if superpose {
            let u_string = self.u_matrix.iter().flat_map(
                |x| x.iter()).map(|&val| format!("{:.4}", val)
            ).collect::<Vec<String>>().join(",");
            let t_string = self.t_matrix.iter().map(|&val| format!("{:.4}", val)
            ).collect::<Vec<String>>().join(",");
            let matching_coordinates = self.matching_coordinates.iter().map(|x| {
                format!("{:.4},{:.4},{:.4}", x.x, x.y, x.z)
            }).collect::<Vec<String>>().join(",");
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
            format_residue_matches(&self.matching_residues, false),
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
fn build_structure_result_columns<'a>(
    qid: String, query_residues: String, chain_sep: bool
) -> HashMap<&'static str, Column<StructureResult<'a>>> {
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
        Column::new("structure_score", "matched^2 x sqrt(IDF) / (1 + RMSD)", |r: &StructureResult| Value::Float(r.structure_score(), DEFAULT_FLOAT_PRECISION)),
        Column::new("min_drmsd", "Min dRMSD", |r: &StructureResult| Value::Float(r.min_drmsd_with_max_match, DEFAULT_FLOAT_PRECISION)),
        Column::new("matching_residues", "Matching residues with RMSD", move |r: &StructureResult| {
            if r.matching_residues_processed.is_empty() {
                "NA".into()
            } else {
                r.matching_residues_processed.iter().map(
                    |(x, y, _, _, _, _, _)| format!("{}:{:.4}", format_residue_matches(x, chain_sep), y)
                ).collect::<Vec<String>>().join(";").into()
            }
        }),
        Column::new("query_residues", "Query residues", move |_r: &StructureResult| {
            query_residues.clone().into()
        }),
    ].into_iter().map(|col| (col.key, col)).collect()
}

/// Build all available columns for MatchResult
fn build_match_result_columns<'a>(
    qid: String, query_residues: String, chain_sep: bool
) -> HashMap<&'static str, Column<MatchResult<'a>>> {
    vec![
        Column::new("qid", "Query structure ID", move |_r: &MatchResult| qid.clone().into()),
        Column::new("tid", "Target structure ID", |r: &MatchResult| r.tid.into()),
        Column::new("nid", "Numeric ID", |r: &MatchResult| (r.nid as i64).into()),
        Column::new("db_key", "Database key", |r: &MatchResult| (r.db_key as u64).into()),
        Column::new("node_count", "Node count", |r: &MatchResult| (r.node_count as u64).into()),
        Column::new("idf", "IDF score", |r: &MatchResult| Value::Float(r.idf, DEFAULT_FLOAT_PRECISION)),
        Column::new("coverage_idf", "IDF x matched residue fraction", |r: &MatchResult| Value::Float(r.coverage_idf(), DEFAULT_FLOAT_PRECISION)),
        Column::new("match_score", "IDF x coverage^2 x TM-score", |r: &MatchResult| Value::Float(r.match_score(), DEFAULT_FLOAT_PRECISION)),
        Column::new("rmsd", "RMSD", |r: &MatchResult| Value::Float(r.rmsd, DEFAULT_FLOAT_PRECISION)),
        Column::new("e_value", "E-value", |r: &MatchResult| Value::ScientificFloat(r.evalue, 4)),
        Column::new("u_matrix", "Rotation matrix", |r: &MatchResult| Value::Float3DMatrix(r.u_matrix, DEFAULT_FLOAT_PRECISION, ",")),
        Column::new("t_vector", "Translation vector", |r: &MatchResult| Value::Float3DVector(r.t_matrix, DEFAULT_FLOAT_PRECISION, ",")),
        Column::new("matching_residues", "Matching residues", move |r: &MatchResult| {
            format_residue_matches(&r.matching_residues, chain_sep).into()
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
        Column::new("drmsd", "Distance-matrix RMSD", |r: &MatchResult| Value::Float(r.metrics.drmsd, DEFAULT_FLOAT_PRECISION)),
        Column::new("max_dist_deviation", "Max internal distance deviation", |r: &MatchResult| Value::Float(r.metrics.max_dist_deviation, DEFAULT_FLOAT_PRECISION)),
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
pub fn structure_result_formatter<'a>(
    column_keys: &[&str], qid: &str, query_residues: &str, chain_sep: bool
) -> TsvFormatter<StructureResult<'a>> {
    let registry = build_structure_result_columns(qid.to_string(), query_residues.to_string(), chain_sep);
    let columns: Vec<Column<StructureResult>> = column_keys.iter()
        .filter_map(|&key| registry.get(key).cloned())
        .collect();
    TsvFormatter::new(columns)
}

/// Create a TsvFormatter for StructureResult with default columns
pub fn structure_result_default_formatter<'a>(
    qid: &str, query_residues: &str, chain_sep: bool
) -> TsvFormatter<StructureResult<'a>> {
    structure_result_formatter(STRUCTURE_RESULT_DEFAULT_COLUMNS, qid, query_residues, chain_sep)
}

/// Default column keys for MatchResult output
pub const MATCH_RESULT_DEFAULT_COLUMNS: &[&str] = &[
    "tid",
    "node_count",
    "idf",
    "rmsd",
    //"e_value",
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



/// E-value from IDF score `x`, index size `m` and query length `l`.
/// Unvalidated fit; not shown in default output.
pub fn evalue_fitting(x: f32, m: f32, l: f32) -> f64 {
    let x_d = x as f64;
    let m_d = m as f64;
    let l_d = l as f64;

    let mu = 4.2161 * (l_d * 0.0489).exp() + 3.6661;
    let lam = 0.2894 * (l_d * -0.0762).exp() + 0.0316;
    
    let ref_db_size = 10546.0; 
    let search_space_ref = ref_db_size;
    
    let k_val = (lam * mu).exp() / search_space_ref;
    let real_search_space = m_d;
    let e_val_raw = k_val * real_search_space * l_d * (-lam * x_d).exp();

    let e_val = (e_val_raw * real_search_space) / (e_val_raw + real_search_space);

    e_val
}

/// Create a TsvFormatter for MatchResult with specified column keys
pub fn match_result_formatter<'a>(
    column_keys: &[&str], qid: &str, query_residues: &str, chain_sep: bool
) -> TsvFormatter<MatchResult<'a>> {
    let registry = build_match_result_columns(qid.to_string(), query_residues.to_string(), chain_sep);
    let columns: Vec<Column<MatchResult>> = column_keys.iter()
        .filter_map(|&key| registry.get(key).cloned())
        .collect();
    TsvFormatter::new(columns)
}

/// Create a TsvFormatter for MatchResult with default columns
pub fn match_result_default_formatter<'a>(
    qid: &str, query_residues: &str, chain_sep: bool
) -> TsvFormatter<MatchResult<'a>> {
    match_result_formatter(MATCH_RESULT_DEFAULT_COLUMNS, qid, query_residues, chain_sep)
}

/// Create a TsvFormatter for MatchResult with superpose columns
pub fn match_result_superpose_formatter<'a>(
    qid: &str, query_residues: &str, chain_sep: bool
) -> TsvFormatter<MatchResult<'a>> {
    match_result_formatter(MATCH_RESULT_SUPERPOSE_COLUMNS, qid, query_residues, chain_sep)
}

/// Sort structure results and write them as TSV to `output_path` (stdout if empty).
pub fn sort_and_print_structure_query_result(
    results: &mut Vec<(usize, StructureResult)>, 
    output_path: &str, qid: &str, query_residues: &str, columns: Option<&[&str]>, header: bool, verbose: bool,
    sort_strategy: StructureSortStrategy, chain_sep: bool,
) {
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
        Some(cols) => structure_result_formatter(cols, qid, query_residues, chain_sep),
        None => structure_result_default_formatter(qid, query_residues, chain_sep),
    };

    if !output_path.is_empty() {
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

/// Sort match results, keep `top_n`, and write them as TSV to `output_path` (stdout if empty).
pub fn sort_and_print_match_query_result(
    results: &mut Vec<(usize, MatchResult)>, top_n: usize, 
    output_path: &str, qid: &str, query_residues: &str, columns: Option<&[&str]>, superpose: bool, header: bool, verbose: bool,
    sort_strategy: MatchSortStrategy, chain_sep: bool,
) {
    if verbose {
        measure_time!(results.par_sort_by(|a, b| {
            sort_strategy.compare(&a.1, &b.1)
        }));
    } else {
        results.par_sort_by(|a, b| {
            sort_strategy.compare(&a.1, &b.1)
        });
    }
    if top_n != usize::MAX {
        if verbose {
            print_log_msg(INFO, &format!("Printing top {} results", top_n));
        }
        results.truncate(top_n);
    }

    let formatter = match columns {
        Some(cols) => match_result_formatter(cols, qid, query_residues, chain_sep),
        None => {
            if superpose {
                match_result_superpose_formatter(qid, query_residues, chain_sep)
            } else {
                match_result_default_formatter(qid, query_residues, chain_sep)
            }
        }
    };

    if !output_path.is_empty() {
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

#[cfg(test)]
mod default_score_tests {
    use super::*;
    use crate::structure::chain_id::ChainId;

    fn match_result(node_count: usize, query_length: usize, idf: f32, tm_score: f32) -> MatchResult<'static> {
        let residues: Vec<ResidueMatch> = (0..query_length)
            .map(|i| (i < node_count).then(|| (ChainId::from_byte(b'A'), i as u64)))
            .collect();
        let mut metrics = StructureSimilarityMetrics::default();
        metrics.tm_score = tm_score;
        MatchResult::new(
            "t", 0, idf, residues, 0.5, [[0.0; 3]; 3], [0.0; 3], Vec::new(), 0, 1000, query_length, metrics,
        )
    }

    /// The default per-match score: IDF x coverage^2 x TM-score.
    #[test]
    fn match_score_weights_coverage_and_geometry() {
        let full = match_result(4, 4, 10.0, 0.5);
        assert!((full.match_score() - 5.0).abs() < 1e-5, "{}", full.match_score());
        // Half the residues matched: coverage^2 quarters the score
        let half = match_result(2, 4, 10.0, 0.5);
        assert!((half.match_score() - 1.25).abs() < 1e-5, "{}", half.match_score());
        assert!((full.coverage_idf() - 10.0).abs() < 1e-5);
        assert!((half.coverage_idf() - 5.0).abs() < 1e-5);
        // A worse superposition ranks below an equal-coverage match
        assert!(match_result(4, 4, 10.0, 0.9).match_score() > full.match_score());
    }

    /// The default per-structure score: matched^2 x sqrt(IDF) / (1 + RMSD).
    #[test]
    fn structure_score_weights_matched_residues_and_rmsd() {
        let mut r = StructureResult::new("t", 0, 5, 3, 4, 16.0, 100, 90.0, 0);
        r.max_matching_node_count = 3;
        r.min_rmsd_with_max_match = 1.0;
        assert!((r.structure_score() - 18.0).abs() < 1e-4, "{}", r.structure_score());
        let mut worse = StructureResult::new("t", 0, 5, 3, 4, 16.0, 100, 90.0, 0);
        worse.max_matching_node_count = 2;
        worse.min_rmsd_with_max_match = 1.0;
        assert!(worse.structure_score() < r.structure_score());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::chain_id::ChainId;

    fn residue(chain: &str, index: u64) -> ResidueMatch {
        Some((ChainId::from_str(chain), index))
    }

    /// Single-letter chains keep the legacy, parser-compatible spelling.
    #[test]
    fn single_letter_chains_keep_the_legacy_spelling() {
        let residues = vec![
            residue("A", 21), residue("A", 23), residue("A", 25),
            residue("A", 27), residue("A", 32),
        ];
        assert_eq!(format_residue_matches(&residues, false), "A21,A23,A25,A27,A32");
        let mixed = vec![residue("B", 57), residue("B", 102), residue("C", 195)];
        assert_eq!(format_residue_matches(&mixed, false), "B57,B102,C195");
    }

    #[test]
    fn unmatched_residues_stay_a_bare_underscore() {
        let residues = vec![residue("A", 21), None, residue("A", 25)];
        assert_eq!(format_residue_matches(&residues, false), "A21,_,A25");
    }

    /// One ambiguous chain switches the whole field to `_` separators.
    #[test]
    fn ambiguous_chains_switch_the_whole_field() {
        let multi_char = vec![residue("AA", 250), residue("AB", 232)];
        assert_eq!(format_residue_matches(&multi_char, false), "AA_250,AB_232");

        let numeric = vec![residue("10", 250), residue("10", 252)];
        assert_eq!(format_residue_matches(&numeric, false), "10_250,10_252");

        let mixed = vec![residue("A", 21), residue("10", 250)];
        assert_eq!(format_residue_matches(&mixed, false), "A_21,10_250");
    }

    #[test]
    fn chain_sep_forces_the_separated_spelling() {
        let residues = vec![residue("A", 21), None, residue("A", 23)];
        assert_eq!(format_residue_matches(&residues, true), "A_21,_,A_23");
    }

    /// Every emitted spelling parses back as a query.
    #[test]
    fn output_round_trips_through_the_query_parser() {
        use crate::controller::query::parse_query_string;
        let cases: Vec<Vec<ResidueMatch>> = vec![
            vec![residue("A", 21), residue("A", 23), residue("A", 32)],
            vec![residue("B", 57), residue("C", 195)],
            vec![residue("AA", 250), residue("AB", 232)],
            vec![residue("10", 250), residue("1", 3)],
        ];
        for residues in cases {
            let expected: Vec<(ChainId, u64)> = residues.iter()
                .map(|r| r.unwrap()).collect();
            for forced in [false, true] {
                let printed = format_residue_matches(&residues, forced);
                let (parsed, _) = parse_query_string(&printed, ChainId::from_byte(b'A'));
                assert_eq!(parsed, expected, "round trip of {printed:?}");
            }
        }
    }
}
