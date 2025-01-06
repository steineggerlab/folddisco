// Result Filtering module
use crate::controller::rank::QueryResult;

// NOTE: Limiting result with top N results is done in CLI
// This module is for filtering single query result

pub struct QueryResultFilter {
    // Filtering parameters that doesn't require residue matching
    pub total_hash_match_count: usize,
    pub total_hash_match_ratio: f32,
    pub node_covered_by_hash_count: usize,
    pub node_covered_by_hash_ratio: f32,
    pub edge_covered_by_hash_count: usize,
    pub edge_covered_by_hash_ratio: f32,
    pub idf: f32,
    pub nres: usize,
    pub plddt: f32,
    // Filtering parameters that require residue matching
    pub node_covered_by_graph_count: usize,
    pub node_covered_by_graph_ratio: f32,
    pub rmsd: f32,
}

impl QueryResultFilter {
    pub fn new() -> Self {
        QueryResultFilter {
            total_hash_match_count: 0,
            total_hash_match_ratio: 0.0,
            node_covered_by_hash_count: 0,
            node_covered_by_hash_ratio: 0.0,
            edge_covered_by_hash_count: 0,
            edge_covered_by_hash_ratio: 0.0,
            idf: 0.0,
            nres: 0,
            plddt: 0.0,
            node_covered_by_graph_count: 0,
            node_covered_by_graph_ratio: 0.0,
            rmsd: 0.0,
        }
    }
    // Filter single query result
    pub fn filter(&self, result: &QueryResult) -> bool {
        // let mut pass = true;
        // if self.total_hash_match_count > 0 {
        //     pass = pass && result.total_hash_match_count >= self.total_hash_match_count;
        // }
        todo!()
    }

}