// Result Filtering module
use super::query_result::{ MatchQueryResult, StructureQueryResult };

// NOTE: Limiting result with top N results is done in CLI
// This module is for filtering single query result

// pub struct QueryResult<'a> {
//     pub id: &'a str,
//     pub nid: usize,
//     pub total_match_count: usize,
//     pub node_count: usize,
//     pub edge_count: usize,
//     pub idf: f32,
//     pub nres: usize,
//     pub plddt: f32,
//     pub node_set: HashSet<usize>,
//     pub edge_set: HashSet<(usize, usize)>,
//     pub matching_residues: Vec<(Vec<ResidueMatch>, f32)>, // Match with connected components
//     pub matching_residues_processed: Vec<(Vec<ResidueMatch>, f32)>, // Match with c-alpha distances
// }

pub struct StructureQueryResultFilter {
    // Filtering parameters that doesn't require residue matching
    pub total_hash_match_count: usize,
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

impl StructureQueryResultFilter {
    pub fn new() -> Self {
        StructureQueryResultFilter {
            total_hash_match_count: 0,
            node_covered_by_hash_count: 0,
            node_covered_by_hash_ratio: 0.0,
            edge_covered_by_hash_count: 0,
            edge_covered_by_hash_ratio: 0.0,
            idf: 0.0,
            nres: usize::MAX,
            plddt: 0.0,
            node_covered_by_graph_count: 0,
            node_covered_by_graph_ratio: 0.0,
            rmsd: 0.0,
        }
    }
    // Filter single query result
    pub fn filter(&self, result: &StructureQueryResult) -> bool {
        let mut pass = true;
        if self.total_hash_match_count > 0 {
            pass = pass && result.total_match_count >= self.total_hash_match_count;
        }
        todo!()
    }

}

// pub struct MatchQueryResult<'a> {
//     pub id: &'a str,
//     pub nid: usize,
//     pub node_count: usize,
//     pub matching_residues: Vec<ResidueMatch>,
//     pub rmsd: f32,
// }

pub struct MatchQueryResultFilter {
    pub node_count: usize,
    pub avg_idf: f32,
    pub rmsd: f32,
}

impl MatchQueryResultFilter {
    pub fn new() -> Self {
        MatchQueryResultFilter {
            node_count: 0,
            avg_idf: 0.0,
            rmsd: 0.0,
        }
    }
    // Filter single query result
    pub fn filter(&self, result: &MatchQueryResult) -> bool {
        // let mut pass = true;
        // if self.node_count > 0 {
        //     pass = pass && result.node_count >= self.node_count;
        // }
        todo!()
    }
}