// Result Filtering module
use super::result::{ MatchResult, StructureResult };

pub struct StructureFilter {
    // Filtering parameters that doesn't require residue matching
    pub total_match_count: usize,
    pub covered_node_count: usize,
    pub covered_node_ratio: f32,
    pub covered_edge_count: usize,
    pub covered_edge_ratio: f32,
    pub idf_score: f32,
    pub nres: usize,
    pub plddt: f32,
    // Filtering parameters that require residue matching
    pub max_matching_node_count: usize,
    pub max_matching_node_ratio: f32,
    pub rmsd: f32,
    // Expected number of residues and nodes
    pub expected_node_count: usize,
    pub expected_edge_count: usize,
}

impl StructureFilter {
    pub fn new(
        total_match_count: usize, covered_node_count: usize, 
        covered_node_ratio: f32, covered_edge_count: usize, covered_edge_ratio: f32,
        idf: f32, nres: usize, plddt: f32,
        max_matching_node_count: usize, max_matching_node_ratio: f32,
        rmsd: f32, expected_node_count: usize, expected_edge_count: usize,
    ) -> Self {
        StructureFilter {
            total_match_count: total_match_count,
            covered_node_count,
            covered_node_ratio,
            covered_edge_count,
            covered_edge_ratio,
            idf_score: idf,
            nres,
            plddt,
            max_matching_node_count: max_matching_node_count,
            max_matching_node_ratio: max_matching_node_ratio,
            rmsd,
            expected_node_count,
            expected_edge_count,
        }
    }

    // No filter
    pub fn none() -> Self {
        StructureFilter {
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.0,
            covered_edge_count: 0,
            covered_edge_ratio: 0.0,
            idf_score: 0.0,
            nres: 0,
            plddt: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            rmsd: 0.0,
            expected_node_count: 0,
            expected_edge_count: 0,
        }
    }
    
    // Default filters
    pub fn default(node_count: usize) -> Self {
        let expected_edge_count = node_count * (node_count - 1);
        StructureFilter {
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.8,
            covered_edge_count: 0,
            covered_edge_ratio: 0.4,
            idf_score: 0.0,
            nres: 0,
            plddt: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            rmsd: 0.0,
            expected_node_count: node_count,
            expected_edge_count: expected_edge_count,
        }
    }
    
    
    // Filter single query result
    #[inline]
    pub fn filter_before_matching(&self, result: &StructureResult) -> bool {
        let mut pass = true;
        if self.total_match_count > 0 {
            pass = pass && result.total_match_count >= self.total_match_count;
        }
        if self.covered_node_count > 0 {
            pass = pass && result.node_count >= self.covered_node_count;
        }
        if self.covered_node_ratio > 0.0 {
            pass = pass && result.node_count as f32 / self.expected_node_count as f32 >= self.covered_node_ratio;
        }
        if self.covered_edge_count > 0 {
            pass = pass && result.edge_count >= self.covered_edge_count;
        }
        if self.covered_edge_ratio > 0.0 {
            pass = pass && result.edge_count as f32 / self.expected_edge_count as f32 >= self.covered_edge_ratio;
        }
        if self.idf_score > 0.0 {
            pass = pass && result.idf >= self.idf_score;
        }
        if self.nres > 0 {
            // Number of residues in the query structure
            // Should be less than or equal to the number of residues in the target structure
            pass = pass && result.nres <= self.nres;
        }
        if self.plddt > 0.0 {
            pass = pass && result.plddt >= self.plddt;
        }
        //
        pass
    }

    #[inline]
    pub fn filter_after_matching(&self, result: &StructureResult) -> bool {
        let mut pass = true;
        if self.max_matching_node_count > 0 {
            pass = pass && result.max_matching_node_count >= self.max_matching_node_count;
        }
        if self.max_matching_node_ratio > 0.0 {
            pass = pass && result.max_matching_node_count as f32 / self.expected_node_count as f32 >= self.max_matching_node_ratio;
        }
        if self.rmsd > 0.0 {
            pass = pass && result.min_rmsd_with_max_match <= self.rmsd;
        }
        //
        pass
    }
    
}

pub struct MatchFilter {
    pub node_count: usize,
    pub node_ratio: f32,
    pub avg_idf: f32,
    pub rmsd: f32,
    // Expected number of nodes
    pub expected_node_count: usize,
}

impl MatchFilter {
    pub fn new(node_count: usize, node_ratio: f32, avg_idf: f32, rmsd: f32, expected_node_count: usize) -> Self {
        MatchFilter {
            node_count,
            node_ratio,
            avg_idf,
            rmsd,
            expected_node_count,
        }
    }

        // 
    pub fn none() -> Self {
        MatchFilter {
            node_count: 0,
            node_ratio: 0.0,
            avg_idf: 0.0,
            rmsd: 0.0,
            expected_node_count: 0,
        }
    }
    pub fn default(node_count: usize) -> Self {
        MatchFilter {
            node_count: 0,
            node_ratio: 0.8,
            avg_idf: 0.0,
            rmsd: 1.0, // Default at 1.0
            expected_node_count: node_count,
        }
    }


    // Filter single query result
    #[inline]
    pub fn filter(&self, result: &MatchResult) -> bool {
        let mut pass = true;
        if self.node_count > 0 {
            pass = pass && result.node_count >= self.node_count;
        }
        if self.node_ratio > 0.0 {
            pass = pass && result.node_count as f32 / self.expected_node_count as f32 >= self.node_ratio;
        }
        if self.avg_idf > 0.0 {
            pass = pass && result.avg_idf >= self.avg_idf;
        }
        if self.rmsd > 0.0 {
            pass = pass && result.rmsd <= self.rmsd;
        }
        //
        pass
    }
}

// TODO: Need testing