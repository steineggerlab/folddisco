// Result Filtering module
use super::result::{ MatchResult, StructureResult };

pub struct StructureFilter {
    // Filtering parameters that doesn't require residue matching
    pub total_match_count: usize,
    pub covered_node_count: usize,
    pub covered_node_ratio: f32,
    pub idf_score: f32,
    pub nres: usize,
    pub plddt: f32,
    // Filtering parameters that require residue matching
    pub max_matching_node_count: usize,
    pub max_matching_node_ratio: f32,
    pub rmsd: f32,
    // Expected number of residues and nodes
    pub expected_node_count: usize,
}

impl StructureFilter {
    pub fn new(
        total_match_count: usize, covered_node_count: usize, 
        covered_node_ratio: f32, idf: f32, nres: usize, plddt: f32,
        max_matching_node_count: usize, max_matching_node_ratio: f32,
        rmsd: f32, expected_node_count: usize,
    ) -> Self {
        StructureFilter {
            total_match_count: total_match_count,
            covered_node_count,
            covered_node_ratio,
            idf_score: idf,
            nres,
            plddt,
            max_matching_node_count: max_matching_node_count,
            max_matching_node_ratio: max_matching_node_ratio,
            rmsd,
            expected_node_count,
        }
    }

    // No filter
    pub fn none() -> Self {
        StructureFilter {
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.0,
            idf_score: 0.0,
            nres: 0,
            plddt: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            rmsd: 0.0,
            expected_node_count: 0,
        }
    }
    
    // Default filters
    pub fn default(node_count: usize) -> Self {
        StructureFilter {
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.8,
            idf_score: 0.0,
            nres: 0,
            plddt: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            rmsd: 0.0,
            expected_node_count: node_count,
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
    // Metrics from StructureSimilarityMetrics
    pub tm_score: f32,
    pub tm_score_strict: f32,
    pub gdt_ts: f32,
    pub gdt_ha: f32,
    pub gdt_strict: f32,
    pub chamfer_distance: f32,
    pub hausdorff_distance: f32,
    pub rmsd: f32,
    // Expected number of nodes
    pub expected_node_count: usize,
}

impl MatchFilter {
    pub fn new(
        node_count: usize, node_ratio: f32, avg_idf: f32, 
        tm_score: f32, tm_score_strict: f32,
        gdt_ts: f32, gdt_ha: f32, gdt_strict: f32,
        chamfer_distance: f32, hausdorff_distance: f32, rmsd: f32,
        expected_node_count: usize
    ) -> Self {
        MatchFilter {
            node_count,
            node_ratio,
            avg_idf,
            tm_score,
            tm_score_strict,
            gdt_ts,
            gdt_ha,
            gdt_strict,
            chamfer_distance,
            hausdorff_distance,
            rmsd,
            expected_node_count,
        }
    }

    // No filter
    pub fn none() -> Self {
        MatchFilter {
            node_count: 0,
            node_ratio: 0.0,
            avg_idf: 0.0,
            tm_score: 0.0,
            tm_score_strict: 0.0,
            gdt_ts: 0.0,
            gdt_ha: 0.0,
            gdt_strict: 0.0,
            chamfer_distance: 0.0,
            hausdorff_distance: 0.0,
            rmsd: 0.0,
            expected_node_count: 0,
        }
    }
    
    pub fn default(node_count: usize) -> Self {
        MatchFilter {
            node_count: 0,
            node_ratio: 0.8,
            avg_idf: 0.0,
            tm_score: 0.0,
            tm_score_strict: 0.0,
            gdt_ts: 0.0,
            gdt_ha: 0.0,
            gdt_strict: 0.0,
            chamfer_distance: 0.0,
            hausdorff_distance: 0.0,
            rmsd: 1.0, // Default at 1.0
            expected_node_count: node_count,
        }
    }


    // Filter single query result
    #[inline]
    pub fn filter(&self, result: &MatchResult) -> bool {
        let mut pass = true;
        
        // Node-based filters
        if self.node_count > 0 {
            pass = pass && result.node_count >= self.node_count;
        }
        if self.node_ratio > 0.0 {
            pass = pass && result.node_count as f32 / self.expected_node_count as f32 >= self.node_ratio;
        }
        if self.avg_idf > 0.0 {
            pass = pass && result.idf >= self.avg_idf;
        }
        
        // Metrics-based filters (using StructureSimilarityMetrics)
        // Higher is better metrics (>=)
        if self.tm_score > 0.0 {
            pass = pass && result.metrics.tm_score >= self.tm_score;
        }
        if self.tm_score_strict > 0.0 {
            pass = pass && result.metrics.tm_score_strict >= self.tm_score_strict;
        }
        if self.gdt_ts > 0.0 {
            pass = pass && result.metrics.gdt_ts >= self.gdt_ts;
        }
        if self.gdt_ha > 0.0 {
            pass = pass && result.metrics.gdt_ha >= self.gdt_ha;
        }
        if self.gdt_strict > 0.0 {
            pass = pass && result.metrics.gdt_strict >= self.gdt_strict;
        }
        
        // Lower is better metrics (<=)
        if self.chamfer_distance > 0.0 {
            pass = pass && result.metrics.chamfer_distance <= self.chamfer_distance;
        }
        if self.hausdorff_distance > 0.0 {
            pass = pass && result.metrics.hausdorff_distance <= self.hausdorff_distance;
        }
        if self.rmsd > 0.0 {
            pass = pass && result.metrics.rmsd <= self.rmsd;
        }
        
        pass
    }
}

// TODO: Need testing