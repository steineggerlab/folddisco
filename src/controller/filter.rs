// Result filters. A cutoff of 0 (or f64::MAX for evalue) disables that check.
use super::result::{ MatchResult, StructureResult };

/// Per-structure cutoffs, applied before and after residue matching.
pub struct StructureFilter {
    // Checked before matching
    pub total_match_count: usize,
    pub covered_node_count: usize,
    pub covered_node_ratio: f32,
    pub idf_per_structure: f32,
    pub nres: usize,
    pub plddt: f32,
    // Checked after matching
    pub max_matching_node_count: usize,
    pub max_matching_node_ratio: f32,
    pub rmsd: f32,
    /// Superposition-free deformation cutoff, for non-rigid matches
    pub drmsd: f32,
    /// Query residue count; denominator of the ratio cutoffs
    pub expected_node_count: usize,
}

impl StructureFilter {
    pub fn new(
        total_match_count: usize, covered_node_count: usize,
        covered_node_ratio: f32, idf_per_structure: f32, nres: usize, plddt: f32,
        max_matching_node_count: usize, max_matching_node_ratio: f32,
        rmsd: f32, drmsd: f32, expected_node_count: usize,
    ) -> Self {
        StructureFilter {
            total_match_count: total_match_count,
            covered_node_count,
            covered_node_ratio,
            idf_per_structure,
            nres,
            plddt,
            max_matching_node_count: max_matching_node_count,
            max_matching_node_ratio: max_matching_node_ratio,
            rmsd,
            drmsd,
            expected_node_count,
        }
    }

    /// No filtering.
    pub fn none() -> Self {
        StructureFilter {
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.0,
            idf_per_structure: 0.0,
            nres: 0,
            plddt: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            rmsd: 0.0,
            drmsd: 0.0,
            expected_node_count: 0,
        }
    }

    /// Require 80% node coverage.
    pub fn default(node_count: usize) -> Self {
        StructureFilter {
            total_match_count: 0,
            covered_node_count: 0,
            covered_node_ratio: 0.8,
            idf_per_structure: 0.0,
            nres: 0,
            plddt: 0.0,
            max_matching_node_count: 0,
            max_matching_node_ratio: 0.0,
            rmsd: 0.0,
            drmsd: 0.0,
            expected_node_count: node_count,
        }
    }
    
    
    /// Checks that need only index counts.
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
        if self.idf_per_structure > 0.0 {
            pass = pass && result.idf >= self.idf_per_structure;
        }
        if self.nres > 0 {
            // Maximum target length (--num-residue)
            pass = pass && result.nres <= self.nres;
        }
        if self.plddt > 0.0 {
            pass = pass && result.plddt >= self.plddt;
        }
        pass
    }

    /// Checks that need residue matching results.
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
        if self.drmsd > 0.0 {
            pass = pass && result.min_drmsd_with_max_match <= self.drmsd;
        }
        pass
    }
    
}

/// Per-match cutoffs.
pub struct MatchFilter {
    pub node_count: usize,
    pub node_ratio: f32,
    pub idf_per_match: f32,
    pub evalue: f64,
    pub rmsd: f32,
    pub tm_score: f32,
    pub gdt_ts: f32,
    pub gdt_ha: f32,
    pub chamfer_distance: f32,
    pub hausdorff_distance: f32,
    /// Superposition-free deformation cutoff, for non-rigid matches
    pub drmsd: f32,
    /// Query residue count; denominator of `node_ratio`
    pub expected_node_count: usize,
}

impl MatchFilter {
    pub fn new(
        node_count: usize, node_ratio: f32, idf_per_match: f32, evalue: f64,
        rmsd: f32, tm_score: f32, gdt_ts: f32, gdt_ha: f32,
        chamfer_distance: f32, hausdorff_distance: f32, drmsd: f32,
        expected_node_count: usize
    ) -> Self {
        MatchFilter {
            node_count,
            node_ratio,
            idf_per_match,
            evalue,
            rmsd,
            tm_score,
            gdt_ts,
            gdt_ha,
            chamfer_distance,
            hausdorff_distance,
            drmsd,
            expected_node_count,
        }
    }

    /// No filtering.
    pub fn none() -> Self {
        MatchFilter {
            node_count: 0,
            node_ratio: 0.0,
            idf_per_match: 0.0,
            evalue: f64::MAX,
            rmsd: 0.0,
            tm_score: 0.0,
            gdt_ts: 0.0,
            gdt_ha: 0.0,
            chamfer_distance: 0.0,
            hausdorff_distance: 0.0,
            drmsd: 0.0,
            expected_node_count: 0,
        }
    }

    /// Require 80% node coverage and RMSD <= 1.0.
    pub fn default(node_count: usize) -> Self {
        MatchFilter {
            node_count: 0,
            node_ratio: 0.8,
            idf_per_match: 0.0,
            evalue: f64::MAX,
            rmsd: 1.0,
            tm_score: 0.0,
            gdt_ts: 0.0,
            gdt_ha: 0.0,
            chamfer_distance: 0.0,
            hausdorff_distance: 0.0,
            drmsd: 0.0,
            expected_node_count: node_count,
        }
    }


    /// True if `result` passes every enabled cutoff.
    #[inline]
    pub fn filter(&self, result: &MatchResult) -> bool {
        let mut pass = true;
        
        if self.node_count > 0 {
            pass = pass && result.node_count >= self.node_count;
        }
        if self.node_ratio > 0.0 {
            pass = pass && result.node_count as f32 / self.expected_node_count as f32 >= self.node_ratio;
        }
        if self.idf_per_match > 0.0 {
            pass = pass && result.idf >= self.idf_per_match;
        }
        if self.evalue < f64::MAX {
            pass = pass && result.evalue <= self.evalue;
        }
        if self.rmsd > 0.0 {
            pass = pass && result.rmsd <= self.rmsd;
        }
        
        // Higher is better
        if self.tm_score > 0.0 {
            pass = pass && result.metrics.tm_score >= self.tm_score;
        }
        if self.gdt_ts > 0.0 {
            pass = pass && result.metrics.gdt_ts >= self.gdt_ts;
        }
        if self.gdt_ha > 0.0 {
            pass = pass && result.metrics.gdt_ha >= self.gdt_ha;
        }
        
        // Lower is better
        if self.chamfer_distance > 0.0 {
            pass = pass && result.metrics.chamfer_distance <= self.chamfer_distance;
        }
        if self.hausdorff_distance > 0.0 {
            pass = pass && result.metrics.hausdorff_distance <= self.hausdorff_distance;
        }
        if self.drmsd > 0.0 {
            pass = pass && result.metrics.drmsd <= self.drmsd;
        }

        pass
    }
}

// TODO: Need testing