// File: sort.rs
// Created: 2025-11-09 15:23:25
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2025 Hyunbin Kim, All rights reserved

//! Sorting strategies for motif match results

use std::cmp::Ordering;

use super::result::{MatchResult, StructureResult};

/// Individual sorting criterion for motif match results
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SortKey {
    /// Node count (number of matched residues)
    NodeCount,
    /// IDF score
    Idf,
    /// E-value (from IDF)
    Evalue,
    /// RMSD
    Rmsd,
    /// TM-score
    TmScore,
    /// GDT-TS (Global Distance Test - Total Score)
    GdtTs,
    /// GDT-HA (Global Distance Test - High Accuracy)
    GdtHa,
    /// Chamfer distance
    ChamferDistance,
    /// Hausdorff distance
    HausdorffDistance,
}

impl SortKey {
    /// Parse from string (case-insensitive)
    /// 
    /// # Valid key names
    /// - `node_count`, `nodes`, `n` -> NodeCount
    /// - `idf` -> Idf
    /// - `rmsd` -> Rmsd
    /// - `tm_score`, `tmscore`, `tm` -> TmScore
    /// - `tm_score_strict`, `tmscore_strict`, `tm_strict` -> TmScoreStrict
    /// - `gdt_ts`, `gdtts`, `gdt` -> GdtTs
    /// - `gdt_ha`, `gdtha` -> GdtHa
    /// - `gdt_strict`, `gdtstrict` -> GdtStrict
    /// - `chamfer`, `chamfer_distance` -> ChamferDistance
    /// - `hausdorff`, `hausdorff_distance` -> HausdorffDistance
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s.trim().to_lowercase().as_str() {
            "node_count" | "node-count" | "nodes" | "node" | "n" => Ok(Self::NodeCount),
            "idf" | "score" => Ok(Self::Idf),
            "evalue" | "e_value" | "e-value" => Ok(Self::Evalue),
            "rmsd" => Ok(Self::Rmsd),
            "tm_score" | "tm-score" | "tmscore" | "tm" => Ok(Self::TmScore),
            "gdt_ts" | "gdt-ts" | "gdtts" | "gdt" => Ok(Self::GdtTs),
            "gdt_ha" | "gdt-ha" | "gdtha" => Ok(Self::GdtHa),
            "chamfer" | "chamfer-distance" | "chamfer_distance" => Ok(Self::ChamferDistance),
            "hausdorff" | "hausdorff-distance" | "hausdorff_distance" => Ok(Self::HausdorffDistance),
            _ => Err(format!(
                "Unknown sort key: '{}'. Valid keys: {}",
                s,
                Self::valid_keys()
            )),
        }
    }

    /// Get all valid key names for help text
    pub fn valid_keys() -> &'static str {
        "node_count, idf, evalue, rmsd, tm_score, tm_score_strict, gdt_ts, gdt_ha, gdt_strict, chamfer_distance, hausdorff_distance"
    }

    /// Get the default sort order for this key
    /// 
    /// Higher is better: NodeCount, IDF, TM-score, GDT scores -> Descending
    /// Lower is better: RMSD, Chamfer, Hausdorff -> Ascending
    pub fn default_order(&self) -> SortOrder {
        match self {
            // Descending order for NodeCount, IDF, TM-score, GDT scores
            Self::NodeCount | Self::Idf | Self::TmScore | Self::GdtTs | Self::GdtHa => SortOrder::Desc,
            // Ascending order for distance metrics: RMSD, Chamfer, Hausdorff, E-value
            Self::Evalue | Self::Rmsd | Self::ChamferDistance | Self::HausdorffDistance => SortOrder::Asc,
        }
    }

    /// Extract value from MatchResult
    fn extract_value(&self, result: &MatchResult) -> f64 {
        match self {
            Self::NodeCount => result.node_count as f64,
            Self::Idf => result.idf as f64,
            Self::Evalue => result.evalue,
            Self::Rmsd => result.rmsd as f64,
            Self::TmScore => result.metrics.tm_score as f64,
            Self::GdtTs => result.metrics.gdt_ts as f64,
            Self::GdtHa => result.metrics.gdt_ha as f64,
            Self::ChamferDistance => result.metrics.chamfer_distance as f64,
            Self::HausdorffDistance => result.metrics.hausdorff_distance as f64,
        }
    }

    /// Compare two MatchResults by this key
    pub fn compare(&self, a: &MatchResult, b: &MatchResult, order: SortOrder) -> Ordering {
        let val_a = self.extract_value(a);
        let val_b = self.extract_value(b);
        // Less comes first, Greater comes last
        // a partial cmp b means: a < b -> Less, a == b -> Equal, a > b -> Greater. 
        match order {
            SortOrder::Asc => val_a.partial_cmp(&val_b).unwrap_or(Ordering::Equal),
            SortOrder::Desc => val_b.partial_cmp(&val_a).unwrap_or(Ordering::Equal),
        }
    }
}

/// Sort direction
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SortOrder {
    /// Ascending (smaller is better)
    Asc,
    /// Descending (larger is better)
    Desc,
}

impl SortOrder {
    /// Parse from string (case-insensitive)
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s.trim().to_lowercase().as_str() {
            "asc" | "ascending" | "a" => Ok(Self::Asc),
            "desc" | "descending" | "d" => Ok(Self::Desc),
            _ => Err(format!(
                "Unknown sort order: '{}'. Use 'asc' or 'desc'",
                s
            )),
        }
    }
}

/// Multi-level sorting strategy for motif match results
#[derive(Debug, Clone)]
pub struct MatchSortStrategy {
    keys: Vec<(SortKey, SortOrder)>,
}

impl MatchSortStrategy {
    /// Create a new empty sorting strategy
    pub fn new() -> Self {
        Self { keys: Vec::new() }
    }

    /// Add a sorting key with custom order
    pub fn then_by(mut self, key: SortKey, order: SortOrder) -> Self {
        self.keys.push((key, order));
        self
    }

    /// Add a sorting key with default order
    pub fn then_by_default(mut self, key: SortKey) -> Self {
        self.keys.push((key, key.default_order()));
        self
    }

    /// Parse from comma-separated string
    ///
    /// Format options:
    /// 1. `"key1,key2,key3"` - Uses default order for each key
    /// 2. `"key1:order1,key2:order2"` - Custom order for each key
    /// 3. Mixed: `"key1,key2:desc,key3"` - Mix default and custom
    ///
    /// # Examples
    /// ```
    /// use folddisco::controller::sort::MatchSortStrategy;
    ///
    /// // Default orders
    /// let strategy = MatchSortStrategy::from_str("node_count,tm_score,rmsd").unwrap();
    ///
    /// // Custom orders
    /// let strategy = MatchSortStrategy::from_str("node_count:desc,tm_score:desc,rmsd:asc").unwrap();
    ///
    /// // Mixed
    /// let strategy = MatchSortStrategy::from_str("node_count,tm_score:asc,rmsd").unwrap();
    /// ```
    pub fn from_str(s: &str) -> Result<Self, String> {
        if s.trim().is_empty() {
            return Ok(Self::default());
        }

        let mut strategy = Self::new();

        for part in s.split(',') {
            let part = part.trim();
            if part.is_empty() {
                continue;
            }

            // Check if it contains order specification (key:order)
            if part.contains(':') {
                let components: Vec<&str> = part.split(':').collect();
                if components.len() != 2 {
                    return Err(format!(
                        "Invalid format: '{}'. Use 'key:order' or just 'key'",
                        part
                    ));
                }
                let key = SortKey::from_str(components[0])?;
                let order = SortOrder::from_str(components[1])?;
                strategy = strategy.then_by(key, order);
            } else {
                // Just key, use default order
                let key = SortKey::from_str(part)?;
                strategy = strategy.then_by_default(key);
            }
        }

        if strategy.keys.is_empty() {
            Ok(Self::default())
        } else {
            Ok(strategy)
        }
    }

    /// Compare two MatchResults using all keys in order
    pub fn compare(&self, a: &MatchResult, b: &MatchResult) -> Ordering {
        for (key, order) in &self.keys {
            match key.compare(a, b, *order) {
                Ordering::Equal => continue,
                other => return other,
            }
        }
        Ordering::Equal
    }

    /// Default strategy: NodeCount (desc) -> RMSD (asc)
    /// This matches the legacy behavior
    pub fn default() -> Self {
        Self::new()
            .then_by_default(SortKey::NodeCount)
            .then_by_default(SortKey::Rmsd)
    }

    /// NodeCount (desc) -> IDF (desc)
    /// Legacy behavior when sort_by_score is true
    pub fn by_node_count_idf() -> Self {
        Self::new()
            .then_by_default(SortKey::NodeCount)
            .then_by_default(SortKey::Idf)
    }
    
    /// IDF (desc)
    /// Pure IDF sorting
    pub fn by_idf() -> Self {
        Self::new()
            .then_by_default(SortKey::Idf)
    }

    /// NodeCount (desc) -> TM-score (desc) -> RMSD (asc)
    /// Recommended for motif matching quality
    pub fn by_node_count_tm_score() -> Self {
        Self::new()
            .then_by_default(SortKey::NodeCount)
            .then_by_default(SortKey::TmScore)
            .then_by_default(SortKey::Rmsd)
    }

    /// NodeCount (desc) -> GDT-TS (desc) -> RMSD (asc)
    /// Alternative quality metric
    pub fn by_node_count_gdt_ts() -> Self {
        Self::new()
            .then_by_default(SortKey::NodeCount)
            .then_by_default(SortKey::GdtTs)
            .then_by_default(SortKey::Rmsd)
    }
}

impl Default for MatchSortStrategy {
    fn default() -> Self {
        Self::default()
    }
}

// Display implementation for MatchSortStrategy
impl std::fmt::Display for MatchSortStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let parts: Vec<String> = self.keys.iter().map(|(key, order)| {
            let order_str = match order {
                SortOrder::Asc => "asc",
                SortOrder::Desc => "desc",
            };
            format!("{}:{}", format!("{:?}", key).to_lowercase(), order_str)
        }).collect();
        write!(f, "{}", parts.join(", "))
    }
}


// ============================================================================
// Structure-level sorting (for StructureResult)
// ============================================================================

/// Individual sorting criterion for structure-level results
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StructureSortKey {
    /// Maximum matching node count
    MaxNodeCount,
    /// Node count (nodes with hashes)
    NodeCount,
    /// IDF score
    Idf,
    /// Minimum RMSD with max match
    MinRmsd,
    /// Total match count
    TotalMatchCount,
    /// Edge count
    EdgeCount,
    /// Number of residues
    Nres,
    /// pLDDT score
    Plddt,
}

impl StructureSortKey {
    /// Parse from string (case-insensitive)
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s.trim().to_lowercase().as_str() {
            "max_node_count" | "max-node-count" | "max_node" | "max-node" | "max_nodes" | "max-nodes" => Ok(Self::MaxNodeCount),
            "node_count" | "node-count" | "nodes" | "node" | "n" => Ok(Self::NodeCount),
            "idf" | "score" => Ok(Self::Idf),
            "min_rmsd" | "min-rmsd" | "rmsd" => Ok(Self::MinRmsd),
            "total_match_count" | "total-match-count" | "total_match" | "total-match" | "matches" | "match" => Ok(Self::TotalMatchCount),
            "edge_count" | "edge-count" | "edges" | "edge" | "e" => Ok(Self::EdgeCount),
            "nres" | "num_residues" | "num-residues" | "length" | "residues" | "residue" | "l" => Ok(Self::Nres),
            "plddt" => Ok(Self::Plddt),
            _ => Err(format!(
                "Unknown structure sort key: '{}'. Valid keys: {}",
                s,
                Self::valid_keys()
            )),
        }
    }

    /// Get all valid key names for help text
    pub fn valid_keys() -> &'static str {
        "max_node_count, node_count, idf, min_rmsd, total_match_count, edge_count, nres, plddt"
    }

    /// Get the default sort order for this key
    pub fn default_order(&self) -> SortOrder {
        match self {
            // Higher is better
            Self::MaxNodeCount | Self::NodeCount | Self::Idf | Self::TotalMatchCount | 
            Self::EdgeCount | Self::Nres | Self::Plddt => SortOrder::Desc,
            // Lower is better
            Self::MinRmsd => SortOrder::Asc,
        }
    }

    /// Extract value from StructureResult
    fn extract_value(&self, result: &StructureResult) -> f32 {
        match self {
            Self::MaxNodeCount => result.max_matching_node_count as f32,
            Self::NodeCount => result.node_count as f32,
            Self::Idf => result.idf,
            Self::MinRmsd => result.min_rmsd_with_max_match,
            Self::TotalMatchCount => result.total_match_count as f32,
            Self::EdgeCount => result.edge_count as f32,
            Self::Nres => result.nres as f32,
            Self::Plddt => result.plddt,
        }
    }

    /// Compare two StructureResults by this key
    pub fn compare(&self, a: &StructureResult, b: &StructureResult, order: SortOrder) -> Ordering {
        let val_a = self.extract_value(a);
        let val_b = self.extract_value(b);
        match order {
            SortOrder::Asc => val_a.partial_cmp(&val_b).unwrap_or(Ordering::Equal),
            SortOrder::Desc => val_b.partial_cmp(&val_a).unwrap_or(Ordering::Equal),
        }
    }
}

/// Multi-level sorting strategy for structure-level results
#[derive(Debug, Clone)]
pub struct StructureSortStrategy {
    keys: Vec<(StructureSortKey, SortOrder)>,
}

impl StructureSortStrategy {
    /// Create a new empty sorting strategy
    pub fn new() -> Self {
        Self { keys: Vec::new() }
    }

    /// Add a sorting key with custom order
    pub fn then_by(mut self, key: StructureSortKey, order: SortOrder) -> Self {
        self.keys.push((key, order));
        self
    }

    /// Add a sorting key with default order
    pub fn then_by_default(mut self, key: StructureSortKey) -> Self {
        self.keys.push((key, key.default_order()));
        self
    }

    /// Parse from comma-separated string
    ///
    /// Format options:
    /// 1. `"key1,key2,key3"` - Uses default order for each key
    /// 2. `"key1:order1,key2:order2"` - Custom order for each key
    /// 3. Mixed: `"key1,key2:desc,key3"` - Mix default and custom
    ///
    /// # Examples
    /// ```
    /// use folddisco::controller::sort::StructureSortStrategy;
    ///
    /// // Default orders
    /// let strategy = StructureSortStrategy::from_str("max_node_count,min_rmsd").unwrap();
    ///
    /// // Custom orders
    /// let strategy = StructureSortStrategy::from_str("max_node_count:desc,idf:desc").unwrap();
    ///
    /// // Mixed
    /// let strategy = StructureSortStrategy::from_str("max_node_count,idf:desc,min_rmsd").unwrap();
    /// ```
    pub fn from_str(s: &str) -> Result<Self, String> {
        if s.trim().is_empty() {
            return Ok(Self::default());
        }

        let mut strategy = Self::new();

        for part in s.split(',') {
            let part = part.trim();
            if part.is_empty() {
                continue;
            }

            // Check if it contains order specification (key:order)
            if part.contains(':') {
                let components: Vec<&str> = part.split(':').collect();
                if components.len() != 2 {
                    return Err(format!(
                        "Invalid format: '{}'. Use 'key:order' or just 'key'",
                        part
                    ));
                }
                let key = StructureSortKey::from_str(components[0])?;
                let order = SortOrder::from_str(components[1])?;
                strategy = strategy.then_by(key, order);
            } else {
                // Just key, use default order
                let key = StructureSortKey::from_str(part)?;
                strategy = strategy.then_by_default(key);
            }
        }

        if strategy.keys.is_empty() {
            Ok(Self::default())
        } else {
            Ok(strategy)
        }
    }

    /// Compare two StructureResults using all keys in order
    pub fn compare(&self, a: &StructureResult, b: &StructureResult) -> Ordering {
        for (key, order) in &self.keys {
            match key.compare(a, b, *order) {
                Ordering::Equal => continue,
                other => return other,
            }
        }
        Ordering::Equal
    }

    /// Default strategy: MaxNodeCount (desc) -> MinRmsd (asc)
    /// This matches the legacy behavior
    pub fn default() -> Self {
        Self::new()
            .then_by_default(StructureSortKey::MaxNodeCount)
            .then_by_default(StructureSortKey::MinRmsd)
    }

    /// MaxNodeCount (desc) -> IDF (desc)
    pub fn by_max_node_count_idf() -> Self {
        Self::new()
            .then_by_default(StructureSortKey::MaxNodeCount)
            .then_by_default(StructureSortKey::Idf)
    }

    /// IDF (desc) only
    pub fn by_idf() -> Self {
        Self::new()
            .then_by_default(StructureSortKey::Idf)
    }
}

impl Default for StructureSortStrategy {
    fn default() -> Self {
        Self::default()
    }
}

// Implementation for printing with formatter

impl std::fmt::Display for StructureSortStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let parts: Vec<String> = self.keys.iter().map(|(key, order)| {
            let order_str = match order {
                SortOrder::Asc => "asc",
                SortOrder::Desc => "desc",
            };
            format!("{}:{}", format!("{:?}", key).to_lowercase(), order_str)
        }).collect();
        write!(f, "{}", parts.join(", "))
    }
}


#[cfg(test)]
mod tests {
    use crate::structure::metrics::StructureSimilarityMetrics;

    use super::*;

    #[test]
    fn test_sort_key_from_str() {
        assert_eq!(SortKey::from_str("node_count").unwrap(), SortKey::NodeCount);
        assert_eq!(SortKey::from_str("nodes").unwrap(), SortKey::NodeCount);
        assert_eq!(SortKey::from_str("n").unwrap(), SortKey::NodeCount);
        assert_eq!(SortKey::from_str("tm_score").unwrap(), SortKey::TmScore);
        assert_eq!(SortKey::from_str("TM_SCORE").unwrap(), SortKey::TmScore);
        assert!(SortKey::from_str("invalid").is_err());
    }

    #[test]
    fn test_sort_order_from_str() {
        assert_eq!(SortOrder::from_str("asc").unwrap(), SortOrder::Asc);
        assert_eq!(SortOrder::from_str("ASC").unwrap(), SortOrder::Asc);
        assert_eq!(SortOrder::from_str("desc").unwrap(), SortOrder::Desc);
        assert!(SortOrder::from_str("invalid").is_err());
    }

    #[test]
    fn test_default_order() {
        assert_eq!(SortKey::NodeCount.default_order(), SortOrder::Desc);
        assert_eq!(SortKey::Rmsd.default_order(), SortOrder::Asc);
        assert_eq!(SortKey::TmScore.default_order(), SortOrder::Desc);
    }

    #[test]
    fn test_strategy_from_str() {
        // Simple
        let strategy = MatchSortStrategy::from_str("node_count").unwrap();
        assert_eq!(strategy.keys.len(), 1);

        // Multiple keys
        let strategy = MatchSortStrategy::from_str("node_count,tm_score,rmsd").unwrap();
        assert_eq!(strategy.keys.len(), 3);

        // With order
        let strategy = MatchSortStrategy::from_str("node_count:desc,tm_score:asc").unwrap();
        assert_eq!(strategy.keys.len(), 2);
        assert_eq!(strategy.keys[1].1, SortOrder::Asc);

        // Mixed
        let strategy = MatchSortStrategy::from_str("node_count,tm_score:asc,rmsd").unwrap();
        assert_eq!(strategy.keys.len(), 3);

        // Empty defaults
        let strategy = MatchSortStrategy::from_str("").unwrap();
        assert_eq!(strategy.keys.len(), 2); // Default has 2 keys
    }
    

    #[test]
    fn test_strategy_with_dummy_data() {
        let result_a = MatchResult {
            tid: "toyA",
            nid: 1,
            db_key: 1,
            node_count: 10,
            idf: 5.0,
            matching_residues: vec![],
            rmsd: 0.1,
            evalue: 0.00001,
            u_matrix: [[0.0; 3]; 3],
            t_matrix: [0.0; 3], 
            matching_coordinates: vec![],
            metrics: StructureSimilarityMetrics {
                tm_score: 0.9,
                gdt_ts: 0.8,
                gdt_ha: 0.4,
                chamfer_distance: 1.2,
                hausdorff_distance: 2.5,
            },
        };
        let result_b = MatchResult {
            tid: "toyB",
            nid: 2,
            db_key: 2,
            node_count: 8,
            idf: 6.0,
            matching_residues: vec![],
            rmsd: 0.24,
            evalue: 0.01,
            u_matrix: [[0.0; 3]; 3],
            t_matrix: [0.0; 3], 
            matching_coordinates: vec![],
            metrics: StructureSimilarityMetrics {
                tm_score: 0.75,
                gdt_ts: 0.9,
                gdt_ha: 0.3,
                chamfer_distance: 2.9,
                hausdorff_distance: 3.0,
            },
        };
        let result_c = MatchResult {
            tid: "toyC",
            nid: 3,
            db_key: 3,
            node_count: 7,
            idf: 7.0,
            matching_residues: vec![],
            rmsd: 0.5,
            evalue: 2.0,
            u_matrix: [[0.0; 3]; 3],
            t_matrix: [0.0; 3], 
            matching_coordinates: vec![],
            metrics: StructureSimilarityMetrics {
                tm_score: 0.96,
                gdt_ts: 0.3,
                gdt_ha: 0.1,
                chamfer_distance: 1.2,
                hausdorff_distance: 8.0,
            },
        };
        
        // NodeCount (desc) -> IDF (desc)
        let strategy = MatchSortStrategy::from_str("node_count:desc,idf:desc").unwrap();
        assert_eq!(strategy.compare(&result_a, &result_b), Ordering::Less); // 10,5.0 comes before 8,6.0
        assert_eq!(strategy.compare(&result_b, &result_c), Ordering::Less); // 8,6.0 comes before 7,7.0
        assert_eq!(strategy.compare(&result_c, &result_a), Ordering::Greater); // 7,7.0 comes after 10,5.0
        println!("Passed NodeCount -> IDF tests.");
        
        // IDF (desc)
        // A: 5.0 | B: 6.0 | C: 7.0
        let strategy = MatchSortStrategy::from_str("idf:desc").unwrap();
        assert_eq!(strategy.compare(&result_a, &result_b), Ordering::Greater); // 5.0 comes after 6.0
        assert_eq!(strategy.compare(&result_b, &result_c), Ordering::Greater); // 6.0 comes after 7.0
        assert_eq!(strategy.compare(&result_a, &result_c), Ordering::Greater); // 5.0 comes after 7.0
        println!("Passed IDF tests.");                
        
        // RMSD (asc)        
        let strategy = MatchSortStrategy::from_str("rmsd:asc").unwrap();
        assert_eq!(strategy.compare(&result_a, &result_b), Ordering::Less); // 0.1 comes before 0.24
        assert_eq!(strategy.compare(&result_b, &result_c), Ordering::Less); // 0.24 comes before 0.5
        assert_eq!(strategy.compare(&result_a, &result_c), Ordering::Less); // 0.1 comes before 0.5
        println!("Passed RMSD tests.");

        // TM-score (desc)
        let strategy = MatchSortStrategy::from_str("tm_score:desc").unwrap();
        assert_eq!(strategy.compare(&result_a, &result_b), Ordering::Less); // 0.9 comes before 0.75 
        assert_eq!(strategy.compare(&result_b, &result_c), Ordering::Greater); // 0.75 comes after 0.96
        assert_eq!(strategy.compare(&result_a, &result_c), Ordering::Greater); // 0.9 comes after 0.96
        println!("Passed TM-score tests.");
                
        // GDT-TS (desc)
        let strategy = MatchSortStrategy::from_str("gdt_ts:desc").unwrap();
        assert_eq!(strategy.compare(&result_a, &result_b), Ordering::Greater); // 0.8 comes after 0.9
        assert_eq!(strategy.compare(&result_b, &result_c), Ordering::Less); // 0.9 comes before 0.3
        assert_eq!(strategy.compare(&result_a, &result_c), Ordering::Less); // 0.8 comes before 0.3
        println!("Passed GDT-TS tests.");
        
        // Chamfer (asc) -> Hausdorff (asc)
        let strategy = MatchSortStrategy::from_str("chamfer_distance:asc,hausdorff_distance:asc").unwrap();
        assert_eq!(strategy.compare(&result_a, &result_b), Ordering::Less); // 1.2,2.5 comes before 2.9,3.0
        assert_eq!(strategy.compare(&result_b, &result_c), Ordering::Greater); // 2.9,3.0 comes after 1.2,8.0
        assert_eq!(strategy.compare(&result_a, &result_c), Ordering::Less); // 1.2,2.5 comes before 1.2,8.0
        println!("Passed Chamfer -> Hausdorff tests.");
    }

    #[test]
    fn test_structure_sort_key_from_str() {
        assert_eq!(StructureSortKey::from_str("max_node_count").unwrap(), StructureSortKey::MaxNodeCount);
        assert_eq!(StructureSortKey::from_str("max-node").unwrap(), StructureSortKey::MaxNodeCount);
        assert_eq!(StructureSortKey::from_str("node_count").unwrap(), StructureSortKey::NodeCount);
        assert_eq!(StructureSortKey::from_str("idf").unwrap(), StructureSortKey::Idf);
        assert_eq!(StructureSortKey::from_str("min_rmsd").unwrap(), StructureSortKey::MinRmsd);
        assert_eq!(StructureSortKey::from_str("rmsd").unwrap(), StructureSortKey::MinRmsd);
        assert!(StructureSortKey::from_str("invalid").is_err());
    }

    #[test]
    fn test_structure_default_order() {
        assert_eq!(StructureSortKey::MaxNodeCount.default_order(), SortOrder::Desc);
        assert_eq!(StructureSortKey::MinRmsd.default_order(), SortOrder::Asc);
        assert_eq!(StructureSortKey::Idf.default_order(), SortOrder::Desc);
    }

    #[test]
    fn test_structure_strategy_from_str() {
        // Simple
        let strategy = StructureSortStrategy::from_str("max_node_count").unwrap();
        assert_eq!(strategy.keys.len(), 1);

        // Multiple keys
        let strategy = StructureSortStrategy::from_str("max_node_count,idf,min_rmsd").unwrap();
        assert_eq!(strategy.keys.len(), 3);

        // With order
        let strategy = StructureSortStrategy::from_str("max_node_count:desc,idf:asc").unwrap();
        assert_eq!(strategy.keys.len(), 2);
        assert_eq!(strategy.keys[1].1, SortOrder::Asc);

        // Empty defaults
        let strategy = StructureSortStrategy::from_str("").unwrap();
        assert_eq!(strategy.keys.len(), 2); // Default has 2 keys
    }

    #[test]
    fn test_structure_strategy_with_dummy_data() {
        let result_a = StructureResult {
            tid: "structA",
            nid: 1,
            db_key: 1,
            total_match_count: 100,
            node_count: 50,
            edge_count: 45,
            idf: 5.0,
            nres: 200,
            plddt: 85.5,
            matching_residues: vec![],
            matching_residues_processed: vec![],
            max_matching_node_count: 10,
            min_rmsd_with_max_match: 0.5,
        };

        let result_b = StructureResult {
            tid: "structB",
            nid: 2,
            db_key: 2,
            total_match_count: 120,
            node_count: 60,
            edge_count: 55,
            idf: 7.0,
            nres: 250,
            plddt: 90.0,
            matching_residues: vec![],
            matching_residues_processed: vec![],
            max_matching_node_count: 8,
            min_rmsd_with_max_match: 0.3,
        };

        let result_c = StructureResult {
            tid: "structC",
            nid: 3,
            db_key: 3,
            total_match_count: 80,
            node_count: 40,
            edge_count: 35,
            idf: 6.0,
            nres: 180,
            plddt: 80.0,
            matching_residues: vec![],
            matching_residues_processed: vec![],
            max_matching_node_count: 10,
            min_rmsd_with_max_match: 0.7,
        };

        // MaxNodeCount (desc) -> MinRmsd (asc)
        let strategy = StructureSortStrategy::from_str("max_node_count:desc,min_rmsd:asc").unwrap();
        assert_eq!(strategy.compare(&result_a, &result_b), Ordering::Less); // 10,0.5 comes before 8,0.3
        assert_eq!(strategy.compare(&result_a, &result_c), Ordering::Less); // 10,0.5 comes before 10,0.7
        assert_eq!(strategy.compare(&result_b, &result_c), Ordering::Greater); // 8,0.3 comes after 10,0.7
        println!("Passed MaxNodeCount -> MinRmsd tests.");

        // IDF (desc)
        let strategy = StructureSortStrategy::from_str("idf:desc").unwrap();
        assert_eq!(strategy.compare(&result_a, &result_b), Ordering::Greater); // 5.0 comes after 7.0
        assert_eq!(strategy.compare(&result_b, &result_c), Ordering::Less); // 7.0 comes before 6.0
        assert_eq!(strategy.compare(&result_a, &result_c), Ordering::Greater); // 5.0 comes after 6.0
        println!("Passed IDF tests.");

        // NodeCount (desc) -> TotalMatchCount (desc)
        let strategy = StructureSortStrategy::from_str("node_count:desc,total_match_count:desc").unwrap();
        assert_eq!(strategy.compare(&result_a, &result_b), Ordering::Greater); // 50,100 comes after 60,120
        assert_eq!(strategy.compare(&result_a, &result_c), Ordering::Less); // 50,100 comes before 40,80
        println!("Passed NodeCount -> TotalMatchCount tests.");
    }
  
    #[test]  
    fn test_printing_strategies() {
        let match_strategy = MatchSortStrategy::from_str("node_count:desc,tm_score:asc,rmsd").unwrap();
        println!("MatchSortStrategy: {}", match_strategy);

        let struct_strategy = StructureSortStrategy::from_str("max_node_count:desc,min_rmsd:asc").unwrap();
        println!("StructureSortStrategy: {}", struct_strategy);
    }
}
