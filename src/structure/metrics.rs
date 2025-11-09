// File: metrics.rs
// Created: 2025-10-29
// Description: Structure similarity metrics implementation
//   - TM-score: Template Modeling score
//   - GDT-TS: Global Distance Test - Total Score.
//   - GDT-HA: Global Distance Test - High Accuracy
//   - Chamfer Distance: Average of nearest neighbor distances
//   - Hausdorff Distance: Maximum of minimum distances. Outlier sensitive.
//   - RMSD: Root Mean Square Deviation. Outlier sensitive.
//
// Performance optimizations:
//   - Uses f64 for internal arithmetic to maintain precision
//   - PrecomputedDistances struct for efficient batch metric calculation
//   - Avoids redundant distance calculations across metrics
//
// Usage example:
// ```rust
// use crate::structure::metrics::*;
//
// // Method 1: Calculate all metrics (simple but slower)
// let metrics = StructureMetrics::calculate_all(&ref_coords, &model_coords);
// println!("TM-score: {:.4}", metrics.tm_score);
//
// // Method 2: Calculate all metrics fast (2x faster, recommended)
// let metrics = StructureMetrics::calculate_all_fast(&ref_coords, &model_coords);
// metrics.print();
//
// // Method 3: Pre-compute distances and calculate individual metrics
// let distances = PrecomputedDistances::new(&ref_coords, &model_coords);
// let tm = tm_score_fast(&distances, None);
// let gdt = gdt_ts_fast(&distances);
// let lddt = lddt_default_fast(&distances);
// let rmsd = rmsd_fast(&distances);
// ```

use core::fmt;

/// Pre-computed distances for efficient metric calculation
/// 
/// This struct stores pre-calculated distances to avoid redundant computations
/// when calculating multiple metrics on the same structure pair.
/// 
/// Memory usage: ~4MB for 1000-residue protein
/// Speedup: 2-3x faster when calculating all metrics together
pub struct PrecomputedDistances {
    /// Squared distances between corresponding atoms: dist_sq(ref[i], model[i])
    /// Used by: TM-score, GDT, RMSD
    pub pairwise_dist: Vec<f32>,
    /// Number of atoms
    pub n: usize,
}

impl PrecomputedDistances {
    /// Pre-calculate all distances between two structures
    /// 
    /// # Arguments
    /// * `reference_coords` - Reference structure coordinates
    /// * `coords` - Model structure coordinates
    /// 
    /// # Returns
    /// PrecomputedDistances struct containing all distance calculations
    pub fn new(reference_coords: &[[f32; 3]], coords: &[[f32; 3]]) -> Self {
        let n = reference_coords.len();
        // Currently only supports equal-length structures. If empty or different lengths, return empty.
        if n == 0 || n != coords.len() {
            return Self {
                pairwise_dist: Vec::new(),
                n: 0,
            };
        }

        // Pre-calculate pairwise squared distances (for TM/GDT/RMSD)
        let mut pairwise_dist: Vec<f32> = Vec::with_capacity(n * n);

        for c in coords {
            for r in reference_coords {
                pairwise_dist.push(dist(*r, *c)); // This is after sqrt
            }
        }
        // Pairwise distances can be retrieved with pairwise_dist[i * n + j] 
        // instead of distances[i][j]
        
        Self {
            pairwise_dist,
            n,
        }
    }

    #[inline(always)]
    pub fn get_distance(&self, i: usize, j: usize) -> f32 {
        self.pairwise_dist[i * self.n + j]
    }

    
}

/// Calculate squared Euclidean distance between two 3D points.
/// This function uses f64 for intermediate calculations to improve precision.
#[inline(always)]
fn dist_sq_as_f64(a: [f32; 3], b: [f32; 3]) -> f64 {
    let dx = a[0] as f64 - b[0] as f64;
    let dy = a[1] as f64 - b[1] as f64;
    let dz = a[2] as f64 - b[2] as f64;
    (dx * dx + dy * dy + dz * dz) as f64
}

/// Calculate Euclidean distance between two 3D points
#[inline(always)]
fn dist(a: [f32; 3], b: [f32; 3]) -> f32 {
    dist_sq_as_f64(a, b).sqrt() as f32
}

/// TM-score normalization function for final evaluation
/// d0(L) = 1.24 * (L - 15)^(1/3) - 1.8 for L > 19
/// d0(L) = 0.5 for L <= 21
#[inline]
fn d0_scale(length: usize) -> f32 {
    if length > 21 {
        1.24 * ((length as f32 - 15.0).powf(1.0 / 3.0)) - 1.8
    } else {
        0.5
    }
}

/// Fast TM-score using precomputed distances
/// 
/// # Arguments
/// * `distances` - Precomputed distance data
/// * `d0` - Optional normalization parameter
/// 
/// # Returns
/// TM-score in range [0, 1]
pub fn tm_score(distances: &PrecomputedDistances, d0: Option<f32>) -> f32 {
    if distances.n == 0 {
        return 0.0;
    }
    
    let d0 = d0.unwrap_or_else(|| d0_scale(distances.n));
    let d0_sq = (d0 * d0) as f64;
    
    let sum: f64 = (0..distances.n)
        .map(|i| {
            let d_sq = distances.get_distance(i, i) as f64;
            1.0 / (1.0 + d_sq / d0_sq)
        })
        .sum();
    
    
    (sum / distances.n as f64) as f32
}

/// Fast GDT using precomputed distances
fn gdt_generic(distances: &PrecomputedDistances, cutoffs: &[f64]) -> f32 {
    if distances.n == 0 || cutoffs.is_empty() {
        return 0.0;
    }
    
    let mut sum = 0.0_f64;
    
    for &cutoff in cutoffs {
        let cutoff_sq = cutoff * cutoff;
        let count = (0..distances.n)
            .filter(|&i| {
                let d_sq = distances.get_distance(i, i) as f64;
                d_sq <= cutoff_sq
            })
            .count();
        
        sum += (count as f64) / (distances.n as f64);
    }
    
    (sum / cutoffs.len() as f64) as f32
}

/// Fast GDT-TS using precomputed distances
pub fn gdt_ts(distances: &PrecomputedDistances) -> f32 {
    const CUTOFFS: [f64; 4] = [1.0, 2.0, 4.0, 8.0];
    gdt_generic(distances, &CUTOFFS)
}

/// Fast GDT-HA using precomputed distances
pub fn gdt_ha(distances: &PrecomputedDistances) -> f32 {
    const CUTOFFS: [f64; 4] = [0.5, 1.0, 2.0, 4.0];
    gdt_generic(distances, &CUTOFFS)
}

/// Fast GDT-strict using precomputed distances
pub fn gdt_strict(distances: &PrecomputedDistances) -> f32 {
    const CUTOFFS: [f64; 4] = [0.25, 0.5, 1.0, 2.0];
    gdt_generic(distances, &CUTOFFS)
}


/// Calculate Chamfer Distance between two point sets
/// 
/// Chamfer Distance is the mean of:
/// - Average nearest neighbor distance from coords to reference_coords
///
/// # Arguments
/// * `distance` - Precomputed distance data
///
/// # Returns
/// Chamfer distance (lower is better, 0 is perfect match)
pub fn chamfer_distance(distance: &PrecomputedDistances) -> f32 {
    if distance.n == 0 {
        return f32::INFINITY;
    }

    // Use f64 for accumulation to maintain precision
    // Average min distance from coords to reference
    let sum_coords_to_ref: f64 = (0..distance.n)
        .map(|i| {
            (0..distance.n)
                .map(|j| distance.get_distance(i, j) as f64)
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap()
        })
        .sum();
    
    (sum_coords_to_ref / distance.n as f64) as f32
}


/// Calculate Hausdorff Distance between two point sets
/// 
/// Hausdorff Distance is the maximum of:
/// - Maximum nearest neighbor distance from coords to reference_coords
///
/// # Arguments
/// * `distance` - Precomputed distance data
///
/// # Returns
/// Hausdorff distance (lower is better, 0 is perfect match)
pub fn hausdorff_distance(distance: &PrecomputedDistances) -> f32 {
    if distance.n == 0 {
        return f32::INFINITY;
    }

    // Max min distance from coords to reference
    let max_coords_to_ref = (0..distance.n)
        .map(|i| {
            (0..distance.n)
                .map(|j| distance.get_distance(i, j))
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap()
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();

    max_coords_to_ref
}


/// Calculate RMSD (Root Mean Square Deviation) between aligned structures
/// 
/// # Arguments
/// * `distance` - Precomputed distance data
///
/// # Returns
/// RMSD value in Ångströms
/// 
pub fn rmsd(distances: &PrecomputedDistances) -> f32 {
    if distances.n == 0 {
        return 0.0;
    }
    
    // Use f64 for accumulation to maintain precision
    let sum_sq: f64 = (0..distances.n)
        .map(|i| (distances.get_distance(i, i) as f64).powi(2))
        .sum();

    ((sum_sq / distances.n as f64).sqrt()) as f32
}


/// Structure similarity metrics calculator
#[derive(Debug, Clone, Default, PartialEq, Copy)]
pub struct StructureSimilarityMetrics {
    pub tm_score: f32,
    pub tm_score_strict: f32,
    pub gdt_ts: f32,
    pub gdt_ha: f32,
    pub gdt_strict: f32,
    pub chamfer_distance: f32,
    pub hausdorff_distance: f32,
    pub rmsd: f32,
}

impl StructureSimilarityMetrics {
    
    /// Calculate all metrics efficiently using precomputed distances
    /// 
    /// This is 2-3x faster than calculate_all() because distances are computed only once.
    /// Recommended for batch processing or when calculating multiple metrics.
    /// 
    /// # Arguments
    /// * `reference_coords` - Reference structure coordinates
    /// * `coords` - Model structure coordinates (should be pre-aligned)
    ///
    /// # Returns
    /// StructureMetrics containing all calculated metrics
    pub fn new() -> Self {      
        Self {
            tm_score: 0.0,
            tm_score_strict: 0.0,
            gdt_ts: 0.0,
            gdt_ha: 0.0,
            gdt_strict: 0.0,
            chamfer_distance: 0.0,
            hausdorff_distance: 0.0,
            rmsd: 0.0,
        }
    }

    pub fn calculate_tm_score(&self, precomputed: &PrecomputedDistances) -> f32 {
        tm_score(precomputed, None)
    }

    pub fn calculate_tm_score_strict(&self, precomputed: &PrecomputedDistances) -> f32 {
        // Using d0 = 0.168 which is used for searching in TM-score source
        tm_score(precomputed, Some(0.168))
    }

    pub fn calculate_gdt_ts(&self, precomputed: &PrecomputedDistances) -> f32 {
        gdt_ts(precomputed)
    }

    pub fn calculate_gdt_ha(&self, precomputed: &PrecomputedDistances) -> f32 {
        gdt_ha(precomputed)
    }
    
    pub fn calculate_gdt_strict(&self, precomputed: &PrecomputedDistances) -> f32 {
        gdt_strict(precomputed)
    }

    pub fn calculate_chamfer_distance(&self, precomputed: &PrecomputedDistances) -> f32 {
        chamfer_distance(precomputed)
    }

    pub fn calculate_hausdorff_distance(&self, precomputed: &PrecomputedDistances) -> f32 {
        hausdorff_distance(precomputed)
    }

    pub fn calculate_rmsd(&self, precomputed: &PrecomputedDistances) -> f32 {
        rmsd(precomputed)
    }

    pub fn calculate_all(&mut self, precomputed: &PrecomputedDistances) {
        self.tm_score = self.calculate_tm_score(precomputed);
        self.tm_score_strict = self.calculate_tm_score_strict(precomputed);
        self.gdt_ts = self.calculate_gdt_ts(precomputed);
        self.gdt_ha = self.calculate_gdt_ha(precomputed);
        self.gdt_strict = self.calculate_gdt_strict(precomputed);
        self.chamfer_distance = self.calculate_chamfer_distance(precomputed);
        self.hausdorff_distance = self.calculate_hausdorff_distance(precomputed);
        self.rmsd = self.calculate_rmsd(precomputed);
    }
    
    /// Print metrics in a formatted way
    pub fn print_in_a_formatted_way(&self) {
        println!("Structure Similarity Metrics:");
        println!("  TM-score:           {:.4}", self.tm_score);
        println!("  TM-score (strict):  {:.4}", self.tm_score_strict);
        println!("  GDT-TS:             {:.4}", self.gdt_ts);
        println!("  GDT-HA:             {:.4}", self.gdt_ha);
        println!("  GDT-strict:         {:.4}", self.gdt_strict);
        println!("  RMSD:               {:.4} Å", self.rmsd);
        println!("  Chamfer Distance:   {:.4} Å", self.chamfer_distance);
        println!("  Hausdorff Distance: {:.4} Å", self.hausdorff_distance);
    }
}

impl fmt::Display for StructureSimilarityMetrics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Print all metrics in a tab-separated format with 4 decimal places
        write!(
            f,
            "{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            self.tm_score,
            self.tm_score_strict,
            self.gdt_ts,
            self.gdt_ha,
            self.gdt_strict,
            self.rmsd,
            self.chamfer_distance,
            self.hausdorff_distance
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::structure::{kabsch::KabschSuperimposer, lms_qcp::LmsQcpSuperimposer};

    use super::*;

    #[test]
    fn test_metrics_calculate_all_with_identical() {
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let precomputed = PrecomputedDistances::new(&coords, &coords);
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_all(&precomputed);

        assert!((metrics.tm_score - 1.0).abs() < 1e-6);
        assert!((metrics.tm_score_strict - 1.0).abs() < 1e-6);
        assert!((metrics.gdt_ts - 1.0).abs() < 1e-6);
        assert!((metrics.gdt_ha - 1.0).abs() < 1e-6);
        assert!((metrics.gdt_strict - 1.0).abs() < 1e-6);
        assert!(metrics.chamfer_distance.abs() < 1e-6);
        assert!(metrics.hausdorff_distance.abs() < 1e-6);
        assert!(metrics.rmsd.abs() < 1e-6);
        
        metrics.print_in_a_formatted_way();
    }
    
    
    #[test]
    fn test_with_real_coordinates() {
        // Read coordinates from PDB files
        use crate::prelude::PDBReader;
        let query_reader = PDBReader::from_file("query/1G2F.pdb").unwrap();
        let query_zinc_structure = query_reader.read_structure().unwrap().to_compact();
        let target_reader = PDBReader::from_file("data/zinc/AF-P36508-F1-model_v6.pdb").unwrap();
        let target_zinc_structure = target_reader.read_structure().unwrap().to_compact();
        // Get reference coordinates: F207,F212,F225,F229
        let reference_indices = vec![
            query_zinc_structure.get_index(&b'F', &207).unwrap(),
            query_zinc_structure.get_index(&b'F', &212).unwrap(),
            query_zinc_structure.get_index(&b'F', &225).unwrap(),
            query_zinc_structure.get_index(&b'F', &229).unwrap(),
        ];
        println!("Reference indices: {:?}", reference_indices);
        let reference_coords = vec![
            query_zinc_structure.get_ca(reference_indices[0]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[0]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[1]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[1]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[2]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[2]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[3]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[3]).unwrap(),
        ];
        println!("Reference coords: {:?}", reference_coords);
        
        // Get target coordinates: A257,A262,A275,A279
        let target_indices = vec![
            target_zinc_structure.get_index(&b'A', &257).unwrap(),
            target_zinc_structure.get_index(&b'A', &262).unwrap(),
            target_zinc_structure.get_index(&b'A', &275).unwrap(),
            target_zinc_structure.get_index(&b'A', &279).unwrap(),
        ];
        println!("Target indices: {:?}", target_indices);
        let target_coords = vec![
            target_zinc_structure.get_ca(target_indices[0]).unwrap(),
            target_zinc_structure.get_cb(target_indices[0]).unwrap(),
            target_zinc_structure.get_ca(target_indices[1]).unwrap(),
            target_zinc_structure.get_cb(target_indices[1]).unwrap(),
            target_zinc_structure.get_ca(target_indices[2]).unwrap(),
            target_zinc_structure.get_cb(target_indices[2]).unwrap(),
            target_zinc_structure.get_ca(target_indices[3]).unwrap(),
            target_zinc_structure.get_cb(target_indices[3]).unwrap(),
        ];
        println!("Target coords: {:?}", target_coords);
        
        let mut kabsch = KabschSuperimposer::new();
        // let mut kabsch = LmsQcpSuperimposer::new();

        kabsch.set_atoms(&reference_coords, &target_coords);
        kabsch.run();

        let precomputed = PrecomputedDistances::new(
            &kabsch.reference_coords.unwrap(), &kabsch.transformed_coords.unwrap()
        );
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_all(&precomputed);

        metrics.print_in_a_formatted_way();
        
    }

    
    #[test]
    fn test_with_long_coordinates() {
        // Read coordinates from PDB files
        use crate::prelude::PDBReader;
        let query_reader = PDBReader::from_file("query/1G2F.pdb").unwrap();
        let query_zinc_structure = query_reader.read_structure().unwrap().to_compact();
        let target_reader = PDBReader::from_file("data/zinc/AF-P36508-F1-model_v6.pdb").unwrap();
        let target_zinc_structure = target_reader.read_structure().unwrap().to_compact();
        // Get reference coordinates: F205-214,F223-232
        let mut reference_indices = (207..213).map(|res_num| {
            query_zinc_structure.get_index(&b'F', &res_num).unwrap()
        }).collect::<Vec<usize>>();
        reference_indices.extend((225..230).map(|res_num| {
            query_zinc_structure.get_index(&b'F', &res_num).unwrap()
        }));
        println!("Reference indices: {:?}", reference_indices);
        let reference_coords = reference_indices.iter().flat_map(|&idx| {
            vec![
                query_zinc_structure.get_n(idx).unwrap(),
                query_zinc_structure.get_ca(idx).unwrap(),
                query_zinc_structure.get_cb(idx).unwrap(),
            ]
        }).collect::<Vec<_>>();
        println!("Reference coords: {:?}", reference_coords);
        
        // Get target coordinates: A255-260,A273-282
        let mut target_indices = (256..262).map(|res_num| {
            target_zinc_structure.get_index(&b'A', &res_num).unwrap()
        }).collect::<Vec<usize>>();
        target_indices.extend((275..280).map(|res_num| {
            target_zinc_structure.get_index(&b'A', &res_num).unwrap()
        }));
        println!("Target indices: {:?}", target_indices);
        let target_coords = target_indices.iter().flat_map(|&idx| {
            vec![
                target_zinc_structure.get_n(idx).unwrap(),
                target_zinc_structure.get_ca(idx).unwrap(),
                target_zinc_structure.get_cb(idx).unwrap(),
            ]
        }).collect::<Vec<_>>();
        println!("Target coords: {:?}", target_coords);
        
        // let mut kabsch = KabschSuperimposer::new();
        let mut kabsch = LmsQcpSuperimposer::new();

        kabsch.set_atoms(&reference_coords, &target_coords);
        kabsch.run();

        let precomputed = PrecomputedDistances::new(
            &kabsch.reference_coords.unwrap(), &kabsch.transformed_coords.unwrap()
        );
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_all(&precomputed);

        metrics.print_in_a_formatted_way();
        
    }
    
    #[test]
    fn test_with_outlier_coordinates() {
        // Read coordinates from PDB files
        use crate::prelude::PDBReader;
        let query_reader = PDBReader::from_file("query/1G2F.pdb").unwrap();
        let query_zinc_structure = query_reader.read_structure().unwrap().to_compact();
        let target_reader = PDBReader::from_file("data/zinc/AF-P36508-F1-model_v6.pdb").unwrap();
        let target_zinc_structure = target_reader.read_structure().unwrap().to_compact();
        // Get reference coordinates: F205,F212,F225,F229
        let reference_indices = vec![
            query_zinc_structure.get_index(&b'F', &205).unwrap(), // Outlier
            query_zinc_structure.get_index(&b'F', &212).unwrap(),
            query_zinc_structure.get_index(&b'F', &225).unwrap(),
            query_zinc_structure.get_index(&b'F', &229).unwrap(),
        ];
        println!("Reference indices: {:?}", reference_indices);
        let reference_coords = vec![
            query_zinc_structure.get_ca(reference_indices[0]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[0]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[1]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[1]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[2]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[2]).unwrap(),
            query_zinc_structure.get_ca(reference_indices[3]).unwrap(),
            query_zinc_structure.get_cb(reference_indices[3]).unwrap(),
        ];
        println!("Reference coords: {:?}", reference_coords);
        
        // Get target coordinates: A257,A262,A275,A279
        let target_indices = vec![
            target_zinc_structure.get_index(&b'A', &257).unwrap(),
            target_zinc_structure.get_index(&b'A', &262).unwrap(),
            target_zinc_structure.get_index(&b'A', &275).unwrap(),
            target_zinc_structure.get_index(&b'A', &279).unwrap(),
        ];
        println!("Target indices: {:?}", target_indices);
        let target_coords = vec![
            target_zinc_structure.get_ca(target_indices[0]).unwrap(),
            target_zinc_structure.get_cb(target_indices[0]).unwrap(),
            target_zinc_structure.get_ca(target_indices[1]).unwrap(),
            target_zinc_structure.get_cb(target_indices[1]).unwrap(),
            target_zinc_structure.get_ca(target_indices[2]).unwrap(),
            target_zinc_structure.get_cb(target_indices[2]).unwrap(),
            target_zinc_structure.get_ca(target_indices[3]).unwrap(),
            target_zinc_structure.get_cb(target_indices[3]).unwrap(),
        ];
        println!("Target coords: {:?}", target_coords);
        
        // let mut kabsch = KabschSuperimposer::new();
        let mut kabsch = LmsQcpSuperimposer::new();
        kabsch.set_atoms(&reference_coords, &target_coords);
        kabsch.run();

        let precomputed = PrecomputedDistances::new(
            &kabsch.reference_coords.unwrap(), &kabsch.transformed_coords.unwrap()
        );
        let mut metrics = StructureSimilarityMetrics::new();
        metrics.calculate_all(&precomputed);

        metrics.print_in_a_formatted_way();
        
    }

}
