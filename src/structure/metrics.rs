// File: metrics.rs
// Created: 2025-10-29
// Description: Structure similarity metrics (TM-score, GDT-TS/HA, Chamfer,
// Hausdorff, RMSD, dRMSD) over pre-superposed coordinates.

use core::fmt;

/// All model-to-reference distances of two equal-length point sets, computed once
/// and shared by every metric.
pub struct PrecomputedDistances {
    /// Row-major `n * n` Euclidean distances: `pairwise_dist[i * n + j] = |model_i - reference_j|`.
    pub pairwise_dist: Vec<f32>,
    /// Number of atoms
    pub n: usize,
}

impl PrecomputedDistances {
    /// Distances between two structures. Empty or unequal-length inputs give `n = 0`.
    pub fn new(reference_coords: &[[f32; 3]], coords: &[[f32; 3]]) -> Self {
        let n = reference_coords.len();
        if n == 0 || n != coords.len() {
            return Self {
                pairwise_dist: Vec::new(),
                n: 0,
            };
        }

        let mut pairwise_dist: Vec<f32> = Vec::with_capacity(n * n);

        for c in coords {
            for r in reference_coords {
                pairwise_dist.push(dist(*r, *c)); // This is after sqrt
            }
        }

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

/// Squared Euclidean distance, accumulated in f64 for precision.
#[inline(always)]
fn dist_sq_as_f64(a: [f32; 3], b: [f32; 3]) -> f64 {
    let dx = a[0] as f64 - b[0] as f64;
    let dy = a[1] as f64 - b[1] as f64;
    let dz = a[2] as f64 - b[2] as f64;
    (dx * dx + dy * dy + dz * dz) as f64
}

#[inline(always)]
fn dist(a: [f32; 3], b: [f32; 3]) -> f32 {
    dist_sq_as_f64(a, b).sqrt() as f32
}

/// TM-score d0: `1.24 * (L - 15)^(1/3) - 1.8` for L > 21, else 0.5.
#[inline]
fn d0_scale(length: usize) -> f32 {
    if length > 21 {
        1.24 * ((length as f32 - 15.0).powf(1.0 / 3.0)) - 1.8
    } else {
        0.5
    }
}

/// TM-score in [0, 1] over corresponding atoms; `d0` defaults to `d0_scale(n)`.
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

/// Mean fraction of corresponding atoms within each cutoff.
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

/// GDT-TS: cutoffs 1, 2, 4, 8 A.
pub fn gdt_ts(distances: &PrecomputedDistances) -> f32 {
    const CUTOFFS: [f64; 4] = [1.0, 2.0, 4.0, 8.0];
    gdt_generic(distances, &CUTOFFS)
}

/// GDT-HA: cutoffs 0.5, 1, 2, 4 A.
pub fn gdt_ha(distances: &PrecomputedDistances) -> f32 {
    const CUTOFFS: [f64; 4] = [0.5, 1.0, 2.0, 4.0];
    gdt_generic(distances, &CUTOFFS)
}

// pub fn gdt_strict(distances: &PrecomputedDistances) -> f32 {
//     const CUTOFFS: [f64; 4] = [0.25, 0.5, 1.0, 2.0];
//     gdt_generic(distances, &CUTOFFS)
// }


/// One-directional Chamfer distance: mean nearest-neighbour distance from model to
/// reference. 0 is a perfect match; empty input gives infinity.
pub fn chamfer_distance(distance: &PrecomputedDistances) -> f32 {
    if distance.n == 0 {
        return f32::INFINITY;
    }

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


/// One-directional Hausdorff distance: largest nearest-neighbour distance from model
/// to reference. 0 is a perfect match; empty input gives infinity.
pub fn hausdorff_distance(distance: &PrecomputedDistances) -> f32 {
    if distance.n == 0 {
        return f32::INFINITY;
    }

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


/// RMSD in Angstroms between already superposed structures.
pub fn rmsd(distances: &PrecomputedDistances) -> f32 {
    if distances.n == 0 {
        return 0.0;
    }

    let sum_sq: f64 = (0..distances.n)
        .map(|i| (distances.get_distance(i, i) as f64).powi(2))
        .sum();

    ((sum_sq / distances.n as f64).sqrt()) as f32
}


/// Superposition-free deformation: `(dRMSD, max |d_ref(i,j) - d_model(i,j)|)` over i<j,
/// where `dRMSD = sqrt(mean((|a_i - a_j| - |b_i - b_j|)^2))`.
///
/// Suited to non-rigid matches (e.g. hinged motifs) that superposition RMSD penalises.
/// Distances come through closures to avoid copying coordinates. `(0.0, 0.0)` for n < 2.
pub fn deformation_stats_indexed(
    n: usize,
    reference_distance: impl Fn(usize, usize) -> f32,
    model_distance: impl Fn(usize, usize) -> f32,
) -> (f32, f32) {
    if n < 2 {
        return (0.0, 0.0);
    }
    let mut sum_sq = 0.0_f64;
    let mut worst = 0.0_f32;
    let mut count = 0usize;
    for i in 0..n - 1 {
        for j in i + 1..n {
            let deviation = reference_distance(i, j) - model_distance(i, j);
            sum_sq += (deviation as f64) * (deviation as f64);
            let magnitude = deviation.abs();
            if magnitude > worst {
                worst = magnitude;
            }
            count += 1;
        }
    }
    if count == 0 {
        return (0.0, 0.0);
    }
    (((sum_sq / count as f64).sqrt()) as f32, worst)
}

/// Similarity metrics for one superposed match.
#[derive(Debug, Clone, Default, PartialEq, Copy)]
pub struct StructureSimilarityMetrics {
    pub tm_score: f32,
    pub gdt_ts: f32,
    pub gdt_ha: f32,
    pub chamfer_distance: f32,
    pub hausdorff_distance: f32,
    /// Distance-matrix RMSD: deformation of the match without superposing it
    pub drmsd: f32,
    /// Worst single internal-distance deviation, in Angstroms
    pub max_dist_deviation: f32,
}

impl StructureSimilarityMetrics {

    /// All metrics zeroed.
    pub fn new() -> Self {
        Self {
            tm_score: 0.0,
            gdt_ts: 0.0,
            gdt_ha: 0.0,
            chamfer_distance: 0.0,
            hausdorff_distance: 0.0,
            drmsd: 0.0,
            max_dist_deviation: 0.0,
        }
    }

    pub fn calculate_tm_score(&self, precomputed: &PrecomputedDistances) -> f32 {
        tm_score(precomputed, None)
    }

    pub fn calculate_gdt_ts(&self, precomputed: &PrecomputedDistances) -> f32 {
        gdt_ts(precomputed)
    }

    pub fn calculate_gdt_ha(&self, precomputed: &PrecomputedDistances) -> f32 {
        gdt_ha(precomputed)
    }

    pub fn calculate_chamfer_distance(&self, precomputed: &PrecomputedDistances) -> f32 {
        chamfer_distance(precomputed)
    }

    pub fn calculate_hausdorff_distance(&self, precomputed: &PrecomputedDistances) -> f32 {
        hausdorff_distance(precomputed)
    }

    /// Fill TM-score, GDT-TS/HA, Chamfer and Hausdorff. dRMSD fields are set by the caller.
    pub fn calculate_all(&mut self, precomputed: &PrecomputedDistances) {
        self.tm_score = self.calculate_tm_score(precomputed);
        self.gdt_ts = self.calculate_gdt_ts(precomputed);
        self.gdt_ha = self.calculate_gdt_ha(precomputed);
        self.chamfer_distance = self.calculate_chamfer_distance(precomputed);
        self.hausdorff_distance = self.calculate_hausdorff_distance(precomputed);
    }

    pub fn print_in_a_formatted_way(&self) {
        println!("Structure Similarity Metrics:");
        println!("  TM-score:           {:.4}", self.tm_score);
        println!("  GDT-TS:             {:.4}", self.gdt_ts);
        println!("  GDT-HA:             {:.4}", self.gdt_ha);
        println!("  Chamfer Distance:   {:.4} Å", self.chamfer_distance);
        println!("  Hausdorff Distance: {:.4} Å", self.hausdorff_distance);
        println!("  dRMSD:              {:.4} Å", self.drmsd);
        println!("  Max dist deviation: {:.4} Å", self.max_dist_deviation);
    }
}

impl fmt::Display for StructureSimilarityMetrics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            self.tm_score,
            self.gdt_ts,
            self.gdt_ha,
            self.chamfer_distance,
            self.hausdorff_distance,
            self.drmsd,
            self.max_dist_deviation
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::structure::chain_id::ChainId;
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
        assert!((metrics.gdt_ts - 1.0).abs() < 1e-6);
        assert!((metrics.gdt_ha - 1.0).abs() < 1e-6);
        assert!(metrics.chamfer_distance.abs() < 1e-6);
        assert!(metrics.hausdorff_distance.abs() < 1e-6);
        
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
            query_zinc_structure.get_index(&ChainId::from_byte(b'F'), &207).unwrap(),
            query_zinc_structure.get_index(&ChainId::from_byte(b'F'), &212).unwrap(),
            query_zinc_structure.get_index(&ChainId::from_byte(b'F'), &225).unwrap(),
            query_zinc_structure.get_index(&ChainId::from_byte(b'F'), &229).unwrap(),
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
            target_zinc_structure.get_index(&ChainId::from_byte(b'A'), &257).unwrap(),
            target_zinc_structure.get_index(&ChainId::from_byte(b'A'), &262).unwrap(),
            target_zinc_structure.get_index(&ChainId::from_byte(b'A'), &275).unwrap(),
            target_zinc_structure.get_index(&ChainId::from_byte(b'A'), &279).unwrap(),
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
            query_zinc_structure.get_index(&ChainId::from_byte(b'F'), &res_num).unwrap()
        }).collect::<Vec<usize>>();
        reference_indices.extend((225..230).map(|res_num| {
            query_zinc_structure.get_index(&ChainId::from_byte(b'F'), &res_num).unwrap()
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
            target_zinc_structure.get_index(&ChainId::from_byte(b'A'), &res_num).unwrap()
        }).collect::<Vec<usize>>();
        target_indices.extend((275..280).map(|res_num| {
            target_zinc_structure.get_index(&ChainId::from_byte(b'A'), &res_num).unwrap()
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
            query_zinc_structure.get_index(&ChainId::from_byte(b'F'), &205).unwrap(), // Outlier
            query_zinc_structure.get_index(&ChainId::from_byte(b'F'), &212).unwrap(),
            query_zinc_structure.get_index(&ChainId::from_byte(b'F'), &225).unwrap(),
            query_zinc_structure.get_index(&ChainId::from_byte(b'F'), &229).unwrap(),
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
            target_zinc_structure.get_index(&ChainId::from_byte(b'A'), &257).unwrap(),
            target_zinc_structure.get_index(&ChainId::from_byte(b'A'), &262).unwrap(),
            target_zinc_structure.get_index(&ChainId::from_byte(b'A'), &275).unwrap(),
            target_zinc_structure.get_index(&ChainId::from_byte(b'A'), &279).unwrap(),
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
