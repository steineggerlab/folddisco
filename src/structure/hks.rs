// This module implements the HKS (Heat Kernel Signature) computation for protein structures.

use crate::structure::distance_graph::{
    AdjacencyGraph,
};
use graphalgs::spec::laplacian_matrix;
use nalgebra::{DMatrix, DVector, SymmetricEigen};

pub struct HksEngine {
    pub eigenvalues: DVector<f64>,
    pub eigenvectors: DMatrix<f64>,
    pub num_nodes: usize,
}

impl HksEngine {
    pub fn from_graph(graph: &AdjacencyGraph) -> Self {
        let n = graph.node_count();

        // 1. Use graphalgs to get the Laplacian matrix (L = D - A)
        // For unweighted adjacency graphs, use unit edge cost.
        let laplacian_nalgebra = laplacian_matrix(graph, |_| 1.0_f64);

        // 3. Eigendecomposition (L is symmetric, so we use SymmetricEigen)
        let decomp = SymmetricEigen::new(laplacian_nalgebra);

        Self {
            eigenvalues: decomp.eigenvalues,
            eigenvectors: decomp.eigenvectors,
            num_nodes: n,
        }
    }

    /// Pointwise HKS: Signature for a single residue (used for node features)
    pub fn compute_hks_pointwise(&self, i: usize, time_scales: &[f64]) -> Vec<f64> {
        time_scales.iter().map(|&t| {
            let mut val = 0.0;
            for k in 0..self.num_nodes {
                let lambda = self.eigenvalues[k];
                let phi_ik = self.eigenvectors[(i, k)];
                val += (-lambda * t).exp() * phi_ik * phi_ik;
            }
            val
        }).collect()
    }

}

#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;
    use crate::structure::distance_graph::convert_distance_graph_to_adjacency_graph;
    use crate::PDBReader;

    #[test]
    fn test_hks_from_real_structure() {
        let pdb_list = vec!["data/1AL1.pdb", "data/4R80.pdb"];
        for pdb in pdb_list {
             let compact = PDBReader::from_file(pdb)
            .unwrap()
            .read_structure()
            .unwrap()
            .to_compact();
            let dist_graph = compact.to_distance_graph();
            let adj_graph = convert_distance_graph_to_adjacency_graph(
                &dist_graph, 5.8
            );
            let adj_graph = convert_distance_graph_to_adjacency_graph(
                &dist_graph, 8.0
            );
            
            let engine = HksEngine::from_graph(&adj_graph);
            let time_scales = vec![0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.8, 1.0, 2.0, 5.0];
            
            println!("num_nodes: {}", engine.num_nodes);
            let eigen_preview: Vec<f64> = engine.eigenvalues.iter().take(10).copied().collect();
            println!("eigenvalues(first up to 10): {:?}", eigen_preview);
            println!("time_scales: {:?}", time_scales);
            
            for i in 0..engine.num_nodes {
                let hks_i = engine.compute_hks_pointwise(i, &time_scales);
                println!("{}, {:?}", i, hks_i);
            }
        }
    }
}