// This module defines a distance graph 
// Will be used to define an adjacency graph after

use crate::structure::core::CompactStructure;

// Using petgraph for graph representation
use petgraph::{graph::{NodeIndex, UnGraph}, visit::EdgeRef};


// Node as usize (index of the residue), edge weight as f32
pub type DistanceGraph = UnGraph<usize, f32>;
pub type AdjacencyGraph = UnGraph<usize, ()>;

impl CompactStructure {
    // Build a C-alpha distance graph from the compact structure
    // Mode can be used to specify different distance metrics (e.g., C-alpha, C-beta, etc.)
    // 2026-03-03 15:58:04 For simplicity, C-alpha distances here.
    pub fn to_distance_graph(&self) -> DistanceGraph {
        let mut graph = DistanceGraph::new_undirected();
        let num_residues = self.num_residues;
        let mode = 0u8;

        let node_vector: Vec<NodeIndex> = (0..num_residues).map(|i| graph.add_node(i)).collect();

        // Add edges with weights based on distance
        for i in 0..num_residues {
            for j in (i + 1)..num_residues {
                
                let dist = if mode == 0 { 
                    self.get_ca_distance(i, j)
                } else {
                    self.get_ca_distance(i, j) // Currently only support C-alpha distance.
                };
                if let Some(dist) = dist {
                    graph.add_edge(node_vector[i], node_vector[j], dist);
                }
            }
        }
        graph
    }
}


// Convert distance graph to adjacency graph based on a distance threshold
pub fn convert_distance_graph_to_adjacency_graph(distance_graph: &DistanceGraph, threshold: f32) -> AdjacencyGraph {
    // Default threshold: let's use 8.0 Å for C-alpha distance, which is a common cutoff for residue contacts
    let mut adj_graph = AdjacencyGraph::new_undirected();

    // Add nodes to the adjacency graph
    for node in distance_graph.node_indices() {
        adj_graph.add_node(*distance_graph.node_weight(node).unwrap());
    }
    
    for edge in distance_graph.edge_references() {
        let weight = edge.weight();
        if *weight <= threshold {
            let source = edge.source();
            let target = edge.target();
            adj_graph.add_edge(source, target, ());
        }
    }
    adj_graph
}


// Test
#[cfg(test)]
mod tests {
    use super::*;
    use crate::PDBReader;

    #[test]
    fn test_distance_graph() {
        let compact = PDBReader::from_file("query/4CHA.pdb").unwrap().read_structure().unwrap().to_compact();
        let dist_graph: DistanceGraph = compact.to_distance_graph();
        let adj_graph = convert_distance_graph_to_adjacency_graph(&dist_graph, 8.0);
        // println!("Distance Graph: {:?}", dist_graph);
        println!("Adjacency Graph: {:?}", adj_graph);
    }
}