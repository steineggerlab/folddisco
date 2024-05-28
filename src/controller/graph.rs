// Graph representation of residue indices to match positions in the structure
//
// 2024-05-09 13:40:20
// Current naive implementation:
// Find both strong and weakly connected components with same node count as query graph

use std::collections::HashMap;

use petgraph::graph::DiGraph;
use crate::geometry::core::GeometricHash;


pub fn create_index_graph(
    ind_vec: &Vec<Vec<(usize, usize)>>, hash_vec: &Vec<GeometricHash>
) -> DiGraph<usize, GeometricHash> {
    let mut graph = DiGraph::<usize, GeometricHash>::new();
    let mut node_indices = HashMap::new();

    for (inner_vec, hash) in ind_vec.iter().zip(hash_vec.iter()) {
        if inner_vec.is_empty() {
            continue;
        }
        inner_vec.iter().for_each(|(i, j)| {
            let node_index_i = *node_indices.entry(*i).or_insert_with(|| graph.add_node(*i));
            let node_index_j = *node_indices.entry(*j).or_insert_with(|| graph.add_node(*j));
            graph.add_edge(node_index_i, node_index_j, hash.clone());
        });
    }
    graph
}

pub fn connected_components_with_given_node_count(
    graph: &DiGraph<usize, GeometricHash>, node_count: usize
) -> Vec<Vec<usize>> {
    // Current implementation.
    // Return both strong and weakly connected components with the same node count
    let mut scc = petgraph::algo::tarjan_scc(graph);
    let undirected_graph = graph.clone().into_edge_type::<petgraph::Undirected>();
    let mut wcc = petgraph::algo::kosaraju_scc(&undirected_graph);
    // Concat
    scc.append(&mut wcc);
    // Filter components with the same node count
    // scc.retain(|component| component.len() == node_count);
    // Filter components with the node count greater than or equal to the given node count
    // scc.retain(|component| component.len() >= 3);

    // Uniqueness cheking. Sort all inner vectors and dedup
    scc.iter_mut().for_each(|component| component.sort());
    scc.sort();
    scc.dedup();
    let output = scc.iter().map(
        |component| component.iter().map(|&node| graph[node]
    ).collect()).collect();
    output
}


#[cfg(test)]
mod tests {
    use petgraph::Undirected;
    use crate::measure_time;
    use crate::prelude::*;
    use super::*;

    #[test]
    fn test_graph() {
        let mut graph = DiGraph::<&str, u32>::new();

        let b102 = graph.add_node("B102");
        println!("{:?}", graph[b102]);
        let b57 = graph.add_node("B57");
        let c195 = graph.add_node("C195");
        let f57 = graph.add_node("F57");
        let f102 = graph.add_node("F102");
        let g195 = graph.add_node("G195");

        graph.add_edge(b102, b57, 109329223);
        graph.add_edge(b102, c195, 116878724);
        graph.add_edge(b57, b102, 271858548);
        graph.add_edge(f57, f102, 271858548);
        graph.add_edge(b57, c195, 284511716);
        graph.add_edge(f57, g195, 284511716);
        graph.add_edge(c195, b102, 506948936);
        graph.add_edge(c195, b57, 512052558);
        graph.add_edge(g195, f102, 512052558);
        let scc = measure_time!(petgraph::algo::tarjan_scc(&graph));
        println!("{:?}", scc);
        
        let undirected_graph = graph.clone().into_edge_type::<Undirected>();
        let weak_cc = measure_time!(petgraph::algo::kosaraju_scc(&undirected_graph));
        println!("{:?}", weak_cc);
    }

    #[test]
    fn test_usize_graph() {
        let ind_vec = vec![
            vec![(57, 102), (57, 195)], 
            vec![(102, 57), (102, 195)],
            vec![],
            vec![(195, 57), (195, 102)], // Triad
            vec![(57, 106), (106, 111), (212, 215), (219, 417), (417, 219)],
            vec![(219, 212), (212, 417)],
            vec![],
        ];
        let hash_vec = vec![GeometricHash::from_u32(1, HashType::PDBTrRosetta); 7];
        let graph = create_index_graph(&ind_vec, &hash_vec);
        println!("{:?}", graph);
        let conn = connected_components_with_given_node_count(&graph, 3);
        println!("{:?}", conn);
    }
    
    
}