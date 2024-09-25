use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Node {
    pub index: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Edge {
    pub source: Node,
    pub target: Node,
    pub edge_weight: u32,
}

#[derive(Debug, Clone)]
pub struct Graph {
    pub nodes: HashMap<usize, Node>,
    pub edges: Vec<Edge>,
    pub adj_list: HashMap<usize, HashSet<usize>>, // adjacency list for the nodes
}

impl Graph {
    pub fn new() -> Self {
        Graph {
            nodes: HashMap::new(),
            edges: Vec::new(),
            adj_list: HashMap::new(),
        }
    }

    pub fn add_node(&mut self, index: usize) -> Node {
        let node = Node { index };
        self.nodes.insert(index, node.clone());
        self.adj_list.entry(index).or_insert_with(HashSet::new);
        node
    }

    pub fn add_edge(&mut self, source: Node, target: Node, edge_weight: u32) {
        self.edges.push(Edge {
            source: source.clone(),
            target: target.clone(),
            edge_weight,
        });
        self.adj_list.get_mut(&source.index).unwrap().insert(target.index);
        self.adj_list.get_mut(&target.index).unwrap().insert(source.index);
    }

    pub fn get_neighbors(&self, index: usize) -> Option<&HashSet<usize>> {
        self.adj_list.get(&index)
    }
}

fn find_maximal_common_subgraph(graph1: &Graph, graph2: &Graph) -> Graph {
    let mut common_graph = Graph::new();
    let mut node_map = HashMap::new();

    // Add common nodes to the new graph
    for node in graph1.nodes.values() {
        if graph2.nodes.contains_key(&node.index) {
            let common_node = common_graph.add_node(node.index);
            node_map.insert(node.index, common_node);
        }
    }

    // Add common edges based on edge weights
    for edge1 in &graph1.edges {
        for edge2 in &graph2.edges {
            if edge1.edge_weight == edge2.edge_weight {
                if let (Some(source), Some(target)) = (node_map.get(&edge1.source.index), node_map.get(&edge1.target.index)) {
                    common_graph.add_edge(source.clone(), target.clone(), edge1.edge_weight);
                }
            }
        }
    }

    common_graph
}

#[cfg(test)]
mod tests {
    use crate::measure_time;

    use super::*;

    #[test]
    fn test_find_maximal_common_subgraph() {
        measure_time!({
        // Create first graph
        let mut graph1 = Graph::new();
        let b102 = graph1.add_node(1);
        let b57 = graph1.add_node(2);
        let c195 = graph1.add_node(3);
        let f57 = graph1.add_node(4);
        let f102 = graph1.add_node(5);
        let g195 = graph1.add_node(6);

        graph1.add_edge(b102.clone(), b57.clone(), 109329223);
        graph1.add_edge(b102.clone(), c195.clone(), 116878724);
        graph1.add_edge(b57.clone(), b102.clone(), 271858548);
        graph1.add_edge(f57.clone(), f102.clone(), 271858548);
        graph1.add_edge(b57.clone(), c195.clone(), 284511716);
        graph1.add_edge(f57.clone(), g195.clone(), 284511716);
        graph1.add_edge(c195.clone(), b102.clone(), 506948936);
        graph1.add_edge(c195.clone(), b57.clone(), 512052558);
        graph1.add_edge(g195.clone(), f102.clone(), 512052558);

        // Create second graph with only B102, B57, C195
        let mut graph2 = Graph::new();
        let b102_2 = graph2.add_node(1);
        let b57_2 = graph2.add_node(2);
        let c195_2 = graph2.add_node(3);

        graph2.add_edge(b102_2.clone(), b57_2.clone(), 109329223);
        graph2.add_edge(b102_2.clone(), c195_2.clone(), 116878724);
        graph2.add_edge(b57_2.clone(), b102_2.clone(), 271858548);
        graph2.add_edge(b57_2.clone(), c195_2.clone(), 284511716);

        // Find maximal common subgraph
        let common_subgraph = find_maximal_common_subgraph(&graph1, &graph2);

        // Expected common subgraph
        let mut expected_common_subgraph = Graph::new();
        let b102_exp = expected_common_subgraph.add_node(1);
        let b57_exp = expected_common_subgraph.add_node(2);
        let c195_exp = expected_common_subgraph.add_node(3);

        expected_common_subgraph.add_edge(b102_exp.clone(), b57_exp.clone(), 109329223);
        expected_common_subgraph.add_edge(b102_exp.clone(), c195_exp.clone(), 116878724);
        expected_common_subgraph.add_edge(b57_exp.clone(), b102_exp.clone(), 271858548);
        expected_common_subgraph.add_edge(b57_exp.clone(), c195_exp.clone(), 284511716);

        // Assert common subgraph
        assert_eq!(common_subgraph.nodes, expected_common_subgraph.nodes);
        assert_eq!(common_subgraph.edges, expected_common_subgraph.edges);
        });
    }
}
