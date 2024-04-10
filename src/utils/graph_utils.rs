use std::collections::HashMap;

#[derive(Debug)]
pub(crate) struct Graph<N> {
    nodes: Vec<N>,
    neighbors: Vec<Vec<i32>>,
}

impl<N> Graph<N>
where
    N: std::cmp::Eq + std::hash::Hash,
{
    pub(crate) fn new() -> Self {
        Self {
            nodes: Vec::new(),
            neighbors: Vec::new(),
        }
    }

    // A public getter for an unmodifiable list of nodes currently held in the graph for iteration.
    pub(crate) fn get_nodes(&self) -> &[N] {
        &self.nodes
    }

    /* directed and private */
    fn add_neighbor(&mut self, from_node: i32, to_node: i32) {
        let from_nodes_neighbors = match self.neighbors.get_mut(from_node as usize) {
            Some(v) => v,
            None => {
                panic!(
                    "len of vector, self.neighbors={} but requested idx={}",
                    self.neighbors.len(),
                    from_node
                );
            }
        };

        if !from_nodes_neighbors.contains(&to_node) {
            from_nodes_neighbors.push(to_node);
        }
    }

    pub(crate) fn add_node(&mut self, singleton: N) -> i32 {
        if let Some(p) = self.nodes.iter().position(|n| n.eq(&singleton)) {
            p as i32
        } else {
            self.nodes.push(singleton);
            self.neighbors.push(vec![]);

            self.nodes.len() as i32 - 1
        }
    }

    /* bidirectional and public */
    pub(crate) fn add_edge(&mut self, left: N, right: N) {
        let lr_equal = left == right;

        let left_index = self.add_node(left);

        if lr_equal {
            return;
        }

        let right_index = self.add_node(right);

        self.add_neighbor(left_index, right_index);
        self.add_neighbor(right_index, left_index);
    }

    // Part of Union-Find with Path Compression that joins to nodes to be part of the same cluster.
    fn join_nodes(grouping: &mut [i32], node_id1: i32, node_id2: i32) {
        let rep_node1 = Self::find_rep_node(grouping, node_id1);
        let rep_node2 = Self::find_rep_node(grouping, node_id2);

        if rep_node1 == rep_node2 {
            return;
        }

        grouping[rep_node1 as usize] = rep_node2;
    }

    // Part of Union-Find with Path Compression to determine the duplicate set a particular UMI belongs to.
    fn find_rep_node(grouping: &mut [i32], mut node_id: i32) -> i32 {
        let mut representative_umi = node_id; // All UMIs of a duplicate set will have the same reprsentativeUmi.

        while representative_umi != grouping[representative_umi as usize] {
            representative_umi = grouping[representative_umi as usize];
        }

        while node_id != representative_umi {
            let new_umi_id = grouping[node_id as usize];
            grouping[node_id as usize] = representative_umi;
            node_id = new_umi_id;
        }

        representative_umi
    }
    pub(crate) fn cluster(&self) -> HashMap<&N, i32> {
        let mut cluster = (0..self.nodes.len() as i32).collect::<Vec<i32>>();

        (0..self.neighbors.len() as i32).for_each(|i| {
            self.neighbors
                .get(i as usize)
                .unwrap()
                .iter()
                .copied()
                .for_each(|j| Self::join_nodes(&mut cluster, j, i))
        });

        // Must call findRepNode() here in case nodes are orphaned and become only transitively equal to roots in their cluster.
        self.nodes
            .iter()
            .enumerate()
            .map(|(node_id, n)| (n, Self::find_rep_node(&mut cluster, node_id as i32)))
            .collect::<HashMap<&N, i32>>()
    }
}

// pub(crate) trait GraphExt<'t, T:'t> {
//     /**
//      * returns the cluster map of connected components
//      *
//      * @return Nodes that point to the same integer are in the same cluster.
//      */
//     fn cluster(&'t self) -> HashMap<T, i32>;
// }

// impl<'s, N> GraphExt<'s, &'s N> for Graph<N>
// where
//     N: std::cmp::Eq + std::hash::Hash,
// {

//     fn cluster(&'s self) -> HashMap<&'s N, i32> {
//         let mut cluster = (0..self.nodes.len() as i32).collect::<Vec<i32>>();

//         (0..self.neighbors.len() as i32).for_each(|i| {
//             self.neighbors
//                 .get(i as usize)
//                 .unwrap()
//                 .iter()
//                 .copied()
//                 .for_each(|j| Self::join_nodes(&mut cluster, j, i))
//         });

//         // Must call findRepNode() here in case nodes are orphaned and become only transitively equal to roots in their cluster.
//         self.nodes
//             .iter()
//             .enumerate()
//             .map(|(node_id, n)| (n, Self::find_rep_node(&mut cluster, node_id as i32)))
//             .collect::<HashMap<&N, i32>>()
//     }
// }

// impl<'s, N> GraphExt<'s, N> for Graph<N>
// where
//     N: std::cmp::Eq + std::hash::Hash + Copy +'s,
// {

//     fn cluster(&'s self) -> HashMap<N, i32> {
//         let mut cluster = (0..self.nodes.len() as i32).collect::<Vec<i32>>();

//         (0..self.neighbors.len() as i32).for_each(|i| {
//             self.neighbors
//                 .get(i as usize)
//                 .unwrap()
//                 .iter()
//                 .copied()
//                 .for_each(|j| Self::join_nodes(&mut cluster, j, i))
//         });

//         // Must call findRepNode() here in case nodes are orphaned and become only transitively equal to roots in their cluster.
//         self.nodes
//             .iter()
//             .copied()
//             .enumerate()
//             .map(|(node_id, n)| (n, Self::find_rep_node(&mut cluster, node_id as i32)))
//             .collect::<HashMap<N, i32>>()
//     }
// }

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn simple_test() {
        let mut graph = Graph::<i32>::new();

        graph.add_edge(5, 6);
        graph.add_edge(5, 8);
        graph.add_edge(7, 9);
        graph.add_edge(7, 8);

        graph.add_edge(4, 2);
        graph.add_edge(1, 2);
        graph.add_edge(3, 1);

        let clusters = graph.cluster();

        // 9 nodes
        assert_eq!(clusters.len(), 9);

        let fives_cluster = clusters.get(&5).unwrap();
        // eprintln!("graph={:#?}", graph);
        // eprintln!("clusters={:#?}", clusters);
        (5..=9).for_each(|i| {
            assert_eq!(clusters.get(&i).unwrap(), fives_cluster);
        });

        let fours_cluster = clusters.get(&4).unwrap();
        (1..=4).for_each(|i| {
            assert_eq!(clusters.get(&i).unwrap(), fours_cluster);
        });

        assert_ne!(fives_cluster, fours_cluster);
    }

    #[test]
    fn second_test() {
        let mut graph = Graph::<i32>::new();

        graph.add_edge(1, 3);
        graph.add_edge(2, 3);
        graph.add_edge(0, 2);
        graph.add_edge(4, 3);

        graph.add_edge(5, 6);
        graph.add_edge(7, 6);
        graph.add_node(8);

        let clusters = graph.cluster();

        // 9 nodes
        assert_eq!(clusters.len(), 9);

        let threes_cluster = clusters.get(&3).unwrap();
        (0..5).for_each(|i| {
            assert_eq!(clusters.get(&i).unwrap(), threes_cluster);
        });

        let fives_cluster = clusters.get(&5).unwrap();
        (5..8).for_each(|i| {
            assert_eq!(clusters.get(&i).unwrap(), fives_cluster);
        });

        let eights_cluster = clusters.get(&8).unwrap();
        (8..9).for_each(|i| {
            assert_eq!(clusters.get(&i).unwrap(), eights_cluster);
        });

        assert_ne!(threes_cluster, fives_cluster);
        assert_ne!(fives_cluster, eights_cluster);
        assert_ne!(threes_cluster, eights_cluster);
    }

    #[test]
    fn pathological_input_ordering_test() {
        let mut graph = Graph::new();

        let [a, b, c, d, e, f, h, g] = (0_i32..=7)
            .into_iter()
            .collect::<Vec<i32>>()
            .try_into()
            .unwrap();
        graph.add_node(a);
        graph.add_node(b);
        graph.add_node(c);
        graph.add_node(d);
        graph.add_node(f);
        graph.add_node(h);
        graph.add_node(e);
        graph.add_node(g);

        // This graph is a hand constructed case where the algorithm produces in its final form an orphaned node that only
        // transitively points to the root for the rest of the graph (node d as it turns out)
        graph.add_edge(a, h);
        graph.add_edge(b, f);
        graph.add_edge(c, h);
        graph.add_edge(d, f);
        graph.add_edge(e, h);
        graph.add_edge(e, f);
        graph.add_edge(g, a);

        let clusters = graph.cluster();

        // 9 nodes
        assert_eq!(clusters.len(), 8);

        // Asserting that this connected graph has every node pointing to the same place
        let rep_cluster = clusters.get(&a).unwrap();
        graph.get_nodes().into_iter().for_each(|i| {
            assert_eq!(clusters.get(i).unwrap(), rep_cluster);
        });
    }
}
