use crate::Graph;
use crate::graph::edge::KmerEdge;
use crate::graph::edge::KmerEdgeHelper;

use std::fmt;
use std::sync::{Arc, Weak};
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering::Acquire;

// type ArcUsize = Arc<AtomicUsize>;
type ArcKmerEdge = Arc<KmerEdge>;
type ArcKmerEdgeHelper = Arc<KmerEdgeHelper>;

type WeakUsize = Weak<AtomicUsize>;
type WeakGraph = Weak<Graph>;

pub struct ProteinVertex {
    key: usize,
    graph: WeakGraph,
    edges_key: Vec<WeakUsize>,
}
impl ProteinVertex {

    // Make a new ProteinVertex
    pub fn new(graph: WeakGraph) -> ProteinVertex {
        ProteinVertex{
            key: 0, 
            graph, 
            edges_key: vec![],
        }
    }

    // Get key of edges vertex is in
    pub fn get_edges_key(&self) -> &Vec<WeakUsize>{
        &self.edges_key
    }

    // Change key info and populate edges_key
    pub fn update(&mut self, key: usize) {
        // Make graph data accessible with temporary strong pointer
        let graph_strong = self.graph.upgrade().unwrap();

        self.key = key;
        self.edges_key = graph_strong.protein_list[key].clone().get_five_hash().iter()
            .map(|k: &u32| Arc::downgrade(&graph_strong.global_edge_keys[*k as usize]))
            .collect();
        
        // Remove temporary strong pointer
        drop(graph_strong);
    }

    // Add key to KmerEdge
    pub fn add_edge(&mut self, edge_key: WeakUsize) {
        // Make graph data accessible with temporary strong pointer
        let graph_strong = self.graph.upgrade().unwrap();

        let length = graph_strong.edges.len();
        let edges_bit_array: Vec<bool> = self.get_edges_bit_array(length);
        if !edges_bit_array[edge_key.upgrade().unwrap().load(Acquire)] {
            self.edges_key.push(edge_key);
        }

        // Remove temporary strong pointer
        drop(graph_strong);
    }

    // Remove keys to KmerEdge
    pub fn remove_edges(&mut self, removed_bit_array: Arc<Vec<bool>>) {
        let mut index = 0usize;
        while index < self.edges_key.len() {
            if self.edges_key[index].strong_count() > 0 {
                if removed_bit_array[self.edges_key[index].upgrade().unwrap().load(Acquire)] {
                    panic!("We didn't get rid of all the ARCs: index {} has {} ARC ptrs",
                        self.edges_key[index].upgrade().unwrap().load(Acquire),
                        self.edges_key[index].strong_count())
                }
            }
            if self.edges_key[index].strong_count() == 0 {
                self.edges_key.remove(index);
            }
            else {
                index += 1;
            }
        }
    }

    // Update ProteinVertex's edges to reflect ProteinVertex's presence
    pub fn update_graph_edges(&self) {
        // Make graph data accessible with temporary strong pointer
        let graph_strong = self.graph.upgrade().unwrap();

        for key in &self.edges_key {
            let loaded_key = key.upgrade().unwrap().load(Acquire);
            let curr_edge: &mut ArcKmerEdge = 
                &mut graph_strong.edges[loaded_key].clone();
            let curr_helper: ArcKmerEdgeHelper =
                graph_strong.edge_update_helpers[loaded_key].clone();

            if curr_helper.can_update_kmer_edge() {
                // Only one thread can access a specific edge at a time, so we are safe
                let curr_edge_raw: &mut KmerEdge = unsafe {
                    Arc::get_mut_unchecked(curr_edge)};
                curr_helper.update_kmer_edge(self.key.clone(), curr_edge_raw);
            }
            else {
                curr_helper.add_to_array(self.key.clone());
            }
        }

        // Remove temporary strong pointer
        drop(graph_strong);
    }

    // Get bit array of ProteinVertex's edges
    pub fn get_edges_bit_array(&self, length: usize) -> Vec<bool> {
        let mut edges_bit_array: Vec<bool> = vec![false; length];
        for key in &self.edges_key {
            edges_bit_array[key.upgrade().unwrap().load(Acquire)] = true;
        }
        edges_bit_array
    }
}

impl fmt::Debug for ProteinVertex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let edge_amt = self.edges_key.len();
        f.debug_struct("Protein")
                .field("key", &self.key)
                .field("size", &edge_amt).finish()
    }
}
