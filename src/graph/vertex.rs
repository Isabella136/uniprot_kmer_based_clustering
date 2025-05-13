use crate::Graph;
use crate::graph::edge::KmerEdge;
// use crate::graph::edge::KmerEdgeHelper;

use std::fmt;
use std::sync::{Arc, Weak};
use std::sync::atomic::Ordering::Acquire;
use std::sync::atomic::{AtomicPtr, AtomicUsize};

// type ArcUsize = Arc<AtomicUsize>;
type ArcKmerEdge = Arc<AtomicPtr<KmerEdge>>;
// type ArcKmerEdgeHelper = Arc<AtomicPtr<KmerEdgeHelper>>;

type WeakUsize = Weak<AtomicUsize>;
type WeakGraph = Weak<AtomicPtr<Graph>>;

pub struct ProteinVertex {
    key: usize,
    graph: WeakGraph,
    edges_key: Vec<WeakUsize>,
}
impl ProteinVertex {

    // Make a new ProteinVertex
    pub fn new(key: usize, graph: WeakGraph) -> ProteinVertex {
        // Pointer is valid, so we are safe
        let graph_ref: &Graph = unsafe {&*graph.upgrade().unwrap().load(Acquire)};
        let edges_key = graph_ref.protein_list[key].clone().get_five_hash().iter()
            .map(|k: &u32| Arc::downgrade(&graph_ref.global_edge_keys[*k as usize]))
            .collect();

        ProteinVertex{
            key, 
            graph, 
            edges_key,
        }
    }

    // Get key of edges vertex is in
    pub fn get_edges_key(&self) -> &Vec<WeakUsize>{
        &self.edges_key
    }

    // Add key to KmerEdge
    pub fn add_edge(&mut self, edge_key: WeakUsize) {
        // Pointer is valid, so we are safe
        let graph_ref: &Graph = unsafe {&*self.graph.upgrade().unwrap().load(Acquire)};

        let length = graph_ref.edges.len();
        let edges_bit_array: Vec<bool> = self.get_edges_bit_array(length);
        if !edges_bit_array[edge_key.upgrade().unwrap().load(Acquire)] {
            self.edges_key.push(edge_key);
        }
    }

    // Remove keys to KmerEdge
    pub fn remove_edges(&mut self) {
        let mut index = 0usize;
        while index < self.edges_key.len() {
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
        // Pointer is valid, so we are safe
        let graph_ref: &Graph = unsafe {&*self.graph.upgrade().unwrap().load(Acquire)};

        for key in &self.edges_key {
            let loaded_key = key.upgrade().unwrap().load(Acquire);
            let aptr_curr_edge: ArcKmerEdge = 
                graph_ref.edges[loaded_key].clone();

            // We won't update kmer edge, so we are safe
            let curr_edge: &KmerEdge = unsafe {
                &*aptr_curr_edge.load(Acquire)};

            curr_edge.add_vertex(self.key);
            
        }
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
