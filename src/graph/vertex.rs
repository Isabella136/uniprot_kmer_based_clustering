use crate::Graph;
use crate::graph::edge::KmerEdge;

use std::fmt;
use std::sync::{Arc, Weak};
use std::sync::atomic::{AtomicPtr, AtomicUsize};
use std::sync::atomic::Ordering::{Acquire, AcqRel};

type ArcUsize = Arc<AtomicUsize>;
type ArcKmerEdge = Arc<AtomicPtr<KmerEdge>>;

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
        let edges_key = vec![];

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

    // Only keep edges if they are included in function parameter
    pub fn keep_specified_edges(&mut self, specified_edge_keys: Arc<Vec<ArcUsize>>) {
        // Pointer is valid, so we are safe
        let graph_ref: &Graph = unsafe {&*self.graph.upgrade().unwrap().load(Acquire)};

        let length = graph_ref.edges.len();
        let old_edges_bit_array: Vec<bool> = self.get_edges_bit_array(length);
        let mut new_edges_key_vec: Vec<WeakUsize> = vec![];

        for index in 0..specified_edge_keys.len() {
            let key = specified_edge_keys[index].load(Acquire);
            if old_edges_bit_array[key] {
                new_edges_key_vec.push(Arc::downgrade(&specified_edge_keys[index]));
            }
        }

        self.edges_key = new_edges_key_vec;
    }

    // Update ProteinVertex's edges to reflect ProteinVertex's presence
    pub fn update_graph_edges(&mut self, edge_count_prefix_sum: Arc<Vec<usize>>, 
            times_kmer_visited: Arc<Vec<ArcUsize>>, kmer_freq: Arc<Vec<usize>>) {
        // Pointer is valid, so we are safe
        let graph_ref: &Graph = unsafe {&*self.graph.upgrade().unwrap().load(Acquire)};
        
        fn update_edge_info<'a>(vertex: &'a mut ProteinVertex, graph_ref: &'a Graph, edge_index: usize, kmer: &'a usize) -> Result<&'a str, &'a str>{
            let edge_ptr: ArcKmerEdge = graph_ref.edges[edge_index].clone();
            
            // add vertex's key to edge
            // we are not modifying ptr, so we are safe
            let edge: &KmerEdge = unsafe {& *edge_ptr.load(Acquire)};
            if !edge.get_kmers().contains(kmer) {
                return Result::Err("Math error yet again")
            }

            let res = edge.add_vertex(vertex.key);
            
            // add edge's key to vertex
            vertex.edges_key.push(Arc::downgrade(&graph_ref.global_edge_keys[edge_index]));

            res
        }

        let mut kmers: Vec<usize> = graph_ref.protein_list[self.key].clone().get_five_hash()
            .iter().map(|x| *x as usize).collect();
        kmers.sort();
        kmers.dedup();        
        
        // for each kmer:
            // number of batches = kmer_freq[kmer] - 1
            // length of batch i = kmer_freq[kmer] - 1 - i
        for kmer in kmers {
            // calculate the left-most edge assigned to the current kmer
            let left_edge = match kmer {
                0 => 0,
                _ => edge_count_prefix_sum[kmer - 1]
            };
            
            // gets the number of times the kmer was already "visited"
            // equal to the index of the batch that is now assigned to the current vertex
            // also equal to the number of already-assigned batches we need to go through
            let visited = times_kmer_visited[kmer].fetch_add(1, AcqRel);
            
            // for each batch
            for batch in 0..visited+1 {
                // offset corresponds to the total amount of edges in previous batches 
                let offset = match batch {
                    0 => 0,
                    _ => (0..batch).map(|x| kmer_freq[kmer] - 1 - x).sum()
                };
                // if current batch is the one assigned to current vertex
                // add vertex to all edges in current batch
                if batch == visited {
                    for j in 0..kmer_freq[kmer] - 1 - batch {
                        let edge_index = left_edge + offset + j;
                        let res = update_edge_info(self, graph_ref, edge_index, &kmer);
                        if res.is_err() {
                            let edge = unsafe {& *graph_ref.edges[edge_index].clone().load(Acquire)};
                            panic!("I did my math wrong again\n\tvertex: {},\tkmer: {},\tvisited: {},\tkmer_freq: {},\tedge: {},\tedge_kmer: {:?},\tprevious vertices: {:?}",
                                self.key, kmer, visited, kmer_freq[kmer], edge_index, 
                                edge.get_kmers(), edge.get_vertices_key());
                        }
                    }
                }
                // else, add vertex to only one edge in current batch
                else {
                    // sub_offset corresponds to total number of edges in current batches
                    // that already have two vertices
                    let sub_offset =  visited - 1 - batch;
                    let edge_index = left_edge + offset + sub_offset;
                    let res = update_edge_info(self, graph_ref, edge_index, &kmer);
                    if res.is_err() {
                        let edge = unsafe {& *graph_ref.edges[edge_index].clone().load(Acquire)};
                        panic!("I did my math wrong again\n\tvertex: {},\tkmer: {},\tvisited: {},\tkmer_freq: {},\tedge: {},\tedge_kmer: {:?},\tprevious vertices: {:?}",
                            self.key, kmer, visited, kmer_freq[kmer], edge_index, 
                            edge.get_kmers(), edge.get_vertices_key());
                    }
                }
            }

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

    // Get id and sequence from original protein struct
    pub fn get_protein_id_and_seq(&self) -> (&String, &String) {
        // Pointer is valid, so we are safe
        let graph_ref: &Graph = unsafe {&*self.graph.upgrade().unwrap().load(Acquire)};
        return graph_ref.protein_list[self.key].get_id_and_seq()
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
