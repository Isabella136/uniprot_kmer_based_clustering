use crate::Graph;
use crate::graph::vertex::ProteinVertex;

use std::fmt;
use std::sync::{Arc, Weak};
use std::sync::atomic::{AtomicPtr, AtomicUsize};
use std::sync::atomic::Ordering::{Acquire, AcqRel, Release};

type ArcUsize = Arc<AtomicUsize>;
type ArcProteinVertex = Arc<AtomicPtr<ProteinVertex>>;

type WeakUsize = Weak<AtomicUsize>;
type WeakGraph = Weak<AtomicPtr<Graph>>;

pub struct KmerEdgeSingle {
    kmer: usize,
    vertices_key: [ArcUsize; 2],
    times_visited: AtomicUsize,
}
impl KmerEdgeSingle {

    // Make a new KmerEdge
    fn new(kmer: usize) -> KmerEdgeSingle {
        KmerEdgeSingle {
            kmer,
            vertices_key: [Arc::new(AtomicUsize::new(0)), 
                Arc::new(AtomicUsize::new(0))],
            times_visited: AtomicUsize::new(0),
        }
    }

    // Add key to ProteinVertex
    fn add_vertex(&self, vertex_key: usize) -> Result<&str, &str> {
        let visited = self.times_visited.fetch_add(1, AcqRel);
        if 0 == visited {
            self.vertices_key[0].store(vertex_key, Release);
        }
        else if 1 == visited {
            self.vertices_key[1].store(vertex_key, Release);
        }
        else {
            return Result::Err("I did my math wrong again")
        }
        Result::Ok("Yeah")
    }
}

pub struct KmerEdgeGroup {
    kmers: Vec<usize>,
    graph: WeakGraph,
    vertices_key: [ArcUsize; 2],
}
impl KmerEdgeGroup {

    // Merge KmerEdge objects into one
    fn new(kmer_edge_keys: &Vec<WeakUsize>, graph: WeakGraph) -> KmerEdgeGroup {
        // Create empty vectors
        let mut kmers: Vec<usize> = vec![];
        let vertices_key = [Arc::new(AtomicUsize::new(0)), 
                Arc::new(AtomicUsize::new(0))];

        let mut keys_are_set = false;
        // Pointer is valid, so we are safe
        let graph_ref: &Graph = unsafe {&*graph.upgrade().unwrap().load(Acquire)};

        // Append each KmerEdge's kmer and vertices info to aforementioned vectors
        for edge_key in kmer_edge_keys {
            let loaded_edge_key = edge_key.upgrade().unwrap().load(Acquire);
            let a_ptr_kmer_edge = graph_ref.edges[loaded_edge_key].clone();

            // We won't update kmer_edge, so we are safe
            let kmer_edge = unsafe{&*a_ptr_kmer_edge.load(Acquire)};

            kmers.append(&mut kmer_edge.get_kmers());
            if !keys_are_set {
                keys_are_set = true;
                let kmer_edge_vertices = kmer_edge.get_vertices_key();
                vertices_key[0].store(kmer_edge_vertices[0], Release);
                vertices_key[1].store(kmer_edge_vertices[1], Release);
            }
        }

        // Make KmerEdgeGroup object
        KmerEdgeGroup {kmers, graph, vertices_key}
    }

    fn get_proteins_ids_and_sequences(&self) -> [(&String, &String); 2] {
        // Pointer is valid, so we are safe
        let graph_ref: &Graph = unsafe {&*self.graph.upgrade().unwrap().load(Acquire)};

        // Get start protein reference
        let start_protein_key: usize = 
            self.vertices_key[0].load(Acquire);
        let start_protein_ptr: ArcProteinVertex = 
            graph_ref.vertices[start_protein_key].clone();
        let start_protein_ref: &ProteinVertex = 
            unsafe {&*start_protein_ptr.load(Acquire)};
        
        
        // Get end protein reference
        let end_protein_key: usize = 
            self.vertices_key[1].load(Acquire);
        let end_protein_ptr: ArcProteinVertex = 
            graph_ref.vertices[end_protein_key].clone();
        let end_protein_ref: &ProteinVertex = 
            unsafe {&*end_protein_ptr.load(Acquire)};

        return [start_protein_ref.get_protein_id_and_seq(), 
            end_protein_ref.get_protein_id_and_seq()]
    }

}

pub enum KmerEdge {
    Single(KmerEdgeSingle),
    Group(KmerEdgeGroup),
}
impl KmerEdge {
    pub fn get_kmers(&self) -> Vec<usize> {
        match self {
            Self::Single(val) => vec![val.kmer],
            Self::Group(val) => val.kmers.clone(),
        }
    }

    pub fn get_vertices_key(&self) -> [usize; 2] {
        match self {
            Self::Single(val) => val.vertices_key.clone().map(|x| x.load(Acquire)),
            Self::Group(val) => val.vertices_key.clone().map(|x| x.load(Acquire)),
        }
    }

    pub fn add_vertex(&self, vertex_key: usize) -> Result<&str, &str> {
        match self {
            Self::Single(val) => val.add_vertex(vertex_key),
            Self::Group(_) => panic!(
                "Can't call add_vertex on kmer group"),
        }
    }

    pub fn new(kmer: usize) -> KmerEdge {
        Self::Single(KmerEdgeSingle::new(kmer))
    }
    
    pub fn merge(kmer_edge_keys: &Vec<WeakUsize>, graph: WeakGraph) -> KmerEdge {
        Self::Group(KmerEdgeGroup::new(kmer_edge_keys, graph))
    }

    pub fn get_proteins_ids_and_sequences(&self) -> [(&String, &String); 2]{
        match self {
            Self::Group(val) => val.get_proteins_ids_and_sequences(),
            _ => panic!("How?!")
        }
    }
    
}

impl fmt::Debug for KmerEdge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Single(val) => {
                let vertices_amt = val.vertices_key.len();
                f.debug_struct("Single Kmer")
                    .field("kmer", &val.kmer)
                    .field("size", &vertices_amt).finish()
            },
            Self::Group(val) => {
                let vertices_amt = val.vertices_key.len();
                f.debug_struct("Kmer Group")
                    .field("kmer", &val.kmers)
                    .field("size", &vertices_amt).finish()
            },
        }
    }
}