use crate::protein::Protein;

use std::fmt;
use std::thread;
use std::sync::{Arc, Mutex, RwLock, RwLockReadGuard};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

struct KmerEdgeSingle {
    kmer: u32,
    vertices_key: Vec<usize>,
    vertices_bit_array: Vec<bool>,
    graph: Arc<Graph>,
}
impl KmerEdgeSingle {
    fn new(graph: Arc<Graph>, kmer: u32, vertices_count: usize) -> KmerEdgeSingle {
        KmerEdgeSingle {
            kmer,
            vertices_key: vec![],
            vertices_bit_array: vec![false; vertices_count],
            graph
        }
    }
    fn add_vertex(&mut self, vertex_key: usize) {
        if ! self.vertices_bit_array[vertex_key] {
            self.vertices_bit_array[vertex_key] = true;
            self.vertices_key.push(vertex_key);
        }
    }
}

struct KmerEdgeGroup {
    kmers: Vec<u32>,
    vertices_key: Vec<usize>,
    vertices_bit_array: Vec<bool>,
    graph: Arc<Graph>,
}
impl KmerEdgeGroup {
    // fn new(kmer_edge_one: KmerEdge, kmer_edge_two: KmerEdge) -> KmerEdgeGroup {
        
    // }
    fn add_vertex(&mut self, vertex_key: usize) {
        if ! self.vertices_bit_array[vertex_key] {
            self.vertices_bit_array[vertex_key] = true;
            self.vertices_key.push(vertex_key);
        }
    }
}

enum KmerEdge {
    Single(KmerEdgeSingle),
    Group(KmerEdgeGroup),
}
impl fmt::Debug for KmerEdge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Single(val) => {
                f.debug_list().entries(val.vertices_key.clone()).finish()
            },
            Self::Group(val) => {
                f.debug_list().entries(val.vertices_key.clone()).finish()
            },
        }
    }
}
impl KmerEdge {
    fn get_kmers(&self) -> Vec<u32> {
        match self {
            Self::Single(val) => vec![val.kmer],
            Self::Group(val) => val.kmers.clone(),
        }
    }
    fn get_vertices_key(&self) -> Vec<usize> {
        match self {
            Self::Single(val) => val.vertices_key.clone(),
            Self::Group(val) => val.vertices_key.clone(),
        }
    }
    fn get_vertices_bit_array(&self) -> Vec<bool> {
        match self {
            Self::Single(val) => val.vertices_bit_array.clone(),
            Self::Group(val) => val.vertices_bit_array.clone(),
        }
    }
    fn get_graph(&self) -> Arc<Graph> {
        match self {
            Self::Single(val) => val.graph.clone(),
            Self::Group(val) => val.graph.clone(),
        }
    }
    fn add_vertex(&mut self, vertex_key: usize) {
        match self {
            Self::Single(val) => val.add_vertex(vertex_key),
            Self::Group(val) => val.add_vertex(vertex_key),
        }
    }
    fn new(graph: Arc<Graph>, kmer: u32, vertices_count: usize) -> KmerEdge {
        Self::Single(KmerEdgeSingle::new(graph, kmer, vertices_count))
    }
}

struct ProteinVertex {
    protein: usize,
    edges_key: Vec<u32>,
    edges_bit_array: Vec<bool>,
    graph: Arc<Graph>,
}
impl ProteinVertex {
    fn new(protein: usize, graph: Arc<Graph>) -> ProteinVertex {
        let protein_obj = graph.protein_list[protein].clone();
        let edges_key = protein_obj.get_five_hash();
        let edges_bit_array = protein_obj.get_five_hash_map();
        ProteinVertex{protein, edges_key, edges_bit_array, graph}
    }
    fn update_graph_edges(&self) {
        for key in &self.edges_key {
            let mut edge = self.graph.edges[*key as usize].write().unwrap();
            edge.add_vertex(self.protein);
        }
    }
}

pub struct Graph {
    vertices: Arc<Vec<Arc<ProteinVertex>>>,
    edges: Arc<Vec<Arc<RwLock<KmerEdge>>>>,
    protein_list: Arc<Vec<Arc<Protein>>>,
}
impl fmt::Debug for Graph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.edges.iter()).finish()

    }
}
impl Graph {
    pub fn new(kmer_count: u32, thread_count: u32, protein_list: Arc<Vec<Arc<Protein>>>) -> Arc<Graph> {
        let protein_count_arc = protein_list.len();
        let graph = Graph {
            vertices: Arc::new(Vec::new()),
            edges: Arc::new(Vec::new()),
            protein_list,
        };
        let mut graph_arc = Arc::new(graph);
        unsafe {
            Arc::get_mut_unchecked(&mut graph_arc).edges = Arc::new((0..kmer_count)
                .into_iter().map(|k| Arc::new(RwLock::new(KmerEdge::new(graph_arc
                .clone(), k.clone(), protein_count_arc)))).collect())
        };
        let vertices: Arc<Vec<ProteinVertex>> = Arc::new((0..protein_count_arc)
            .into_iter().map(|protein| ProteinVertex::new(protein, graph_arc.clone()))
            .collect());
        let vertices_index = Arc::new(AtomicUsize::new(0usize));
        let mut handles = vec![];
        for _ in 0..thread_count {
            let vertices = vertices.clone();
            let vertices_index = vertices_index.clone();

            handles.push(thread::spawn(move || {
                loop {
                    let curr_vertices_index = vertices_index.fetch_add(
                        1, Ordering::AcqRel);
                    if curr_vertices_index >= 10619 {
                        break;
                    }
                    vertices[curr_vertices_index].update_graph_edges();
                }
            }))
        }
        for handle in handles {
            handle.join().unwrap();
        }
        graph_arc
    }
    pub fn combine_edges(mut graph: Arc<Graph>, thread_count: u32) {
        let mut key_reassignment: Vec<Option<u32>> = (0..graph.edges.len()).into_iter()
            .map(|k| Some(k as u32)).collect();
        let mut size_ordered_edges: Vec<usize> = (0..graph.edges.len()).into_iter()
            .collect();
        size_ordered_edges.sort_by(|&a, &b| {
            let edge_a: Arc<RwLock<KmerEdge>> = graph.edges[a].clone();
            let edge_b: Arc<RwLock<KmerEdge>> = graph.edges[b].clone();
            let edge_a_size = edge_a.read().unwrap().get_vertices_key().len();
            let edge_b_size = edge_b.read().unwrap().get_vertices_key().len();
            edge_a_size.cmp(&edge_b_size)});
        for edge_key in size_ordered_edges {
            if key_reassignment[edge_key].is_none() {
                continue;
            }
            let edge = graph.edges[edge_key].clone();
            let vertices_keys = Arc::new(edge.read().unwrap().get_vertices_key());
            let vertices_index = Arc::new(AtomicUsize::new(1usize));
            let other_edges_keys = Arc::new(Mutex::new(
                graph.vertices[vertices_keys[0]].edges_key.clone()));
            let other_edges_bitarray: Arc<Vec<AtomicBool>> = Arc::new(
                graph.vertices[vertices_keys[0]].edges_bit_array.iter()
                .map(|b| AtomicBool::new(*b)).collect());
            let mut handles = vec![];
            for _ in 0..thread_count {
                let graph = graph.clone();
                let vertices_keys = vertices_keys.clone();
                let vertices_index = vertices_index.clone();
                let other_edges_keys = other_edges_keys.clone();
                let other_edges_bitarray = other_edges_bitarray.clone();
                
                handles.push(thread::spawn(move || {
                    loop {
                        let curr_vertices_index = vertices_index.fetch_add(
                            1, Ordering::AcqRel);
                        let curr_vertex_key = vertices_keys[curr_vertices_index];
                        let curr_edges_keys = &graph.vertices[curr_vertex_key].edges_key;
                        for curr_key in curr_edges_keys {
                            if !other_edges_bitarray[*curr_key as usize].swap(
                                    true, Ordering::AcqRel) {
                                other_edges_keys.lock().unwrap().push(*curr_key);
                            }
                        }
                    }
                }))
            }
        }
    }
}