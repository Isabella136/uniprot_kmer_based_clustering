use crate::protein::Protein;

use std::fmt;
use std::rc::Rc;
use std::ops::Deref;
use std::thread::spawn;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, MutexGuard, RwLock, RwLockReadGuard};


type ArcGraph = Arc<Graph>;
type ArcVec<T> = Arc<Vec<T>>;
type ArcProtein = Arc<Protein>;
type ArcMutex<T> = Arc<Mutex<T>>;
type ArcUsize = Arc<AtomicUsize>;
type ArcKmerEdge = Arc<RwLock<KmerEdge>>;
type ArcProteinVertex = Arc<RwLock<ProteinVertex>>;

type ReadGuard<'a, T> = RwLockReadGuard<'a, T>;

struct KmerEdgeSingle {
    kmer: usize,
    graph: ArcGraph,
    vertices_key: Vec<usize>,
    vertices_bit_array: Vec<bool>,
}
impl KmerEdgeSingle {

    // Make a new KmerEdge
    fn new(graph: ArcGraph, kmer: usize, vertices_count: usize) -> KmerEdgeSingle {
        KmerEdgeSingle {
            kmer,
            graph,
            vertices_key: vec![],
            vertices_bit_array: vec![false; vertices_count],
        }
    }

    // Add key to ProteinVertex
    fn add_vertex(&mut self, vertex_key: usize) {
        if ! self.vertices_bit_array[vertex_key] {
            self.vertices_bit_array[vertex_key] = true;
            self.vertices_key.push(vertex_key);
        }
    }
}

struct KmerEdgeGroup {
    kmers: Vec<usize>,
    graph: ArcGraph,
    vertices_key: Vec<usize>,
    vertices_bit_array: Vec<bool>,
}
impl KmerEdgeGroup {

    // Merge KmerEdge objects into one
    fn new(kmer_edge_keys: &Vec<ArcUsize>, graph: ArcGraph) -> KmerEdgeGroup {
        // Create empty vectors
        let mut kmers: Vec<usize> = vec![];
        let mut vertices_key: Vec<usize> = vec![];
        let mut vertices_bit_array: Vec<bool> = vec![];

        // Append each KmerEdge's kmer and vertices info to aforementioned vectors
        for edge_key in kmer_edge_keys {
            let loaded_edge_key = edge_key.load(Ordering::Acquire);
            let kmer_edge = graph.edges[loaded_edge_key].read().unwrap();

            kmers.append(&mut kmer_edge.get_kmers());
            vertices_key.append(&mut kmer_edge.get_vertices_key());
            vertices_bit_array.append(&mut kmer_edge.get_vertices_bit_array());
        }

        // Make KmerEdgeGroup object
        KmerEdgeGroup {kmers, graph, vertices_key, vertices_bit_array}
    }

    // Add key to ProteinVertex
    fn add_vertex(&mut self, vertex_key: usize) {
        if ! self.vertices_bit_array[vertex_key] {
            self.vertices_bit_array[vertex_key] = true;
            self.vertices_key.push(vertex_key);
        }
    }

    // Update proteins to reflect changed in graph
    fn update_protein_vertices(&self,
            old_keys: &Vec<ArcUsize>, new_key: ArcUsize, 
            total_edges: usize, thread_count: u32) {

        // Define bit array of edge keys to remove    
        let mut removed_bit_array: Vec<bool> = vec![false; total_edges];
        for key in old_keys {
            removed_bit_array[key.load(Ordering::Acquire)] = true;
        }
        removed_bit_array[new_key.load(Ordering::Acquire)] = false;

        // Make ARC pointers of shared variables        
        let vertices_key: Arc<&Vec<usize>> = Arc::new(&self.vertices_key);
        let vertices_index: ArcUsize = Arc::new(AtomicUsize::new(0));
        let old_keys: Arc<&Vec<ArcUsize>> = Arc::new(old_keys);
        let mut handles = vec![];

        for _ in 0..thread_count {
            let vertices_key: Arc<&Vec<usize>> = vertices_key.clone();
            let vertices_index: ArcUsize = vertices_index.clone();
            let old_keys: Arc<&Vec<ArcUsize>> = old_keys.clone();
            let graph: ArcGraph = self.graph.clone();
            handles.push(spawn(move || {
                loop {
                    
                }
            }))

        }
        for handle in handles {
            handle.join().unwrap();
        }
    }
}

enum KmerEdge {
    Single(KmerEdgeSingle),
    Group(KmerEdgeGroup),
}
impl KmerEdge {
    fn get_kmers(&self) -> Vec<usize> {
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
    fn get_graph(&self) -> ArcGraph {
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
    fn new(graph: ArcGraph, kmer: usize, vertices_count: usize) -> KmerEdge {
        Self::Single(KmerEdgeSingle::new(graph, kmer, vertices_count))
    }
    // find the relative frequency of vertices in other that are also in self
    fn merge(kmer_edge_keys: &Vec<ArcUsize>, graph: ArcGraph) -> KmerEdge {
        Self::Group(KmerEdgeGroup::new(kmer_edge_keys, graph))
    }
    fn compare(&self, other: &Self) -> f64 {
        let self_bit_array = self.get_vertices_bit_array();
        let other_keys = other.get_vertices_key();
        let total = other_keys.len();
        let mut abs_freq = 0usize;
        for key in other_keys {
            if self_bit_array[key] {
                abs_freq += 1;
            }
        }
        (abs_freq as f64) / (total as f64)
    }
}

struct ProteinVertex {
    key: usize,
    graph: ArcGraph,
    edges_key: Vec<ArcUsize>,
}
impl ProteinVertex {

    // Make a new ProteinVertex
    fn new(key: usize, graph: ArcGraph) -> ProteinVertex {
        let protein_obj: ArcProtein = graph.protein_list[key].clone();
        let edges_key: Vec<ArcUsize> = protein_obj.get_five_hash().iter()
            .map(|k: &u32| graph.global_edge_keys[*k as usize].clone()).collect();
        ProteinVertex{key, graph, edges_key}
    }

    // Update ProteinVertex's edges to reflect ProteinVertex's presence
    fn update_graph_edges(&self) {
        for key in &self.edges_key {
            let loaded_key = key.load(Ordering::Acquire);
            let mut edge = self.graph.edges[loaded_key].write().unwrap();
            edge.add_vertex(self.key);
        }
    }

    // Get bit array of ProteinVertex's edges
    fn get_edges_bit_array(&self, length: usize) -> Vec<bool> {
        let mut edges_bit_array: Vec<bool> = vec![false; length];
        for key in &self.edges_key {
            edges_bit_array[key.load(Ordering::Acquire)] = true;
        }
        edges_bit_array
    }
}

pub struct Graph {
    edges: ArcVec<ArcKmerEdge>,
    protein_list: ArcVec<ArcProtein>,
    vertices: ArcVec<ArcProteinVertex>,
    global_edge_keys: ArcVec<ArcUsize>,
}
impl Graph {

    // Create a graph and return an ARC pointer to said graph
    pub fn new(kmer_count: usize, thread_count: u32, 
            protein_list: ArcVec<ArcProtein>) -> ArcGraph {

        // Create Graph object skeleton
        let global_edge_keys: ArcVec<ArcUsize> = Arc::new((0..kmer_count).into_iter()
            .map(|k: usize| Arc::new(AtomicUsize::new(k))).collect());
        let vertices_count: usize = protein_list.len();
        let graph: Graph = Graph {
            vertices: Arc::new(Vec::new()),
            edges: Arc::new(Vec::new()),
            protein_list,
            global_edge_keys,
        };

        // Add empty edges
        let mut graph_arc: ArcGraph = Arc::new(graph);
        let edges: ArcVec<ArcKmerEdge> = Arc::new(
            (0..kmer_count).into_iter().map(|kmer: usize| Arc::new(RwLock::new(
                KmerEdge::new(graph_arc.clone(), kmer.clone(), vertices_count))
            )).collect()
        );
        unsafe {Arc::get_mut_unchecked(&mut graph_arc).edges = edges};
        
        // Make vertices
        let vertices: ArcVec<ArcProteinVertex> = Arc::new(
            (0..vertices_count).into_iter().map(|k: usize| Arc::new(RwLock::new(
                ProteinVertex::new(k, graph_arc.clone()))
            )).collect()
        );
        let vertices_index: ArcUsize = Arc::new(AtomicUsize::new(0usize));
        let mut handles  = vec![];

        // Fill up edges with vertices
        for _ in 0..thread_count {
            let vertices: ArcVec<ArcProteinVertex> = vertices.clone();
            let vertices_index: ArcUsize = vertices_index.clone();

            handles.push(spawn(move || {
                loop {
                    let curr_vertices_index = vertices_index.fetch_add(1, Ordering::AcqRel);
                    if curr_vertices_index >= 10619 {
                        break;
                    }
                    vertices[curr_vertices_index].read().unwrap().update_graph_edges();
                }
            }))
        }
        for handle in handles {
            handle.join().unwrap();
        }
        unsafe {Arc::get_mut_unchecked(&mut graph_arc).vertices = vertices};
        graph_arc
    }

    // Combine edges with lots of overlap
    pub fn combine_edges(mut graph: ArcGraph, thread_count: u32) {

        type OptionArcUsize = Rc<Option<ArcUsize>>;

        // Get current keys to edges
        let mut unordered_edge_keys: Vec<OptionArcUsize> = (0..graph.edges.len())
            .into_iter().map(|i| Rc::new(Some(graph.global_edge_keys[i].clone())))
            .collect();

        // Make a second vector of current keys, ordered by edge size
        let mut ordered_edge_keys: Vec<OptionArcUsize> = unordered_edge_keys.clone();
        ordered_edge_keys.sort_by(|a: &OptionArcUsize, b: &OptionArcUsize| {
            let cloned_a: Option<ArcUsize> = a.as_ref().clone();
            let cloned_b: Option<ArcUsize> = b.as_ref().clone();
            let key_a: usize = cloned_a.unwrap().load(Ordering::Acquire);
            let key_b: usize = cloned_b.unwrap().load(Ordering::Acquire);
            let edge_a: ArcKmerEdge = graph.edges[key_a].clone();
            let edge_b: ArcKmerEdge = graph.edges[key_b].clone();
            let edge_a_size: usize = edge_a.read().unwrap().get_vertices_key().len();
            let edge_b_size: usize = edge_b.read().unwrap().get_vertices_key().len();
            edge_a_size.cmp(&edge_b_size)});
        
        // Traverse through each edge
        for edge_key in ordered_edge_keys {
            // If edge was merged, skip
            if edge_key.is_none() {
                continue;
            }

            // Get actual pointer to KmerEdge
            let unwraped_edge_key: ArcUsize = edge_key.as_ref().clone().unwrap();
            let loaded_edge_key: usize = unwraped_edge_key.load(Ordering::Acquire);
            let edge: ArcKmerEdge = graph.edges[loaded_edge_key].clone();

            // Make ARC pointers of shared variables   
            let vertices_index: ArcUsize = 
                Arc::new(AtomicUsize::new(1));
            let vertices_keys: ArcVec<usize> = 
                Arc::new(edge.read().unwrap().get_vertices_key());
            let first_vertex: ReadGuard<ProteinVertex> = 
                graph.vertices[vertices_keys[0]].read().unwrap();
            let other_edges_bitarray: ArcVec<AtomicBool> = 
                Arc::new(first_vertex.get_edges_bit_array(graph.edges.len())
                    .iter().map(|b: &bool| AtomicBool::new(*b)).collect());
            let other_edges_keys: ArcMutex<Vec<ArcUsize>> = 
                Arc::new(Mutex::new(first_vertex.edges_key.clone()));
            let mut handles  = vec![];

            // Retrieve other edges that share a vertex with current edge
            for _ in 0..thread_count {
                let graph: ArcGraph = graph.clone();
                let vertices_index: ArcUsize = vertices_index.clone();
                let vertices_keys: ArcVec<usize> = vertices_keys.clone();
                let other_edges_keys: ArcMutex<Vec<ArcUsize>> = other_edges_keys.clone();
                let other_edges_bitarray: ArcVec<AtomicBool> = other_edges_bitarray.clone();
                
                handles.push(spawn(move || {
                    loop {
                        let curr_vertices_index = vertices_index.fetch_add(
                            1, Ordering::AcqRel
                        );
                        if curr_vertices_index >= vertices_keys.len() {
                            break;
                        }

                        let curr_vertex_key = vertices_keys[curr_vertices_index];
                        let curr_vertex: &ArcProteinVertex = 
                            &graph.vertices[curr_vertex_key];
                        let more_edges_keys: &Vec<ArcUsize> = 
                            &curr_vertex.read().unwrap().edges_key;

                        for key in more_edges_keys {
                            let loaded_key = key.load(Ordering::Acquire);
                            if !other_edges_bitarray[loaded_key].swap(true, Ordering::AcqRel) {
                                other_edges_keys.lock().unwrap().push(key.clone());
                            }
                        }
                    }
                }))
            }
            for handle in handles {
                handle.join().unwrap();
            }
            drop(vertices_index);
            drop(other_edges_bitarray);

            // Gather keys to other edges
            let mut other_edges_keys: MutexGuard<Vec<ArcUsize>> = 
                other_edges_keys.lock().unwrap();
            if other_edges_keys.len() == 1 {
                continue;
            }

            // Sort keys
            other_edges_keys.sort_by(|a: &ArcUsize, b: &ArcUsize| 
                a.load(Ordering::Acquire).cmp(&b.load(Ordering::Acquire)));

            // Make ARC pointers of shared variables 
            let other_edges_index: ArcUsize = 
                Arc::new(AtomicUsize::new(0));
            let other_edges_keys: ArcVec<ArcUsize> = 
                Arc::new(other_edges_keys.clone());
            let edges_to_merge_with: ArcMutex<Vec<ArcUsize>> = 
                Arc::new(Mutex::new(Vec::new()));
            let mut handles= vec![];

            // Find which of the other edges can be merged with current edge
            for _ in 0..thread_count {
                let graph: ArcGraph = graph.clone();
                let edge: ArcKmerEdge = edge.clone();
                let other_edges_index: ArcUsize = other_edges_index.clone();
                let other_edges_keys: ArcVec<ArcUsize> = other_edges_keys.clone();
                let edges_to_merge_with: ArcMutex<Vec<ArcUsize>> = edges_to_merge_with.clone();

                handles.push(spawn(move || {
                    loop {
                        let other_edges_index = other_edges_index.fetch_add(
                            1, Ordering::AcqRel);
                        if other_edges_index >= other_edges_keys.len() {
                            break;
                        }

                        let other_key: &ArcUsize = &other_edges_keys[other_edges_index];
                        let loaded_other_key: usize = other_key.load(Ordering::Acquire);
                        let other_edge: &ArcKmerEdge = &graph.edges[loaded_other_key];

                        let other_cmp_edge: f64 = edge.read().unwrap()
                            .compare(&*other_edge.read().unwrap());
                        let edge_cmp_other: f64 = other_edge.read().unwrap()
                            .compare(&*edge.read().unwrap());
                        if other_cmp_edge == 1.0 {
                            edges_to_merge_with.lock().unwrap().push(other_key.clone());
                        }
                        else if edge_cmp_other == 1.0 {
                            edges_to_merge_with.lock().unwrap().push(other_key.clone());
                        }
                        else if edge_cmp_other >= 0.8 && other_cmp_edge >= 0.8 {
                            edges_to_merge_with.lock().unwrap().push(other_key.clone());
                        }
                    }
                }))
            }
            for handle in handles {
                handle.join().unwrap();
            }
            drop(other_edges_keys);
            drop(other_edges_index);

            // Gather keys to edges we're merging with
            let edges_to_merge_with: MutexGuard<Vec<ArcUsize>> = 
                edges_to_merge_with.lock().unwrap();
            if edges_to_merge_with.len() == 1 {
                continue;
            }

            // Merge edges into one
            let kmer_edge_group: KmerEdge = KmerEdge::merge(
                edges_to_merge_with.deref(), graph.clone());
        }
    }
}



impl fmt::Debug for Graph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.edges.iter()).finish()
    }
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