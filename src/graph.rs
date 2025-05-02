use crate::protein::Protein;
//use crate::FiveMer;

use std::fmt;
use std::rc::Rc;
use std::thread::spawn;
use std::sync::{Arc, mpsc, Weak};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize};
use std::sync::atomic::Ordering::{Acquire, AcqRel, SeqCst};

type ArcGraph = Arc<Graph>;
type ArcU64 = Arc<AtomicU64>;
type ArcVec<T> = Arc<Vec<T>>;
type ArcProtein = Arc<Protein>;
type ArcBool = Arc<AtomicBool>;
type ArcUsize = Arc<AtomicUsize>;
type ArcKmerEdge = Arc<KmerEdge>;
type ArcProteinVertex = Arc<ProteinVertex>;
type ArcKmerEdgeHelper = Arc<KmerEdgeHelper>;

type WeakGraph = Weak<Graph>;

type OptionArcUsize = Rc<Option<ArcUsize>>;

struct KmerEdgeHelper {
    // thread that gets to switch this to false will be allowed to write to KmerEdge
    update_helper: ArcBool,

    // const-sized array containing vertex keys to add to KmerEdge
    array: [ArcUsize; 1000],

    //  1 + 21 + 21 + 21
    //  zero_one switch ->  if zero, will look through first 21
    //                      if one, will look through second 21
    //  third 21 ->         if zero, helper will increment by one while looking through array until it reaches num in first 21
    //                      if one, helper will increment by one while looking through array until it reaches num in second 21
    //  second 21 ->        if zero, threads will increment by one and retrieve previous; if third 21 <= second 21, wait
    //  first 21 ->         if one, threads will increment by one and retrieve previous; if third 21 <= first 21, wait
    filled_indices_tracker: ArcU64,
}

impl KmerEdgeHelper {
    fn new() -> KmerEdgeHelper {
        KmerEdgeHelper { 
            update_helper: Arc::new(AtomicBool::new(false)), 
            array: [0usize; 1000].map(|x| Arc::new(AtomicUsize::new(x))),
            filled_indices_tracker: Arc::new(AtomicU64::new(0u64)),
        }
    }

    fn can_update_kmer_edge(&self) -> bool {
        self.update_helper.swap(false, AcqRel)
    }
}


pub struct KmerEdgeSingle {
    kmer: usize,
    vertices_key: Vec<usize>,
    vertices_bit_array: Vec<ArcBool>,
}
impl KmerEdgeSingle {

    // Make a new KmerEdge
    fn new(kmer: usize, vertices_count: usize) -> KmerEdgeSingle {
        KmerEdgeSingle {
            kmer,
            vertices_key: vec![],
            vertices_bit_array: vec![false; vertices_count]
                .iter().map(|b| Arc::new(AtomicBool::new(*b))).collect(),
        }
    }

    // Change kmer info
    fn update_kmer(&mut self, kmer:usize) {
        self.kmer = kmer;
    }

    // Add key to ProteinVertex
    fn add_vertex(&mut self, vertex_key: usize) {
        if ! self.vertices_bit_array[vertex_key].swap(true, AcqRel) {
            self.vertices_key.push(vertex_key);
        }
    }
}

pub struct KmerEdgeGroup {
    kmers: Vec<usize>,
    graph: WeakGraph,
    vertices_key: Vec<usize>,
    vertices_bit_array: Vec<ArcBool>,
}
impl KmerEdgeGroup {

    // Merge KmerEdge objects into one
    fn new(kmer_edge_keys: &Vec<ArcUsize>, graph: WeakGraph) -> KmerEdgeGroup {
        // Create empty vectors
        let mut kmers: Vec<usize> = vec![];
        let mut vertices_key: Vec<usize> = vec![];
        let mut vertices_bit_array: Vec<ArcBool> = vec![];

        // Make graph data accessible with temporary strong pointer
        let graph_strong = graph.upgrade().unwrap();

        // Append each KmerEdge's kmer and vertices info to aforementioned vectors
        for edge_key in kmer_edge_keys {
            let loaded_edge_key = edge_key.load(Acquire);
            let kmer_edge = &*graph_strong.edges[loaded_edge_key];

            kmers.append(&mut kmer_edge.get_kmers());
            if vertices_bit_array.is_empty() {
                vertices_key.append(&mut kmer_edge.get_vertices_key());
                vertices_bit_array.append(&mut kmer_edge.get_vertices_bit_array());
            }
            else {
                for key in kmer_edge.get_vertices_key() {
                    if !vertices_bit_array[key].swap(true, AcqRel) {
                        vertices_key.push(key);
                    }
                }
            }
        }

        // Remove temporary strong pointer
        drop(graph_strong);

        // Make KmerEdgeGroup object
        KmerEdgeGroup {kmers, graph, vertices_key, vertices_bit_array}
    }

    // Add key to ProteinVertex
    fn add_vertex(&mut self, vertex_key: usize) {
        if ! self.vertices_bit_array[vertex_key].swap(true, AcqRel) {
            self.vertices_key.push(vertex_key);
        }
    }

    // Update proteins to reflect changes in graph
    fn update_protein_vertices(&self,
            new_key: ArcUsize, thread_count: u32) {

        // Make ARC pointers of shared variables 
        let vertices_key: Arc<Vec<usize>> = Arc::new(self.vertices_key.clone()); 
        let vertices_index: ArcUsize = Arc::new(AtomicUsize::new(0));
        let mut handles = vec![];

        for _ in 0..thread_count {
            let vertices_key: Arc<Vec<usize>> = vertices_key.clone();
            let vertices_index: ArcUsize = vertices_index.clone();
            let graph: ArcGraph = self.graph.upgrade().unwrap();
            let new_key: ArcUsize = new_key.clone();
            handles.push(spawn(move || {
                loop {
                    let curr_vertex_index = vertices_index.fetch_add(
                        1, AcqRel
                    );
                    if curr_vertex_index >= vertices_key.len() {
                        break;
                    }
                    let curr_key = vertices_key[curr_vertex_index];
                    let curr_vertex: &mut ArcProteinVertex = 
                        &mut graph.vertices[curr_key].clone();

                    // No two threads will ever traverse the same vertex, so we are safe
                    let curr_vertex_raw: &mut ProteinVertex = unsafe {
                         Arc::get_mut_unchecked(curr_vertex)};
                    curr_vertex_raw.add_edge(new_key.clone());

                }
            }))
        }
        for handle in handles {
            handle.join().unwrap();
        }
    }
}

pub enum KmerEdge {
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

    fn get_vertices_bit_array(&self) -> Vec<ArcBool> {
        match self {
            Self::Single(val) => val.vertices_bit_array.clone(),
            Self::Group(val) => val.vertices_bit_array.clone(),
        }
    }

    fn add_vertex(&mut self, vertex_key: usize) {
        match self {
            Self::Single(val) => val.add_vertex(vertex_key),
            Self::Group(val) => val.add_vertex(vertex_key),
        }
    }

    fn update_kmer(&mut self, kmer: usize) {
        match self {
            Self::Single(val) => val.update_kmer(kmer),
            Self::Group(_) => panic!("Shouldn't be able to call update_kmer"),
        }
    }

    fn update_protein_vertices(&self,
            new_key: ArcUsize, thread_count: u32) {
        match self {
            Self::Single(_) => panic!(
                "Can't call update_protein_vertices on single kmer"),
            Self::Group(val) => val.update_protein_vertices(
                new_key, thread_count),
        }
    }
    fn new(kmer: usize, vertices_count: usize) -> KmerEdge {
        Self::Single(KmerEdgeSingle::new(kmer, vertices_count))
    }
    
    fn merge(kmer_edge_keys: &Vec<ArcUsize>, graph: WeakGraph) -> KmerEdge {
        Self::Group(KmerEdgeGroup::new(kmer_edge_keys, graph))
    }
    // find the relative frequency of vertices in other that are also in self
    fn compare(&self, other: &Self) -> f64 {
        let self_bit_array = self.get_vertices_bit_array();
        let other_keys = other.get_vertices_key();
        let total = other_keys.len();
        let mut abs_freq = 0usize;
        for key in other_keys {
            if self_bit_array[key].load(Acquire) {
                abs_freq += 1;
            }
        }
        (abs_freq as f64) / (total as f64)
    }
}

struct ProteinVertex {
    key: usize,
    graph: WeakGraph,
    edges_key: Vec<ArcUsize>,
}
impl ProteinVertex {

    // Make a new ProteinVertex
    fn new(key: usize, graph: WeakGraph) -> ProteinVertex {
        // Make graph data accessible with temporary strong pointer
        let graph_strong = graph.upgrade().unwrap();

        let protein_obj: ArcProtein = graph_strong.protein_list[key].clone();
        let edges_key: Vec<ArcUsize> = protein_obj.get_five_hash().iter()
            .map(|k: &u32| graph_strong.global_edge_keys[*k as usize].clone())
            .collect();
        // Remove temporary strong pointer
        drop(graph_strong);

        ProteinVertex{key, graph, edges_key}
    }

    // Add key to KmerEdge
    fn add_edge(&mut self, edge_key: ArcUsize) {
        // Make graph data accessible with temporary strong pointer
        let graph_strong = self.graph.upgrade().unwrap();

        let length = graph_strong.edges.len();
        let edges_bit_array: Vec<bool> = self.get_edges_bit_array(length);
        if !edges_bit_array[edge_key.load(Acquire)] {
            self.edges_key.push(edge_key);
        }

        // Remove temporary strong pointer
        drop(graph_strong);
    }

    // Remove keys to KmerEdge
    fn remove_edges(&mut self, removed_bit_array: Arc<Vec<bool>>) {
        let mut index = 0usize;
        while index < self.edges_key.len() {
            if removed_bit_array[self.edges_key[index].load(Acquire)] {
                self.edges_key.remove(index);
            }
            else {
                index += 1;
            }
        }
    }

    // Update ProteinVertex's edges to reflect ProteinVertex's presence
    fn update_graph_edges(&self) {
        // Make graph data accessible with temporary strong pointer
        let graph_strong = self.graph.upgrade().unwrap();

        for key in &self.edges_key {
            let loaded_key = key.load(Acquire);

            let curr_edge: &mut ArcKmerEdge = 
                &mut graph_strong.edges[loaded_key].clone();

            // Two threads may try to access the same edge -> how do we deal with this?
            let curr_edge_raw: &mut KmerEdge = unsafe {
                Arc::get_mut_unchecked(curr_edge)};
            curr_edge_raw.add_vertex(self.key);
        }

        // Remove temporary strong pointer
        drop(graph_strong);
    }

    // Get bit array of ProteinVertex's edges
    fn get_edges_bit_array(&self, length: usize) -> Vec<bool> {
        let mut edges_bit_array: Vec<bool> = vec![false; length];
        for key in &self.edges_key {
            edges_bit_array[key.load(Acquire)] = true;
        }
        edges_bit_array
    }
}

pub struct Graph {
    pub edges: ArcVec<ArcKmerEdge>,
    protein_list: ArcVec<ArcProtein>,
    vertices: ArcVec<ArcProteinVertex>,
    global_edge_keys: ArcVec<ArcUsize>,
    edge_update_helpers : ArcVec<ArcKmerEdgeHelper>,
}
impl Graph {

    // Create a graph and return an ARC pointer to said graph
    pub fn new(kmer_count: usize, thread_count: u32, 
            protein_list: ArcVec<ArcProtein>) -> ArcGraph {

        // Create Graph object skeleton
        let global_edge_keys: ArcVec<ArcUsize> = Arc::new((0..kmer_count).into_iter()
            .map(|k: usize| Arc::new(AtomicUsize::new(k))).collect());
        let vertices_count: usize = protein_list.len();

        // Create empty edges
        let edges: ArcVec<ArcKmerEdge> = Arc::new(Vec::new());
        let edges_index: ArcUsize = Arc::new(AtomicUsize::new(0usize));

        let arc_kmer_count: ArcUsize = Arc::new(AtomicUsize::new(kmer_count));
        let arc_vertices_count: ArcUsize = Arc::new(AtomicUsize::new(vertices_count));

        let mut handles  = vec![];
        for _ in 0..thread_count {
            let mut edges: ArcVec<ArcKmerEdge> = edges.clone();
            let edges_index: ArcUsize = edges_index.clone();
            let arc_kmer_count: ArcUsize = arc_kmer_count.clone();
            let arc_vertices_count: ArcUsize = arc_vertices_count.clone();
            handles.push(spawn(move || {
                loop {
                    let curr_edge_index = edges_index.load(Acquire);
                    let kmer_count = arc_kmer_count.load(Acquire);
                    if curr_edge_index >= kmer_count {
                        break;
                    }
                    let mut curr_edge: ArcKmerEdge = Arc::new(KmerEdge::new(
                        0, arc_vertices_count.load(Acquire)));
                    let curr_edge_index: Result<usize, usize> = edges_index
                        .fetch_update(SeqCst, SeqCst, |x| {
                            if x >= kmer_count {
                                None
                            }
                            else {
                                // Part of atomic operation, so we are safe
                                unsafe { Arc::get_mut_unchecked(&mut edges)
                                    .push(curr_edge.clone()); }
                                Some(x+1)
                            }
                        });
                    if curr_edge_index.is_err() {
                        break;
                    }
                    // Only current thread will access this edge, so we are safe
                    unsafe { Arc::get_mut_unchecked(&mut curr_edge)
                        .update_kmer(curr_edge_index.unwrap());
                    }
                }
                
            }));
        }

        // Create edge update helpers
        let edge_update_helpers: ArcVec<ArcKmerEdgeHelper> = Arc::new(Vec::new());
        let edge_helpers_index: ArcUsize = Arc::new(AtomicUsize::new(0usize));
        
        let arc_kmer_count: ArcUsize = Arc::new(AtomicUsize::new(kmer_count));

        let mut handles  = vec![];
        for _ in 0..thread_count {
            let mut edge_helpers: ArcVec<ArcKmerEdgeHelper> = edge_update_helpers.clone();
            let edge_helpers_index: ArcUsize = edge_helpers_index.clone();
            let arc_kmer_count: ArcUsize = arc_kmer_count.clone();
            handles.push(spawn(move || {
                loop {
                    let curr_helper_index = edge_helpers_index.load(Acquire);
                    let kmer_count = arc_kmer_count.load(Acquire);
                    if curr_helper_index >= kmer_count {
                        break;
                    }
                    let curr_helper: ArcKmerEdgeHelper = Arc::new(KmerEdgeHelper::new());
                    let curr_helper_index: Result<usize, usize> = edges_index
                        .fetch_update(SeqCst, SeqCst, |x| {
                            if x >= kmer_count {
                                None
                            }
                            else {
                                // Part of atomic operation, so we are safe
                                unsafe { Arc::get_mut_unchecked(&mut edge_helpers)
                                    .push(curr_helper.clone()); }
                                Some(x+1)
                            }
                        });
                    if curr_helper_index.is_err() {
                        break;
                    }
                }
            }));
        }

        let graph: Graph = Graph {
            vertices: Arc::new(Vec::new()),
            edges,
            protein_list,
            global_edge_keys,
            edge_update_helpers,
        };
        
        let mut graph_arc: ArcGraph = Arc::new(graph);

        // Make vertices
        let vertices: ArcVec<ArcProteinVertex> = Arc::new(
            (0..vertices_count).into_iter().map(|k: usize| Arc::new(
                ProteinVertex::new(k, Arc::downgrade(&graph_arc))
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
                    let curr_vertices_index = vertices_index.fetch_add(1, AcqRel);
                    if curr_vertices_index >= 10619 {
                        break;
                    }

                    // CONCURRENCY
                    vertices[curr_vertices_index].update_graph_edges();
                }
            }))
        }
        for handle in handles {
            handle.join().unwrap();
        }
        // Currently in a single threaded environment, so we are safe
        unsafe {Arc::get_mut_unchecked(&mut graph_arc).vertices = vertices};
        graph_arc
    }

    // Combine edges with lots of overlap
    pub fn combine_edges(mut graph: ArcGraph, thread_count: u32) {
        // Get current keys to edges
        let mut unordered_edge_keys: Vec<OptionArcUsize> = (0..graph.edges.len())
            .into_iter().map(|i| Rc::new(Some(graph.global_edge_keys[i].clone())))
            .collect();

        // Make a second vector of current keys, ordered by edge size
        let mut ordered_edge_keys: Vec<OptionArcUsize> = unordered_edge_keys.clone();
        ordered_edge_keys.sort_by(|a: &OptionArcUsize, b: &OptionArcUsize| {
            let cloned_a: Option<ArcUsize> = a.as_ref().clone();
            let cloned_b: Option<ArcUsize> = b.as_ref().clone();
            let key_a: usize = cloned_a.unwrap().load(Acquire);
            let key_b: usize = cloned_b.unwrap().load(Acquire);
            let edge_a: ArcKmerEdge = graph.edges[key_a].clone();
            let edge_b: ArcKmerEdge = graph.edges[key_b].clone();
            let edge_a_size: usize = edge_a.get_vertices_key().len();
            let edge_b_size: usize = edge_b.get_vertices_key().len();
            edge_b_size.cmp(&edge_a_size)});
        
        // Traverse through each edge
        for edge_key in ordered_edge_keys {
            // If edge was merged, skip
            if edge_key.is_none() {
                continue;
            }

            // Get actual pointer to KmerEdge
            let unwraped_edge_key: ArcUsize = edge_key.as_ref().clone().unwrap();
            let loaded_edge_key: usize = unwraped_edge_key.load(Acquire);
            let edge: ArcKmerEdge = graph.edges[loaded_edge_key].clone();
            eprintln!("Currently at edge {loaded_edge_key} of size {}", 
                edge.get_vertices_key().len());

            // Make ARC pointers of shared variables   
            let vertices_index: ArcUsize = 
                Arc::new(AtomicUsize::new(1));
            let vertices_keys: ArcVec<usize> = 
                Arc::new(edge.get_vertices_key());
            let first_vertex: ArcProteinVertex = graph.vertices[vertices_keys[0]].clone();
            let other_edges_bitarray: ArcVec<AtomicBool> = 
                Arc::new(first_vertex.get_edges_bit_array(graph.edges.len())
                    .iter().map(|b: &bool| AtomicBool::new(*b)).collect());
            let mut other_edges_keys: Vec<ArcUsize> = 
                first_vertex.edges_key.clone();

            // Make channels and thread handles
            let (send, recv) = mpsc::channel::<ArcUsize>();
            let mut handles  = vec![];

            // Retrieve other edges that share a vertex with current edge
            for _ in 0..thread_count {
                let graph: ArcGraph = graph.clone();
                let send: mpsc::Sender<ArcUsize> = send.clone();
                let vertices_index: ArcUsize = vertices_index.clone();
                let vertices_keys: ArcVec<usize> = vertices_keys.clone();
                let other_edges_bitarray: ArcVec<AtomicBool> = other_edges_bitarray.clone();
                
                handles.push(spawn(move || {
                    loop {
                        let curr_vertices_index = vertices_index.fetch_add(
                            1, AcqRel
                        );
                        if curr_vertices_index >= vertices_keys.len() {
                            break;
                        }

                        let curr_vertex_key = vertices_keys[curr_vertices_index];
                        let curr_vertex: &ArcProteinVertex = 
                            &graph.vertices[curr_vertex_key];

                        let more_edges_keys: &Vec<ArcUsize> = &curr_vertex.edges_key;

                        for key in more_edges_keys {
                            let loaded_key = key.load(Acquire);
                            if !other_edges_bitarray[loaded_key].swap(true, AcqRel) {
                                send.send(key.clone()).unwrap();
                            }
                        }
                    }
                }))
            }
            drop(send);
            let mut key = recv.recv();
            while key.is_ok() {
                other_edges_keys.push(key.unwrap());
                key = recv.recv();
            }

            for handle in handles {
                handle.join().unwrap();
            }
            drop(vertices_index);
            drop(other_edges_bitarray);

            // Gather keys to other edges
            if other_edges_keys.len() == 1 {
                continue;
            }

            // Sort keys
            other_edges_keys.sort_by(|a: &ArcUsize, b: &ArcUsize| 
                a.load(Acquire).cmp(&b.load(Acquire)));

            // Make ARC pointers of shared variables 
            let other_edges_index: ArcUsize = 
                Arc::new(AtomicUsize::new(0));
            let other_edges_keys: ArcVec<ArcUsize> = 
                Arc::new(other_edges_keys.clone());
            let mut edges_to_merge_with: Vec<ArcUsize> = Vec::new();

            // Make channels and thread handles
            let (send, recv) = mpsc::channel::<ArcUsize>();
            let mut handles= vec![];

            // Find which of the other edges can be merged with current edge
            for _ in 0..thread_count {
                let graph: ArcGraph = graph.clone();
                let edge: ArcKmerEdge = edge.clone();
                let send: mpsc::Sender<ArcUsize> = send.clone();
                let other_edges_index: ArcUsize = other_edges_index.clone();
                let other_edges_keys: ArcVec<ArcUsize> = other_edges_keys.clone();

                handles.push(spawn(move || {
                    loop {
                        let other_edges_index = other_edges_index.fetch_add(
                            1, AcqRel);
                        if other_edges_index >= other_edges_keys.len() {
                            break;
                        }

                        let other_key: &ArcUsize = &other_edges_keys[other_edges_index];
                        let loaded_other_key: usize = other_key.load(Acquire);
                        let other_edge: &ArcKmerEdge = &graph.edges[loaded_other_key];

                        let other_cmp_edge: f64 = edge.compare(other_edge);
                        let edge_cmp_other: f64 = other_edge.compare(&edge);
                        if other_cmp_edge == 1.0 {
                            send.send(other_key.clone()).unwrap();
                        }
                        else if edge_cmp_other >= 0.8 && other_cmp_edge >= 0.8 {
                            send.send(other_key.clone()).unwrap();
                        }
                    }
                }))
            }
            drop(send);
            let mut key = recv.recv();
            while key.is_ok() {
                edges_to_merge_with.push(key.unwrap());
                key = recv.recv();
            }

            for handle in handles {
                handle.join().unwrap();
            }
            drop(other_edges_keys);
            drop(other_edges_index);
            if edges_to_merge_with.len() == 1 {
                continue;
            }

            // Merge edges into one
            let mut kmer_edge_group: KmerEdge = KmerEdge::merge(
                &edges_to_merge_with, Arc::downgrade(&graph));
            kmer_edge_group.update_protein_vertices(unwraped_edge_key, thread_count);

            // Currently in a single threaded environment, so we are safe
            unsafe {
                let graph_mut = Arc::get_mut_unchecked(&mut graph);
                Arc::get_mut_unchecked(&mut graph_mut.edges)[loaded_edge_key] = 
                    Arc::new(kmer_edge_group)
            };
            
            for old_key in &*edges_to_merge_with {
                let loaded = old_key.load(Acquire);
                if loaded != loaded_edge_key {
                    unsafe {
                        *Rc::get_mut_unchecked(&mut unordered_edge_keys[loaded]) = None;
                    }
                }
            }
        }

        // Remove edges that have been merged with a bigger one
        let removed_bit_array = Self::remove_edges_marked_for_deletions(
            graph.clone(), unordered_edge_keys);

        // Remove weak pointers to non-existing edges
        let mut handles= vec![];
        let vertices_index: ArcUsize = 
            Arc::new(AtomicUsize::new(0));
        let removed_bit_array: ArcVec<bool> = Arc::new(removed_bit_array);

        for _ in 0..thread_count {
            let vertices_index: ArcUsize = vertices_index.clone();
            let removed_bit_array: ArcVec<bool> = removed_bit_array.clone();
            let mut vertices: ArcVec<ArcProteinVertex> = Arc::new((*graph.vertices).clone());

            handles.push(spawn(move || {
                loop {
                    let curr_vertices_index = vertices_index.fetch_add(1, AcqRel);
                    if curr_vertices_index >= 10619 {
                        break;
                    }

                    // No two threads will ever traverse the same vertex, so we are safe
                    unsafe { 
                        let vertices_mut: &mut Vec<Arc<ProteinVertex>>  = 
                            Arc::get_mut_unchecked(&mut vertices);
                        Arc::get_mut_unchecked(&mut vertices_mut[curr_vertices_index])
                            .remove_edges(removed_bit_array.clone()) 
                    };
                }
            }))
        }
    }

    fn remove_edges_marked_for_deletions(
            mut graph: ArcGraph, unordered_edge_keys: Vec<OptionArcUsize>) -> Vec<bool>{

        let edge_num = graph.edges.len();
        let mut removed_bit_array = vec![false; edge_num];
        let mut substractive = 0usize;
        for (index, edge_key) in unordered_edge_keys.iter().enumerate() {
            if edge_key.is_none() {
                // Currently in a single threaded environment, so we are safe
                unsafe {
                    let graph_mut = Arc::get_mut_unchecked(&mut graph);
                    Arc::get_mut_unchecked(&mut graph_mut.edges)
                        .remove(index - substractive);
                    Arc::get_mut_unchecked(&mut graph_mut.global_edge_keys)
                        .remove(index - substractive);
                }
                substractive += 1;
                removed_bit_array[index] = true;
            }
            else if substractive > 0 {
                // Currently in a single threaded environment, so we are safe
                unsafe {
                    let graph_mut = Arc::get_mut_unchecked(&mut graph);
                    let edge_keys = 
                        Arc::get_mut_unchecked(&mut graph_mut.global_edge_keys);
                    let edge_key = 
                        Arc::get_mut_unchecked(&mut edge_keys[index-substractive]);
                    edge_key.fetch_sub(substractive, AcqRel);
                }
            }
        }
        removed_bit_array
    }

    // Remove edges without diverging AMR labels
    // pub fn remove_uninteresting_edges(mut graph: ArcGraph, thread_count: u32) {
    //     // Get current keys to edges
    //     let mut unordered_edge_keys: Vec<OptionArcUsize> = (0..graph.edges.len())
    //         .into_iter().map(|i| Rc::new(Some(graph.global_edge_keys[i].clone())))
    //         .collect();
    // }
}



impl fmt::Debug for Graph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Graph")
            .field("Kmers", &self.edges)
            .field("Proteins", &self.vertices).finish()
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