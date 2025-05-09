use crate::Graph;
use crate::graph::vertex::ProteinVertex;

use std::fmt;
use std::thread::spawn;
use std::sync::{Arc, Weak};
use std::sync::atomic::Ordering::{Acquire, AcqRel, SeqCst, Release};
use std::sync::atomic::{AtomicBool, AtomicPtr, AtomicU64, AtomicUsize};

type ArcU64 = Arc<AtomicU64>;
type ArcBool = Arc<AtomicBool>;
type ArcUsize = Arc<AtomicUsize>;
type ArcProteinVertex = Arc<AtomicPtr<ProteinVertex>>;

type WeakUsize = Weak<AtomicUsize>;
type WeakGraph = Weak<AtomicPtr<Graph>>;

const ONE_SWITCH: u64 = 
    0b1000000000000000000000000000000000000000000000000000000000000000;
const INCR_FIRST: u64 = 
    0b0000000000000000000001000000000000000000000000000000000000000000;
const FIRST_BLOCK: u64 = 
    0b0111111111111111111111000000000000000000000000000000000000000000;
const FIRST_EMPTY: u64 =
    0b1000000000000000000000111111111111111111111111111111111111111111;
const INCR_SECOND: u64 = 
    0b0000000000000000000000000000000000000000001000000000000000000000;
const SECOND_BLOCK: u64 = 
    0b0000000000000000000000111111111111111111111000000000000000000000;
const SECOND_EMPTY: u64 =
    0b1111111111111111111111000000000000000000000111111111111111111111;
const INCR_THIRD: u64 = 
    0b0000000000000000000000000000000000000000000000000000000000000001;

pub struct KmerEdgeHelper {
    // thread that gets to switch this from true to false can update KmerEdge
    update_helper: ArcBool,

    // const-sized arrays containing vertex keys to add to KmerEdge
    array_0: Vec<ArcUsize>,
    array_1: Vec<ArcUsize>,

    // bit-array informing whether other threads have populated vertex key array
    bit_array_0: Vec<ArcBool>,
    bit_array_1: Vec<ArcBool>,

    //  1 + 21 + 21 + 21
    //  zero_one switch ->  if zero, will look through first 21
    //                      if one, will look through second 21
    //  third 21 ->         number of times that tracker was updated (avoids ABA problem)
    //  second 21 ->        if zero, threads will increment by one and retrieve previous
    //  first 21 ->         if one, threads will increment by one and retrieve previous
    filled_indices_tracker: ArcU64,
}

impl KmerEdgeHelper {
    pub fn new() -> KmerEdgeHelper {
        KmerEdgeHelper { 
            update_helper: Arc::new(AtomicBool::new(true)), 
            array_0: vec![0usize; 4000].into_iter()
                .map(|x| Arc::new(AtomicUsize::new(x))).collect(),
            array_1: vec![0usize; 4000].into_iter()
                .map(|x| Arc::new(AtomicUsize::new(x))).collect(),
            bit_array_0: vec![false; 4000].into_iter()
                .map(|x| Arc::new(AtomicBool::new(x))).collect(),
            bit_array_1: vec![false; 4000].into_iter()
                .map(|x| Arc::new(AtomicBool::new(x))).collect(),
            filled_indices_tracker: Arc::new(AtomicU64::new(0u64)),
        }
    }

    // If returns true, then current thread is allowed to update KmerEdge
    pub fn can_update_kmer_edge(&self) -> bool {
        self.update_helper.swap(false, AcqRel)
    }

    // Add all keys in array as well as key in argument list to KmerEdge
    pub fn update_kmer_edge(&self, key: usize, edge: &mut KmerEdge) {
        let curr_index_info = self.filled_indices_tracker.fetch_update(
            SeqCst, SeqCst, |x| {
                if x & ONE_SWITCH > 0 {
                    Some((x - ONE_SWITCH + INCR_THIRD) & SECOND_EMPTY)
                }
                else {
                    Some((x + ONE_SWITCH + INCR_THIRD) & FIRST_EMPTY)
                }
        }).unwrap();

        // Add key in argument list to array
        let (curr_switch, curr_index) = Self::parse_index_info(curr_index_info);
        self.push(key, &curr_switch, &curr_index);

        // In case there are threads that are still running, come back to them later
        let mut leftovers = vec![];

        // Add keys in array 
        for i in 0..(curr_index+1) {
            if !curr_switch {
                if !self.bit_array_0[i].load(Acquire) {
                    leftovers.push(i);
                }
                else {
                    edge.add_vertex(self.array_0[i].load(Acquire));
                    self.bit_array_0[i].swap(false, AcqRel);
                }
            }
            else {
                if !self.bit_array_1[i].load(Acquire) {
                    leftovers.push(i);
                }
                else {
                    edge.add_vertex(self.array_1[i].load(Acquire));
                    self.bit_array_1[i].store(false, Release);
                }
            }
        }

        // Add keys in array if not previously added
        for i in leftovers {
            if !curr_switch {
                while !self.bit_array_0[i].load(Acquire) {
                    continue;
                }
                edge.add_vertex(self.array_0[i].load(Acquire));
                self.bit_array_0[i].store(false, Release);
            }
            else {
                while !self.bit_array_1[i].load(Acquire) {
                    continue;
                }
                edge.add_vertex(self.array_1[i].load(Acquire));
                self.bit_array_1[i].store(false, Release);
            }
        }

        // Now another thread can update KmerEdge if needed
        self.update_helper.store(true, Release);
    }

    // During last run-through, make a final update to KmerEdge if some keys are leftover
    pub fn update_kmer_edge_final(&self, edge: &mut KmerEdge) {
        let (curr_switch, curr_index) = Self::parse_index_info(
            self.filled_indices_tracker.load(Acquire));

        for i in 0..curr_index {
            if !curr_switch {
                edge.add_vertex(self.array_0[i].load(Acquire));
                self.bit_array_0[i].store(false, Release);
            }
            else {
                edge.add_vertex(self.array_1[i].load(Acquire));
                self.bit_array_1[i].store(false, Release);
            }
        }
        self.update_helper.store(true, Release);
    }

    // Get switch orientation and current index of updating array
    fn parse_index_info(curr_index_info: u64) -> (bool, usize) {
        let curr_switch = curr_index_info & ONE_SWITCH > 0;

        let curr_index = {
            if curr_switch {
                let second_block = curr_index_info & SECOND_BLOCK;
                second_block as usize >> 21
            }
            else {
                let first_block = curr_index_info & FIRST_BLOCK;
                first_block as usize >> 42
            }
        };
        return (curr_switch, curr_index)
    }

    // Push key in correct array
    fn push(&self, key: usize, curr_switch: &bool, curr_index: &usize) {
        if !curr_switch {
            let curr_arc: ArcUsize = self.array_0[*curr_index].clone();
            curr_arc.swap(key, AcqRel);
            let curr_bit_arc: ArcBool = self.bit_array_0[*curr_index].clone();
            curr_bit_arc.swap(true, AcqRel);
        }
        else {
            let curr_arc: ArcUsize = self.array_1[*curr_index].clone();
            curr_arc.swap(key, AcqRel);
            let curr_bit_arc: ArcBool = self.bit_array_1[*curr_index].clone();
            curr_bit_arc.swap(true, AcqRel);
        }
    }

    // If thread can't update KmerEdge, add key to array
    pub fn add_to_array(&self, key: usize) {
        let curr_index_info = self.filled_indices_tracker.fetch_update(
            SeqCst, SeqCst, |x| {
                if x & ONE_SWITCH > 0 {
                    Some(x + INCR_SECOND + INCR_THIRD)
                }
                else {
                    Some(x + INCR_FIRST + INCR_THIRD)
                }
            }).unwrap();
        
        let (curr_switch, curr_index) = Self::parse_index_info(curr_index_info);
        self.push(key, &curr_switch, &curr_index);
    }
}


pub struct KmerEdgeSingle {
    kmer: usize,
    vertices_key: Vec<usize>,
    vertices_bit_array: Vec<bool>,
}
impl KmerEdgeSingle {

    // Make a new KmerEdge
    fn new(kmer: usize, vertices_count: usize) -> KmerEdgeSingle {
        KmerEdgeSingle {
            kmer,
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

pub struct KmerEdgeGroup {
    kmers: Vec<usize>,
    graph: WeakGraph,
    vertices_key: Vec<usize>,
    vertices_bit_array: Vec<bool>,
}
impl KmerEdgeGroup {

    // Merge KmerEdge objects into one
    fn new(kmer_edge_keys: &Vec<WeakUsize>, graph: WeakGraph) -> KmerEdgeGroup {
        // Create empty vectors
        let mut kmers: Vec<usize> = vec![];
        let mut vertices_key: Vec<usize> = vec![];
        let mut vertices_bit_array: Vec<bool> = vec![];

        // Pointer is valid, so we are safe
        let graph_ref: &Graph = unsafe {&*graph.upgrade().unwrap().load(Acquire)};

        // Append each KmerEdge's kmer and vertices info to aforementioned vectors
        for edge_key in kmer_edge_keys {
            let loaded_edge_key = edge_key.upgrade().unwrap().load(Acquire);
            let a_ptr_kmer_edge = graph_ref.edges[loaded_edge_key].clone();

            // We won't update kmer_edge, so we are safe
            let kmer_edge = unsafe{&*a_ptr_kmer_edge.load(Acquire)};

            kmers.append(&mut kmer_edge.get_kmers());
            if vertices_bit_array.is_empty() {
                vertices_key.append(&mut kmer_edge.get_vertices_key());
                vertices_bit_array.append(&mut kmer_edge.get_vertices_bit_array());
            }
            else {
                for key in kmer_edge.get_vertices_key() {
                    if !vertices_bit_array[key] {
                        vertices_bit_array[key] = true;
                        vertices_key.push(key);
                    }
                }
            }
        }

        // Make KmerEdgeGroup object
        KmerEdgeGroup {kmers, graph, vertices_key, vertices_bit_array}
    }

    // Add key to ProteinVertex
    fn add_vertex(&mut self, vertex_key: usize) {
        if !self.vertices_bit_array[vertex_key] {
            self.vertices_bit_array[vertex_key] = true;
            self.vertices_key.push(vertex_key);
        }
    }

    // Update proteins to reflect changes in graph
    fn update_protein_vertices(&self,
            new_key: WeakUsize, thread_count: u32) {

        // Pointer is valid, so we are safe
        let graph_ref: &Graph = unsafe {&*self.graph.upgrade().unwrap().load(Acquire)};

        // Make ARC pointers of shared variables 
        let vertices_key: Arc<Vec<usize>> = Arc::new(self.vertices_key.clone());
        let vertices_index: ArcUsize = Arc::new(AtomicUsize::new(0));
        let mut handles = vec![];

        for _ in 0..thread_count {
            let vertices: Arc<Vec<ArcProteinVertex>> = graph_ref.vertices.clone();
            let vertices_key: Arc<Vec<usize>> = vertices_key.clone();
            let vertices_index: ArcUsize = vertices_index.clone();
            let new_key: WeakUsize = new_key.clone();
            handles.push(spawn(move || {
                loop {
                    let curr_vertex_index = vertices_index.fetch_add(
                        1, AcqRel
                    );
                    if curr_vertex_index >= vertices_key.len() {
                        break;
                    }
                    let curr_key = vertices_key[curr_vertex_index];
                    let a_ptr_curr_vertex: &mut ArcProteinVertex = 
                        &mut vertices[curr_key].clone();

                    // No two threads will ever traverse the same vertex, so we are safe
                    let curr_vertex: &mut ProteinVertex = unsafe {
                        &mut *a_ptr_curr_vertex.load(Acquire)};
                    curr_vertex.add_edge(new_key.clone());

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
    pub fn get_kmers(&self) -> Vec<usize> {
        match self {
            Self::Single(val) => vec![val.kmer],
            Self::Group(val) => val.kmers.clone(),
        }
    }

    pub fn get_vertices_key(&self) -> Vec<usize> {
        match self {
            Self::Single(val) => val.vertices_key.clone(),
            Self::Group(val) => val.vertices_key.clone(),
        }
    }

    pub fn get_vertices_bit_array(&self) -> Vec<bool> {
        match self {
            Self::Single(val) => val.vertices_bit_array.clone(),
            Self::Group(val) => val.vertices_bit_array.clone(),
        }
    }

    pub fn add_vertex(&mut self, vertex_key: usize) {
        match self {
            Self::Single(val) => val.add_vertex(vertex_key),
            Self::Group(val) => val.add_vertex(vertex_key),
        }
    }

    pub fn update_protein_vertices(&self,
        new_key: WeakUsize, thread_count: u32) {
        match self {
            Self::Single(_) => panic!(
                "Can't call update_protein_vertices on single kmer"),
            Self::Group(val) => val.update_protein_vertices(
                new_key, thread_count),
        }
    }

    pub fn new(kmer: usize, vertices_count: usize) -> KmerEdge {
        Self::Single(KmerEdgeSingle::new(kmer, vertices_count))
    }
    
    pub fn merge(kmer_edge_keys: &Vec<WeakUsize>, graph: WeakGraph) -> KmerEdge {
        Self::Group(KmerEdgeGroup::new(kmer_edge_keys, graph))
    }
    
    // Find the relative frequency of vertices in other that are also in self
    pub fn compare(&self, other: &Self) -> f64 {
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