use crate::Graph;
// use crate::graph::vertex::ProteinVertex;

use std::fmt;
use std::sync::{Arc, Weak};
use std::sync::atomic::Ordering::{Acquire, AcqRel, Release};
use std::sync::atomic::{AtomicBool, AtomicPtr, AtomicU64, AtomicUsize};

// type ArcU64 = Arc<AtomicU64>;
// type ArcBool = Arc<AtomicBool>;
type ArcUsize = Arc<AtomicUsize>;
// type ArcProteinVertex = Arc<AtomicPtr<ProteinVertex>>;

type WeakUsize = Weak<AtomicUsize>;
type WeakGraph = Weak<AtomicPtr<Graph>>;

// const ONE_SWITCH: u64 = 
//     0b1000000000000000000000000000000000000000000000000000000000000000;
// const INCR_FIRST: u64 = 
//     0b0000000000000000000001000000000000000000000000000000000000000000;
// const FIRST_BLOCK: u64 = 
//     0b0111111111111111111111000000000000000000000000000000000000000000;
// const FIRST_EMPTY: u64 =
//     0b1000000000000000000000111111111111111111111111111111111111111111;
// const INCR_SECOND: u64 = 
//     0b0000000000000000000000000000000000000000001000000000000000000000;
// const SECOND_BLOCK: u64 = 
//     0b0000000000000000000000111111111111111111111000000000000000000000;
// const SECOND_EMPTY: u64 =
//     0b1111111111111111111111000000000000000000000111111111111111111111;
// const INCR_THIRD: u64 = 
//     0b0000000000000000000000000000000000000000000000000000000000000001;

// pub struct KmerEdgeHelper {
//     // thread that gets to switch this from true to false can update KmerEdge
//     update_helper: ArcBool,

//     // const-sized arrays containing vertex keys to add to KmerEdge
//     array_0: Vec<ArcUsize>,
//     array_1: Vec<ArcUsize>,

//     // bit-array informing whether other threads have populated vertex key array
//     bit_array_0: Vec<ArcBool>,
//     bit_array_1: Vec<ArcBool>,

//     //  1 + 21 + 21 + 21
//     //  zero_one switch ->  if zero, will look through first 21
//     //                      if one, will look through second 21
//     //  third 21 ->         number of times that tracker was updated (avoids ABA problem)
//     //  second 21 ->        if zero, threads will increment by one and retrieve previous
//     //  first 21 ->         if one, threads will increment by one and retrieve previous
//     filled_indices_tracker: ArcU64,
// }

// impl KmerEdgeHelper {
//     pub fn new() -> KmerEdgeHelper {
//         KmerEdgeHelper { 
//             update_helper: Arc::new(AtomicBool::new(true)), 
//             array_0: vec![0usize; 4000].into_iter()
//                 .map(|x| Arc::new(AtomicUsize::new(x))).collect(),
//             array_1: vec![0usize; 4000].into_iter()
//                 .map(|x| Arc::new(AtomicUsize::new(x))).collect(),
//             bit_array_0: vec![false; 4000].into_iter()
//                 .map(|x| Arc::new(AtomicBool::new(x))).collect(),
//             bit_array_1: vec![false; 4000].into_iter()
//                 .map(|x| Arc::new(AtomicBool::new(x))).collect(),
//             filled_indices_tracker: Arc::new(AtomicU64::new(0u64)),
//         }
//     }

//     // If returns true, then current thread is allowed to update KmerEdge
//     pub fn can_update_kmer_edge(&self) -> bool {
//         self.update_helper.swap(false, AcqRel)
//     }

//     // Add all keys in array as well as key in argument list to KmerEdge
//     pub fn update_kmer_edge(&self, key: usize, edge: &mut KmerEdge) {
//         let curr_index_info = self.filled_indices_tracker.fetch_update(
//             SeqCst, SeqCst, |x| {
//                 if x & ONE_SWITCH > 0 {
//                     Some((x - ONE_SWITCH + INCR_THIRD) & SECOND_EMPTY)
//                 }
//                 else {
//                     Some((x + ONE_SWITCH + INCR_THIRD) & FIRST_EMPTY)
//                 }
//         }).unwrap();

//         // Add key in argument list to array
//         let (curr_switch, curr_index) = Self::parse_index_info(curr_index_info);
//         self.push(key, &curr_switch, &curr_index);

//         // In case there are threads that are still running, come back to them later
//         let mut leftovers = vec![];

//         // Add keys in array 
//         for i in 0..(curr_index+1) {
//             if !curr_switch {
//                 if !self.bit_array_0[i].load(Acquire) {
//                     leftovers.push(i);
//                 }
//                 else {
//                     edge.add_vertex(self.array_0[i].load(Acquire));
//                     self.bit_array_0[i].swap(false, AcqRel);
//                 }
//             }
//             else {
//                 if !self.bit_array_1[i].load(Acquire) {
//                     leftovers.push(i);
//                 }
//                 else {
//                     edge.add_vertex(self.array_1[i].load(Acquire));
//                     self.bit_array_1[i].store(false, Release);
//                 }
//             }
//         }

//         // Add keys in array if not previously added
//         for i in leftovers {
//             if !curr_switch {
//                 while !self.bit_array_0[i].load(Acquire) {
//                     continue;
//                 }
//                 edge.add_vertex(self.array_0[i].load(Acquire));
//                 self.bit_array_0[i].store(false, Release);
//             }
//             else {
//                 while !self.bit_array_1[i].load(Acquire) {
//                     continue;
//                 }
//                 edge.add_vertex(self.array_1[i].load(Acquire));
//                 self.bit_array_1[i].store(false, Release);
//             }
//         }

//         // Now another thread can update KmerEdge if needed
//         self.update_helper.store(true, Release);
//     }

//     // During last run-through, make a final update to KmerEdge if some keys are leftover
//     pub fn update_kmer_edge_final(&self, edge: &mut KmerEdge) {
//         let (curr_switch, curr_index) = Self::parse_index_info(
//             self.filled_indices_tracker.load(Acquire));

//         for i in 0..curr_index {
//             if !curr_switch {
//                 edge.add_vertex(self.array_0[i].load(Acquire));
//                 self.bit_array_0[i].store(false, Release);
//             }
//             else {
//                 edge.add_vertex(self.array_1[i].load(Acquire));
//                 self.bit_array_1[i].store(false, Release);
//             }
//         }
//         self.update_helper.store(true, Release);
//     }

//     // Get switch orientation and current index of updating array
//     fn parse_index_info(curr_index_info: u64) -> (bool, usize) {
//         let curr_switch = curr_index_info & ONE_SWITCH > 0;

//         let curr_index = {
//             if curr_switch {
//                 let second_block = curr_index_info & SECOND_BLOCK;
//                 second_block as usize >> 21
//             }
//             else {
//                 let first_block = curr_index_info & FIRST_BLOCK;
//                 first_block as usize >> 42
//             }
//         };
//         return (curr_switch, curr_index)
//     }

//     // Push key in correct array
//     fn push(&self, key: usize, curr_switch: &bool, curr_index: &usize) {
//         if !curr_switch {
//             let curr_arc: ArcUsize = self.array_0[*curr_index].clone();
//             curr_arc.swap(key, AcqRel);
//             let curr_bit_arc: ArcBool = self.bit_array_0[*curr_index].clone();
//             curr_bit_arc.swap(true, AcqRel);
//         }
//         else {
//             let curr_arc: ArcUsize = self.array_1[*curr_index].clone();
//             curr_arc.swap(key, AcqRel);
//             let curr_bit_arc: ArcBool = self.bit_array_1[*curr_index].clone();
//             curr_bit_arc.swap(true, AcqRel);
//         }
//     }

//     // If thread can't update KmerEdge, add key to array
//     pub fn add_to_array(&self, key: usize) {
//         let curr_index_info = self.filled_indices_tracker.fetch_update(
//             SeqCst, SeqCst, |x| {
//                 if x & ONE_SWITCH > 0 {
//                     Some(x + INCR_SECOND + INCR_THIRD)
//                 }
//                 else {
//                     Some(x + INCR_FIRST + INCR_THIRD)
//                 }
//             }).unwrap();
        
//         let (curr_switch, curr_index) = Self::parse_index_info(curr_index_info);
//         self.push(key, &curr_switch, &curr_index);
//     }
// }


pub struct KmerEdgeSingle {
    kmer: usize,
    vertices_key: [ArcUsize; 2],
    first_key_is_set: AtomicBool,
}
impl KmerEdgeSingle {

    // Make a new KmerEdge
    fn new(kmer: usize) -> KmerEdgeSingle {
        KmerEdgeSingle {
            kmer,
            vertices_key: [Arc::new(AtomicUsize::new(0)), 
                Arc::new(AtomicUsize::new(0))],
            first_key_is_set: AtomicBool::new(false),
        }
    }

    // Add key to ProteinVertex
    fn add_vertex(&self, vertex_key: usize) {
        if ! self.first_key_is_set.swap(true, AcqRel) {
            self.vertices_key[0].store(vertex_key, Release);
        }
        else {
            self.vertices_key[1].store(vertex_key, Release);
        }
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

    // fn new_vertex_pair(pair: [usize; 2], graph: WeakGraph) -> KmerEdgeGroup {
    //     // Create empty vectors
    //     let mut kmers: Vec<usize> = vec![];
    //     let vertices_key: Vec<usize> = pair.as_slice().to_vec();
        

    //     // Pointer is valid, so we are safe
    //     let graph_ref: &Graph = unsafe {&*graph.upgrade().unwrap().load(Acquire)};

    //     let mut vertices_bit_array: Vec<bool> = vec![false; graph_ref.vertices.len()];
    //     vertices_bit_array[pair[0]] = true;
    //     vertices_bit_array[pair[1]] = true;

    //     let pair =

    //     // Append each KmerEdge's kmer and vertices info to aforementioned vectors
    //     for edge_key in kmer_edge_keys {
    //         let loaded_edge_key = edge_key.upgrade().unwrap().load(Acquire);
    //         let a_ptr_kmer_edge = graph_ref.edges[loaded_edge_key].clone();

    //         // We won't update kmer_edge, so we are safe
    //         let kmer_edge = unsafe{&*a_ptr_kmer_edge.load(Acquire)};

    //         kmers.append(&mut kmer_edge.get_kmers());
    //         if vertices_bit_array.is_empty() {
    //             vertices_key.append(&mut kmer_edge.get_vertices_key());
    //             vertices_bit_array.append(&mut kmer_edge.get_vertices_bit_array());
    //         }
    //         else {
    //             for key in kmer_edge.get_vertices_key() {
    //                 if !vertices_bit_array[key] {
    //                     vertices_bit_array[key] = true;
    //                     vertices_key.push(key);
    //                 }
    //             }
    //         }
    //     }

    //     // Make KmerEdgeGroup object
    //     KmerEdgeGroup {kmers, graph, vertices_key, vertices_bit_array}
    // }

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

    pub fn add_vertex(&self, vertex_key: usize) {
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

    // pub fn make_vertex_pair(pair: [usize; 2], graph: WeakGraph) -> KmerEdge {
    //     Self::Group(KmerEdgeGroup::new_vertex_pair(pair, graph))
    // }
    
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