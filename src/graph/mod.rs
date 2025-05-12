mod edge;
mod vertex;

use crate::protein::Protein;
use crate::graph::edge::KmerEdge;
use crate::graph::edge::KmerEdgeHelper;
use crate::graph::vertex::ProteinVertex;

use std::fmt;
use std::ptr;
use std::rc::Rc;
use std::boxed::Box;
use std::thread::spawn;
use std::time::Instant;
use std::sync::{Arc, mpsc, Weak};
use std::sync::atomic::{AtomicBool, AtomicPtr, AtomicUsize};
use std::sync::atomic::Ordering::{Acquire, AcqRel, Release};

type ArcVec<T> = Arc<Vec<T>>;
type ArcBool = Arc<AtomicBool>;
type ArcProtein = Arc<Protein>;
type ArcUsize = Arc<AtomicUsize>;
type ArcGraph = Arc<AtomicPtr<Graph>>;
type ArcKmerEdge = Arc<AtomicPtr<KmerEdge>>;
type ArcProteinVertex = Arc<AtomicPtr<ProteinVertex>>;
type ArcKmerEdgeHelper = Arc<AtomicPtr<KmerEdgeHelper>>;

type WeakUsize = Weak<AtomicUsize>;
type WeakGraph = Weak<AtomicPtr<Graph>>;

type OptionWeakUsize = Rc<Option<WeakUsize>>;

pub struct Graph {
    pub edges: ArcVec<ArcKmerEdge>,
    protein_list: ArcVec<ArcProtein>,
    vertices: ArcVec<ArcProteinVertex>,
    global_edge_keys: ArcVec<ArcUsize>,
    helpers : ArcVec<ArcKmerEdgeHelper>,
}
impl Graph {
    // Create a graph and return an ARC pointer to said graph
    pub fn new(kmer_count: usize, thread_count: usize, 
            protein_list: ArcVec<ArcProtein>) -> ArcGraph {

        // Create Graph object skeleton
        let now = Instant::now();
        let global_edge_keys: ArcVec<ArcUsize> = Arc::new((0..kmer_count)
            .into_iter().map(|k: usize| Arc::new(AtomicUsize::new(k)))
            .collect());
        let elapsed = now.elapsed().as_nanos();
        eprintln!("global_edge_keys vector construction time: {} nanosecs", 
            elapsed);
        
        let helpers: ArcVec<ArcKmerEdgeHelper> = Arc::new((0..kmer_count)
            .into_iter().map(|_| Arc::new(AtomicPtr::new(ptr::null_mut())))
            .collect());

        let vertices: ArcVec<ArcProteinVertex> = Arc::new((0..protein_list.len())
            .into_iter().map(|_| Arc::new(AtomicPtr::new(ptr::null_mut())))
            .collect());

        let edges: ArcVec<ArcKmerEdge> = Arc::new((0..kmer_count)
            .into_iter().map(|_| Arc::new(AtomicPtr::new(ptr::null_mut())))
            .collect());

        eprintln!("Empty vectors created");

        // Make ARC pointers of shared variables 
        let arc_kmer_count: ArcUsize = Arc::new(AtomicUsize::new(kmer_count));
        let arc_protein_count: ArcUsize = Arc::new(AtomicUsize::new(protein_list.len()));

        // Concurrently make new KmerEdge and KmerEdgeHelper structs
        let now = Instant::now();
        let edge_index: ArcUsize = Arc::new(AtomicUsize::new(0usize));
        let mut handles  = vec![];
        for _ in 0..thread_count {
            let edge_index: ArcUsize = edge_index.clone();
            let edges: ArcVec<ArcKmerEdge> = edges.clone();
            let arc_kmer_count: ArcUsize = arc_kmer_count.clone();
            let helpers: ArcVec<ArcKmerEdgeHelper> = helpers.clone();
            let arc_protein_count: ArcUsize = arc_protein_count.clone();
            handles.push(spawn(move || {
                let vertices_count = arc_protein_count.load(Acquire);
                let kmer_count = arc_kmer_count.load(Acquire);
                loop {

                    // If we already have enough edges, break
                    let curr_index = edge_index.fetch_add(1, AcqRel);
                    if curr_index >= kmer_count {
                        break;
                    }

                    // Get atomic ptrs to edge and helper
                    let a_ptr_edge: ArcKmerEdge = edges[curr_index].clone();
                    let a_ptr_helper: ArcKmerEdgeHelper = helpers[curr_index].clone();

                    // Create new KmerEdge and KmerEdgeHelper structs
                    let curr_edge: KmerEdge = KmerEdge::new(
                        curr_index, vertices_count);
                    let curr_helper: KmerEdgeHelper = KmerEdgeHelper::new();

                    // Get raw ptrs to new KmerEdge and KmerEdgeHelper structs
                    let ptr_edge = Box::into_raw(Box::new(curr_edge));
                    let ptr_helper = Box::into_raw(Box::new(curr_helper));

                    // Store raw ptrs into atomic ptrs
                    a_ptr_edge.store(ptr_edge, Release);
                    a_ptr_helper.store(ptr_helper, Release);
                }
                
            }));
        }

        // Wait for threads to end
        for handle in handles {
            handle.join().unwrap();
        }
        
        let elapsed = now.elapsed().as_nanos();
        eprintln!("edges vector construction time: {} nanosecs", 
            elapsed);

        // Construct graph
        let graph: Graph = Graph {
            vertices,
            edges,
            protein_list,
            global_edge_keys,
            helpers,
        };

        let graph_ptr: *mut Graph = Box::into_raw(Box::new(graph));
        let graph_arc: ArcGraph = Arc::new(AtomicPtr::new(graph_ptr));

        let graph: &Graph = unsafe{& *graph_arc.load(Acquire)};

        // Concurrently make new ProteinVertex structs and update KmerEdge structs
        let now = Instant::now();
        let vertex_index: ArcUsize = Arc::new(AtomicUsize::new(0usize));
        let mut handles  = vec![];
        for _ in 0..thread_count {
            let vertices: ArcVec<ArcProteinVertex> = graph.vertices.clone();
            let arc_protein_count: ArcUsize = arc_protein_count.clone();
            let vertex_index: ArcUsize = vertex_index.clone();
            let graph_arc: ArcGraph = graph_arc.clone();
            handles.push(spawn(move || {
                let protein_count: usize = arc_protein_count.load(Acquire);
                loop {

                    // If we already have enough vertices, break
                    let curr_index = vertex_index.fetch_add(1, AcqRel);
                    if curr_index >= protein_count {
                        break;
                    }

                    // Get atomic ptrs to vertex
                    let a_ptr_vertex: ArcProteinVertex = vertices[curr_index].clone();

                    // Create new ProteinVertex
                    let curr_vertex: ProteinVertex = ProteinVertex::new(
                        curr_index, Arc::downgrade(&graph_arc));
                    curr_vertex.update_graph_edges();

                    // Get raw ptrs to new KmerEdge and KmerEdgeHelper structs
                    let ptr_vertex = Box::into_raw(Box::new(curr_vertex));

                    // Store raw ptrs into atomic ptrs
                    a_ptr_vertex.store(ptr_vertex, Release);
                }
            }));
        }

        // Wait for threads to end
        for handle in handles {
            handle.join().unwrap();
        }

        // Concurrently update KmerEdge structs one last time
        let edge_index: ArcUsize = Arc::new(AtomicUsize::new(0usize));
        let mut handles  = vec![];
        for _ in 0..thread_count {
            let edge_helpers: ArcVec<ArcKmerEdgeHelper> = 
                graph.helpers.clone();
            let edges: ArcVec<ArcKmerEdge> = graph.edges.clone();
            let arc_kmer_count: ArcUsize = arc_kmer_count.clone();
            let edge_index: ArcUsize = edge_index.clone();
            
            handles.push(spawn(move || {
                let kmer_count = arc_kmer_count.load(Acquire);
                loop {

                    // If we have gone through all edges, break
                    let curr_index = edge_index.fetch_add(1, AcqRel);
                    if curr_index >= kmer_count {
                        break;
                    }

                    // Retrieve current helper and edge
                    let a_ptr_curr_helper: ArcKmerEdgeHelper = edge_helpers[curr_index].clone();
                    let a_ptr_curr_edge: ArcKmerEdge = edges[curr_index].clone();

                    // Only one thread can access a specific edge at a time, so we are safe
                    let curr_edge: &mut KmerEdge = unsafe {
                        &mut *a_ptr_curr_edge.load(Acquire)};
                    let curr_helper: &KmerEdgeHelper = unsafe {
                        & *a_ptr_curr_helper.load(Acquire)};

                    curr_helper.update_kmer_edge_final(curr_edge);
                }
            }));
        }
        
        // Wait for threads to end
        for handle in handles {
            handle.join().unwrap();
        }
        let elapsed = now.elapsed().as_nanos();
        eprintln!("vertices vector construction time: {} nanosecs", 
            elapsed);

        graph_arc
    }

    // Combine edges with lots of overlap
    pub fn combine_edges(&self, graph_weak: WeakGraph, thread_count: u32) {
        // Get current keys to edges
        let mut unordered_edge_keys: Vec<OptionWeakUsize> = (0..self.edges.len())
            .into_iter().map(|i| Rc::new(Some(Arc::downgrade(
            &self.global_edge_keys[i])))).collect();

        // Make a second vector of current keys, ordered by edge size
        let mut ordered_edge_keys: Vec<OptionWeakUsize> = unordered_edge_keys.clone();
        ordered_edge_keys.sort_by(|a: &OptionWeakUsize, b: &OptionWeakUsize| {
            let cloned_a: Option<WeakUsize> = a.as_ref().clone();
            let cloned_b: Option<WeakUsize> = b.as_ref().clone();
            let key_a: usize = cloned_a.unwrap().upgrade().unwrap().load(Acquire);
            let key_b: usize = cloned_b.unwrap().upgrade().unwrap().load(Acquire);
            let edge_a_ptr: ArcKmerEdge = self.edges[key_a].clone();
            let edge_b_ptr: ArcKmerEdge = self.edges[key_b].clone();
            let edge_a: &KmerEdge = unsafe {& *edge_a_ptr.load(Acquire)};
            let edge_b: &KmerEdge = unsafe {& *edge_b_ptr.load(Acquire)};
            let edge_a_size: usize = edge_a.get_vertices_key().len();
            let edge_b_size: usize = edge_b.get_vertices_key().len();
            edge_b_size.cmp(&edge_a_size)});
        
        // Traverse through each edge
        for edge_key in ordered_edge_keys {
            // If edge was merged, skip
            if edge_key.is_none() {
                continue;
            }
            let now = Instant::now();

            // Get actual pointer to KmerEdge
            let unwraped_edge_key: WeakUsize = edge_key.as_ref().clone().unwrap();
            let loaded_edge_key: usize = unwraped_edge_key.upgrade().unwrap().load(Acquire);
            let edge_ptr: ArcKmerEdge = self.edges[loaded_edge_key].clone();

            // No one is updating edge, so we are safe
            let edge = unsafe {& *edge_ptr.load(Acquire)};

            // Make ARC pointers of shared variables   
            let vertices_index: ArcUsize = 
                Arc::new(AtomicUsize::new(1));
            let vertices_keys: ArcVec<usize> = 
                Arc::new(edge.get_vertices_key());
            let first_vertex_ptr: ArcProteinVertex = 
                self.vertices[vertices_keys[0]].clone();

            // No one is updating vertex, so we are safe
            let first_vertex: &ProteinVertex = unsafe {
                & *first_vertex_ptr.load(Acquire)};

            let other_edges_bitarray: ArcVec<AtomicBool> = 
                Arc::new(first_vertex.get_edges_bit_array(self.edges.len())
                    .iter().map(|b: &bool| AtomicBool::new(*b)).collect());
            let mut other_edges_keys: Vec<WeakUsize> = 
                first_vertex.get_edges_key().clone();

            // Make channels and thread handles
            let (send, recv) = mpsc::channel::<ArcVec<WeakUsize>>();
            let mut handles  = vec![];

            // Retrieve other edges that share a vertex with current edge
            for _ in 0..thread_count {
                let vertices_index: ArcUsize = vertices_index.clone();
                let vertices_keys: ArcVec<usize> = vertices_keys.clone();
                let send: mpsc::Sender<ArcVec<WeakUsize>> = send.clone();
                let vertices: ArcVec<ArcProteinVertex> = self.vertices.clone();
                let mut other_edges_keys_split: ArcVec<WeakUsize> = Arc::new(vec![]);
                let other_edges_bitarray: ArcVec<AtomicBool> = other_edges_bitarray.clone();
                handles.push(spawn(move || {
                    loop {
                        let curr_vertices_index = vertices_index.fetch_add(
                            1, AcqRel
                        );
                        if curr_vertices_index >= vertices_keys.len() {
                            let _ = send.send(other_edges_keys_split);
                            break;
                        }

                        let curr_vertex_key = vertices_keys[curr_vertices_index];
                        let curr_vertex_ptr: &ArcProteinVertex = 
                            &vertices[curr_vertex_key];

                        // No one is updating vertex, so we are safe
                        let curr_vertex: &ProteinVertex = unsafe {
                            & *curr_vertex_ptr.load(Acquire)};

                        let more_edges_keys: &Vec<WeakUsize> = 
                            &curr_vertex.get_edges_key();
                        for key in more_edges_keys {
                            let loaded_key = key.upgrade().unwrap().load(Acquire);
                            if !other_edges_bitarray[loaded_key].swap(true, AcqRel) {
                                Arc::get_mut(&mut other_edges_keys_split)
                                    .unwrap().push(key.clone());
                            }
                        }
                    }
                }))
            }
            drop(send);
            
            for handle in handles {
                handle.join().unwrap();
            }

            for _ in 0..thread_count {
                let mut other_edges_keys_split: ArcVec<WeakUsize> = recv.recv().unwrap();
                other_edges_keys.append(
                    Arc::get_mut(&mut other_edges_keys_split).unwrap());
            }
            drop(recv);
            drop(vertices_index);
            drop(other_edges_bitarray);

            // Gather keys to other edges
            if other_edges_keys.len() == 1 {
                let elapsed = now.elapsed().as_nanos();
                eprintln!("edge # {} traversal time: {} nanosecs", 
                    loaded_edge_key,
                    elapsed);
                continue;
            }

            // Sort keys
            other_edges_keys.sort_by(|a: &WeakUsize, b: &WeakUsize| 
                a.upgrade().unwrap().load(Acquire)
                .cmp(&b.upgrade().unwrap().load(Acquire)));

            // Make ARC pointers of shared variables 
            let other_edges_index: ArcUsize = Arc::new(AtomicUsize::new(0));
            let other_edges_keys: ArcVec<WeakUsize> = Arc::new(other_edges_keys);
            let mut edges_to_merge_with: Vec<WeakUsize> = Vec::new();

            // Make channels and thread handles
            let (send, recv) = mpsc::channel::<ArcVec<WeakUsize>>();
            let mut handles= vec![];

            // Find which of the other edges can be merged with current edge
            for _ in 0..thread_count {
                let edge_ptr: ArcKmerEdge = edge_ptr.clone();
                let edges: ArcVec<ArcKmerEdge> = self.edges.clone();
                let send: mpsc::Sender<ArcVec<WeakUsize>> = send.clone();
                let other_edges_index: ArcUsize = other_edges_index.clone();
                let other_edges_keys: ArcVec<WeakUsize> = other_edges_keys.clone();
                let mut edges_to_merge_with_split: ArcVec<WeakUsize> = Arc::new(vec![]);
                handles.push(spawn(move || {
                    // No one is updating edge, so we are safe
                    let edge = unsafe {& *edge_ptr.load(Acquire)};
                    loop {
                        let other_edges_index = other_edges_index.fetch_add(
                            1, AcqRel);
                        if other_edges_index >= other_edges_keys.len() {
                            let _ = send.send(edges_to_merge_with_split);
                            break;
                        }

                        let other_key: &WeakUsize = &other_edges_keys[other_edges_index];
                        let loaded_other_key: usize = other_key.upgrade()
                            .unwrap().load(Acquire);
                        let other_edge_ptr: ArcKmerEdge = edges[loaded_other_key].clone();
                        
                        // No one is updating other edge, so we are safe
                        let other_edge = unsafe {& *other_edge_ptr.load(Acquire)};

                        let other_cmp_edge: f64 = edge.compare(other_edge);
                        // let edge_cmp_other: f64 = other_edge.compare(edge);
                        if other_cmp_edge == 1.0 {
                            Arc::get_mut(&mut edges_to_merge_with_split)
                                .unwrap().push(other_key.clone());
                        }
                        // else if edge_cmp_other >= 0.8 && other_cmp_edge >= 0.8 {
                        //     Arc::get_mut(&mut edges_to_merge_with_split)
                        //         .unwrap().push(other_key.clone());
                        // }
                    }
                }))
            }
            drop(send);
            
            for handle in handles {
                handle.join().unwrap();
            }

            for _ in 0..thread_count {
                let mut edges_to_merge_with_split: ArcVec<WeakUsize> = recv.recv().unwrap();
                edges_to_merge_with.append(
                    Arc::get_mut(&mut edges_to_merge_with_split).unwrap());
            }

            drop(recv);
            drop(other_edges_keys);
            drop(other_edges_index);
            if edges_to_merge_with.len() == 1 {
                let elapsed = now.elapsed().as_nanos();
                eprintln!("edge # {} traversal time: {} nanosecs", 
                    loaded_edge_key,
                    elapsed);
                continue;
            }

            // Merge edges into one
            let kmer_edge_group: KmerEdge = KmerEdge::merge(
                &edges_to_merge_with, graph_weak.clone());
            kmer_edge_group.update_protein_vertices(unwraped_edge_key, thread_count);

            // Get raw ptrs to new KmerEdge and KmerEdgeHelper structs
            let edge_raw: *mut KmerEdge = Box::into_raw(Box::new(kmer_edge_group));

            // Swap and deallocate+delete previous
            let edge_raw = edge_ptr.swap(edge_raw, AcqRel);
            drop(unsafe { Box::from_raw(edge_raw) });
            
            for old_key in &*edges_to_merge_with {
                let loaded = old_key.upgrade().unwrap().load(Acquire);
                if loaded != loaded_edge_key {
                    unsafe {
                        *Rc::get_mut_unchecked(&mut unordered_edge_keys[loaded]) = None;
                    }
                }
            }
            let elapsed = now.elapsed().as_nanos();
            eprintln!("edge # {} traversal time: {} nanosecs", 
                loaded_edge_key,
                elapsed);
        }

        // Remove edges that have been merged with a bigger one
        let mut_self = unsafe {&mut *graph_weak.upgrade().unwrap().load(Acquire)};
        let mut substractive = 0usize;
        for (index, edge_key) in unordered_edge_keys.iter().enumerate() {
            mut_self.remove_edges_marked_for_deletions(
                edge_key.is_none(), &index, &mut substractive);
        }

        // Remove weak pointers to non-existing edges
        let mut handles= vec![];
        let vertices_index: ArcUsize = Arc::new(AtomicUsize::new(0));

        for _ in 0..thread_count {
            let vertices_index: ArcUsize = vertices_index.clone();
            let vertices: ArcVec<ArcProteinVertex> = mut_self.vertices.clone();

            handles.push(spawn(move || {
                loop {
                    let curr_vertices_index = vertices_index.fetch_add(1, AcqRel);
                    if curr_vertices_index >= 10619 {
                        break;
                    }

                    // No two threads will ever traverse the same vertex, so we are safe
                    let curr_vertex: &mut ProteinVertex = unsafe {
                        &mut *vertices[curr_vertices_index].load(Acquire)};
                    curr_vertex.remove_edges();
                }
            }))
        }

        for handle in handles {
            handle.join().unwrap();
        }
    }

    fn remove_edges_marked_for_deletions(&mut self, 
            remove_key: bool, index: &usize, substractive: &mut usize){
        if remove_key {
            // Currently in a single threaded environment, so we are safe
            unsafe {
                // must remove and dealloc
                let edge = Arc::get_mut_unchecked(&mut self.edges)
                    .remove(index - *substractive);
                let edge_ptr = edge.load(Acquire);
                drop(Box::from_raw(edge_ptr));
                drop(edge);

                let helper = Arc::get_mut_unchecked(&mut self.helpers)
                    .remove(index - *substractive);
                let helper_ptr = helper.load(Acquire);
                drop(Box::from_raw(helper_ptr));
                drop(helper);

                let key = Arc::get_mut_unchecked(&mut self.global_edge_keys)
                    .remove(index - *substractive);
                drop(key);
            }
            *substractive += 1;
        }
        else if *substractive > 0 {
            self.global_edge_keys[index - *substractive]
                .fetch_sub(*substractive, AcqRel);
        }
    }

    // Remove edges without diverging AMR labels
    pub fn remove_uninteresting_edges(&self, graph_weak: WeakGraph, thread_count: u32) {
        // Get current keys to edges
        let bitarray_keys_vec: ArcVec<ArcBool> = Arc::new((0..self.edges.len())
            .into_iter().map(|_| Arc::new(AtomicBool::new(true))).collect());

        // Make ARC pointers of shared variables
        let bitarray_vec_index: ArcUsize = Arc::new(AtomicUsize::new(0));
        let mut handles= vec![];

        // Find edges that only represent one amr class
        for _ in 0..thread_count {
            let bitarray_keys_vec: ArcVec<ArcBool> = bitarray_keys_vec.clone();
            let proteins: ArcVec<ArcProtein> = self.protein_list.clone();
            let bitarray_vec_index: ArcUsize = bitarray_vec_index.clone();
            let edges: ArcVec<ArcKmerEdge> = self.edges.clone();
            handles.push(spawn(move || {
                let edge_count = edges.len();
                loop {
                    let curr_index = bitarray_vec_index.fetch_add(1, AcqRel);
                    if curr_index >= edge_count {
                        break;
                    }
                    
                    // No two threads will ever traverse the same edge, so we are safe
                    let edge_ptr: ArcKmerEdge = edges[curr_index].clone();
                    let edge_ref: &KmerEdge = unsafe {& *edge_ptr.load(Acquire)};
                    let vertices_key: Vec<usize> = edge_ref.get_vertices_key();

                    let mut at_least_two_classes: bool = false;
                    let first_amr_class: &str = proteins[vertices_key[0]].get_amr_class();

                    // Traverse through vertices in edge to see 
                    // if edge represents at least two amr classes
                    for vertex_key in &vertices_key[1..] {
                        if *first_amr_class == *proteins[*vertex_key].get_amr_class() {
                            at_least_two_classes = true;
                            break;
                        }
                    }

                    // Set bit to 0 so we know to remove edge later on
                    if !at_least_two_classes {
                        bitarray_keys_vec[curr_index].store(false, Release);
                    }

                }
            }));
        }

        for handle in handles {
            handle.join().unwrap();
        }

        // Remove edges that represent only one amr class
        let mut_self = unsafe {&mut *graph_weak.upgrade().unwrap().load(Acquire)};
        let mut substractive = 0usize;
        for (index, edge_key) in bitarray_keys_vec.iter().enumerate() {
            mut_self.remove_edges_marked_for_deletions(
                edge_key.load(Acquire), &index, &mut substractive);
        }

        // Remove weak pointers to non-existing edges
        let mut handles= vec![];
        let vertices_index: ArcUsize = Arc::new(AtomicUsize::new(0));

        for _ in 0..thread_count {
            let vertices_index: ArcUsize = vertices_index.clone();
            let vertices: ArcVec<ArcProteinVertex> = mut_self.vertices.clone();

            handles.push(spawn(move || {
                loop {
                    let curr_vertices_index = vertices_index.fetch_add(1, AcqRel);
                    if curr_vertices_index >= 10619 {
                        break;
                    }

                    // No two threads will ever traverse the same vertex, so we are safe
                    let curr_vertex: &mut ProteinVertex = unsafe {
                        &mut *vertices[curr_vertices_index].load(Acquire)};
                    curr_vertex.remove_edges();
                }
            }))
        }

        for handle in handles {
            handle.join().unwrap();
        }
        
    }
}

impl fmt::Debug for Graph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Graph")
            .field("Kmers", &self.edges)
            .field("Proteins", &self.vertices).finish()
    }
}
