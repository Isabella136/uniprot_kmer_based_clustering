mod edge;
mod vertex;

use crate::protein::Protein;
use crate::graph::edge::KmerEdge;
// use crate::graph::edge::KmerEdgeHelper;
use crate::graph::vertex::ProteinVertex;

use std::fmt;
use std::ptr;
use std::rc::Rc;
use std::boxed::Box;
use std::thread::spawn;
use std::time::Instant;
use std::sync::{Arc, mpsc, Weak};
use std::sync::atomic::Ordering::{Acquire, AcqRel, Release};
use std::sync::atomic::{AtomicBool, AtomicPtr, AtomicUsize, AtomicU64};

type ArcU64 = Arc<AtomicU64>;
type ArcVec<T> = Arc<Vec<T>>;
type ArcBool = Arc<AtomicBool>;
type ArcProtein = Arc<Protein>;
type ArcUsize = Arc<AtomicUsize>;
type ArcGraph = Arc<AtomicPtr<Graph>>;
type ArcKmerEdge = Arc<AtomicPtr<KmerEdge>>;
type ArcProteinVertex = Arc<AtomicPtr<ProteinVertex>>;
// type ArcKmerEdgeHelper = Arc<AtomicPtr<KmerEdgeHelper>>;

type WeakUsize = Weak<AtomicUsize>;
type WeakGraph = Weak<AtomicPtr<Graph>>;

type OptionWeakUsize = Rc<Option<WeakUsize>>;

pub struct Graph {
    pub edges: ArcVec<ArcKmerEdge>,
    protein_list: ArcVec<ArcProtein>,
    vertices: ArcVec<ArcProteinVertex>,
    global_edge_keys: ArcVec<ArcUsize>,
}
impl Graph {
    // Create a graph and return an ARC pointer to said graph
    pub fn new(kmer_freq: Vec<usize>, thread_count: usize, 
            protein_list: ArcVec<ArcProtein>) -> ArcGraph {

        // Calculate edge count
        let kmer_count = kmer_freq.len();
        let edges_per_kmer: Vec<usize> = (0..kmer_count)
            .map(|i| kmer_freq[i] * (kmer_freq[i] - 1) / 2).collect();
        let edge_count_prefix_sum: ArcVec<usize> = Arc::new((0..kmer_count)
            .map(|i| edges_per_kmer[0..i+1].iter().sum()).collect());
        let edge_count = edge_count_prefix_sum.last().unwrap();

        // Create Graph object skeleton
        let now = Instant::now();
        let global_edge_keys: ArcVec<ArcUsize> = Arc::new((0..*edge_count)
            .map(|k: usize| Arc::new(AtomicUsize::new(k))).collect());
        let elapsed = now.elapsed().as_nanos();
        eprintln!("global_edge_keys vector construction time: {} nanosecs", 
            elapsed);
        
        let vertices: ArcVec<ArcProteinVertex> = Arc::new((0..protein_list.len())
            .map(|_| Arc::new(AtomicPtr::new(ptr::null_mut()))).collect());

        let edges: ArcVec<ArcKmerEdge> = Arc::new((0..*edge_count)
            .map(|_| Arc::new(AtomicPtr::new(ptr::null_mut()))).collect());

        eprintln!("Empty vectors created");

        // Make ARC pointers of shared variables 
        let arc_edge_count: ArcUsize = Arc::new(AtomicUsize::new(*edge_count));
        let arc_protein_count: ArcUsize = Arc::new(AtomicUsize::new(protein_list.len()));

        // Concurrently make new KmerEdge structs
        let now: Instant = Instant::now();

        // Left 32 bits =   current kmer
        // Right 32 bits =  current edge
        let kmer_edge_indices: ArcU64 = Arc::new(AtomicU64::new(0));
        let mut handles  = vec![];
        for _ in 0..thread_count {
            let edges: ArcVec<ArcKmerEdge> = edges.clone();
            let arc_edge_count: ArcUsize = arc_edge_count.clone();
            let kmer_edge_indices: ArcU64 = kmer_edge_indices.clone();
            let edge_count_prefix_sum: ArcVec<usize> = edge_count_prefix_sum.clone();
            handles.push(spawn(move || {
                let edge_count = arc_edge_count.load(Acquire);
                loop {
                    // If we already have enough edges, break
                    let curr_kmer_edge_index = kmer_edge_indices.fetch_add(1, AcqRel);
                    let curr_edge_index = (curr_kmer_edge_index % (1 << 32)) as usize;
                    if curr_edge_index >= edge_count {
                        break;
                    }

                    // Find the actual current kmer
                    let mut curr_kmer = (curr_kmer_edge_index >> 32) as usize;
                    if edge_count_prefix_sum[curr_kmer] == curr_edge_index {
                        kmer_edge_indices.fetch_add(1 << 32, AcqRel);
                        curr_kmer += 1;
                    }
                    while edge_count_prefix_sum[curr_kmer] <= curr_edge_index {
                        curr_kmer += 1;
                    }

                    // Get atomic ptrs to edge
                    let a_ptr_edge: ArcKmerEdge = edges[curr_edge_index].clone();

                    // Create new KmerEdge struct
                    let curr_edge: KmerEdge = KmerEdge::new(curr_kmer);

                    // Get raw ptr to new KmerEdge struct
                    let ptr_edge = Box::into_raw(Box::new(curr_edge));

                    // Store raw ptr into atomic ptr
                    a_ptr_edge.store(ptr_edge, Release);
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
        };

        let graph_ptr: *mut Graph = Box::into_raw(Box::new(graph));
        let graph_arc: ArcGraph = Arc::new(AtomicPtr::new(graph_ptr));

        let graph: &Graph = unsafe{& *graph_arc.load(Acquire)};

        // Concurrently make new ProteinVertex structs and update KmerEdge structs
        let now = Instant::now();

        // Goal: keep track of which edges have been assigned the first vertex
        // and which edges have also been assigned the second vertex


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
        let elapsed = now.elapsed().as_nanos();
        eprintln!("vertices vector construction time: {} nanosecs", 
            elapsed);

        graph_arc
    }

    // Combine edges with lots of overlap
    pub fn combine_edges(&self, graph_weak: WeakGraph, thread_count: u32) {
        // Get current keys to edges
        let mut option_edge_keys: Vec<OptionWeakUsize> = (0..self.edges.len())
            .into_iter().map(|i| Rc::new(Some(Arc::downgrade(
            &self.global_edge_keys[i])))).collect();
        
        // Traverse through each edge
        for loaded_edge_key in 0..self.edges.len() {
            // If edge was merged, skip
            if option_edge_keys[loaded_edge_key].is_none() {
                continue;
            }
            let now = Instant::now();

            // Get actual pointer to KmerEdge
            let edge_ptr: ArcKmerEdge = self.edges[loaded_edge_key].clone();

            // No one is updating edge, so we are safe
            let edge = unsafe {& *edge_ptr.load(Acquire)};

            // Get vertices associated with this edge
            let vertices_keys: [usize; 2] = edge.get_vertices_key();

            // Get edges containing first vertex
            // No one is updating vertex, so we are safe
            let first_vertex_ptr: ArcProteinVertex = 
                self.vertices[vertices_keys[0]].clone();
            let first_vertex: &ProteinVertex = unsafe {
                & *first_vertex_ptr.load(Acquire)};
            if first_vertex.get_edges_key().clone().len() == 1 {
                let elapsed = now.elapsed().as_nanos();
                eprintln!("edge # {} traversal time: {} nanosecs", 
                    loaded_edge_key,
                    elapsed);
                continue;
            }
            let first_vertex_edges_bitarray: ArcVec<bool> = 
                Arc::new(first_vertex.get_edges_bit_array(self.edges.len()));

            // Get edges containing second vertex
            // No one is updating vertex, so we are safe
            let second_vertex_ptr: ArcProteinVertex = 
                self.vertices[vertices_keys[1]].clone();
            let second_vertex: &ProteinVertex = unsafe {
                & *second_vertex_ptr.load(Acquire)};
            let second_vertex_edges_keys: ArcVec<WeakUsize> = 
                Arc::new(second_vertex.get_edges_key().clone());
            if second_vertex_edges_keys.len() == 1 {
                let elapsed = now.elapsed().as_nanos();
                eprintln!("edge # {} traversal time: {} nanosecs", 
                    loaded_edge_key,
                    elapsed);
                continue;
            }

            // Make ARC pointers of shared variables
            let second_vertex_edge_index: ArcUsize = Arc::new(AtomicUsize::new(0));

            // Make channels and thread handles
            let (send, recv) = mpsc::channel::<ArcVec<WeakUsize>>();
            let mut handles  = vec![];

            // Retrieve other edges that share both vertices with current edge
            for _ in 0..thread_count {
                let send: mpsc::Sender<ArcVec<WeakUsize>> = send.clone();
                let keys_index: ArcUsize = second_vertex_edge_index.clone();
                let keys: ArcVec<WeakUsize> = second_vertex_edges_keys.clone();
                let bitarray: ArcVec<bool> = first_vertex_edges_bitarray.clone();
                let mut dedup_edge_keys_split: ArcVec<WeakUsize> = Arc::new(vec![]);
                
                handles.push(spawn(move || {
                    loop {
                        let curr_index = keys_index.fetch_add(
                            1, AcqRel
                        );
                        if curr_index >= keys.len() {
                            let _ = send.send(dedup_edge_keys_split);
                            break;
                        }

                        let curr_key = keys[curr_index].clone();
                        let loaded_key = curr_key.upgrade().unwrap().load(Acquire);
                        if bitarray[loaded_key] {
                            Arc::get_mut(&mut dedup_edge_keys_split)
                                .unwrap().push(curr_key.clone());
                        }
                    }
                }))
            }
            drop(send);
            
            for handle in handles {
                handle.join().unwrap();
            }

            let mut edges_to_merge_with: Vec<WeakUsize> = vec![];

            for _ in 0..thread_count {
                let mut other_edges_keys_split: ArcVec<WeakUsize> = recv.recv().unwrap();
                edges_to_merge_with.append(
                    Arc::get_mut(&mut other_edges_keys_split).unwrap());
            }
            drop(recv);
            drop(second_vertex_edge_index);
            drop(first_vertex_edges_bitarray);

            // Sort keys
            edges_to_merge_with.sort_by(|a: &WeakUsize, b: &WeakUsize| 
                a.upgrade().unwrap().load(Acquire)
                .cmp(&b.upgrade().unwrap().load(Acquire)));

            
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

            // Get raw ptrs to new KmerEdge and KmerEdgeHelper structs
            let edge_raw: *mut KmerEdge = Box::into_raw(Box::new(kmer_edge_group));

            // Swap and deallocate+delete previous
            let edge_raw = edge_ptr.swap(edge_raw, AcqRel);
            drop(unsafe { Box::from_raw(edge_raw) });
            
            for old_key in &*edges_to_merge_with {
                let loaded = old_key.upgrade().unwrap().load(Acquire);
                if loaded != loaded_edge_key {
                    unsafe {
                        *Rc::get_mut_unchecked(&mut option_edge_keys[loaded]) = None;
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
        for (index, edge_key) in option_edge_keys.iter().enumerate() {
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
        // Get bitarray representing edges to remove
        let bitarray_keys_vec: ArcVec<ArcBool> = Arc::new((0..self.edges.len())
            .into_iter().map(|_| Arc::new(AtomicBool::new(false))).collect());

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
                    let vertices_key: [usize; 2] = edge_ref.get_vertices_key();

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

                    // Set bit to 1 so we know to remove edge later on
                    if !at_least_two_classes {
                        bitarray_keys_vec[curr_index].store(true, Release);
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

    pub fn get_divergent_protein_pairs(&self, thread_count: u32) {
        // Get bitarray representing vertex pairs we've already visited
        let vertex_count = self.vertices.len();
        let pair_count = vertex_count*(vertex_count-1)/2;

        let bitarray_edge_visited: ArcVec<ArcBool> = Arc::new((0..self.edges.len())
            .into_iter().map(|_| Arc::new(AtomicBool::new(false))).collect());
        // let divergent_pair: ArcVec<ArcKmerEdge> = Arc::new((0..pair_count)
        //     .into_iter().map(|_| Arc::new(AtomicBool::new(false))).collect());

        // Make ARC pointers of shared variables
        let edge_index: ArcUsize = Arc::new(AtomicUsize::new(0));
        // let mut handles= vec![];

        // Find protein pairs
        for _ in 0..thread_count {
            // let divergent_pair: ArcVec<ArcKmerEdge> = divergent_pair.clone();
            let bitarray_edge_visited: ArcVec<ArcBool> = bitarray_edge_visited.clone();
            let proteins: ArcVec<ArcProtein> = self.protein_list.clone();
            let edges: ArcVec<ArcKmerEdge> = self.edges.clone();
            let edge_index: ArcUsize = edge_index.clone();
            // handles.push(spawn(move || {
            //     let edge_count = edges.len();
            //     loop {
            //         let curr_index = edge_index.fetch_add(1, AcqRel);
            //         if curr_index >= edge_count {
            //             break;
            //         }
                    
            //         // No two threads will ever traverse the same edge, so we are safe
            //         let edge_ptr: ArcKmerEdge = edges[curr_index].clone();
            //         let edge_ref: &KmerEdge = unsafe {& *edge_ptr.load(Acquire)};
            //         let vertices_key: Vec<usize> = edge_ref.get_vertices_key();

            //         let class_separated: [Vec<usize>; 30];
            //         let first_amr_class: &str = proteins[vertices_key[0]].get_amr_class();

            //         // Traverse through vertices in edge to see 
            //         // if edge represents at least two amr classes
            //         for vertex_key in &vertices_key[1..] {
            //             if *first_amr_class == *proteins[*vertex_key].get_amr_class() {
            //                 at_least_two_classes = true;
            //                 break;
            //             }
            //         }

            //         // Set bit to 1 so we know to remove edge later on
            //         if !at_least_two_classes {
            //             bitarray_keys_vec[curr_index].store(true, Release);
            //         }

            //     }
            // }));
        }

        // for handle in handles {
        //     handle.join().unwrap();
        // }
        
    }
}

impl fmt::Debug for Graph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let edges: Vec<&KmerEdge> = self.edges.iter().map(|x| unsafe{& *x.load(Acquire)}).collect();
        let vertices: Vec<&ProteinVertex> = self.vertices.iter().map(|x| unsafe{& *x.load(Acquire)}).collect();
        f.debug_struct("Graph")
            .field("Kmers", &edges)
            .field("Proteins", &vertices).finish()
    }
}
