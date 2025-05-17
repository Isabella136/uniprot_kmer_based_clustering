mod edge;
mod vertex;

use crate::protein::Protein;
use crate::graph::edge::KmerEdge;
use crate::graph::vertex::ProteinVertex;

use std::fmt;
use std::ptr;
use std::boxed::Box;
use std::thread::spawn;
use std::time::Instant;
use std::sync::{Arc, mpsc, Weak};
use std::sync::atomic::Ordering::{Acquire, AcqRel, Release};
use std::sync::atomic::{AtomicPtr, AtomicUsize, AtomicU64};

type ArcU64 = Arc<AtomicU64>;
type ArcVec<T> = Arc<Vec<T>>;
type ArcProtein = Arc<Protein>;
type ArcUsize = Arc<AtomicUsize>;
type ArcGraph = Arc<AtomicPtr<Graph>>;
type ArcKmerEdge = Arc<AtomicPtr<KmerEdge>>;
type ArcProteinVertex = Arc<AtomicPtr<ProteinVertex>>;

type WeakUsize = Weak<AtomicUsize>;
type WeakGraph = Weak<AtomicPtr<Graph>>;

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
        let edges_per_kmer: ArcVec<usize> = Arc::new((0..kmer_count)
            .map(|i| kmer_freq[i] * (kmer_freq[i] - 1) / 2).collect());
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
        let arc_kmer_freq: ArcVec<usize> = Arc::new(kmer_freq);

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
        let vertex_index: ArcUsize = Arc::new(AtomicUsize::new(0usize));
        let times_kmer_visited: ArcVec<ArcUsize> = Arc::new((0..kmer_count)
            .map(|_: usize| Arc::new(AtomicUsize::new(0usize))).collect());
        let mut handles  = vec![];
        for _ in 0..thread_count {
            let edge_count_prefix_sum: ArcVec<usize> = edge_count_prefix_sum.clone();
            let times_kmer_visited: ArcVec<ArcUsize> = times_kmer_visited.clone();
            let vertices: ArcVec<ArcProteinVertex> = graph.vertices.clone();
            let arc_protein_count: ArcUsize = arc_protein_count.clone();
            let arc_kmer_freq: ArcVec<usize> = arc_kmer_freq.clone();
            let vertex_index: ArcUsize = vertex_index.clone();
            let graph_arc: ArcGraph = graph_arc.clone();
            handles.push(spawn(move || {
                let protein_count: usize = arc_protein_count.load(Acquire);
                loop {
                    // let thread_now: Instant = Instant::now();
                    // If we already have enough vertices, break
                    let curr_protein = vertex_index.fetch_add(1, AcqRel);
                    if curr_protein >= protein_count {
                        break;
                    }

                    // Get atomic ptrs to vertex
                    let a_ptr_vertex: ArcProteinVertex = vertices[curr_protein].clone();

                    // Create new ProteinVertex
                    let mut curr_vertex: ProteinVertex = ProteinVertex::new(
                        curr_protein, Arc::downgrade(&graph_arc));
                    curr_vertex.update_graph_edges(edge_count_prefix_sum.clone(), 
                        times_kmer_visited.clone(), arc_kmer_freq.clone());

                    // Get raw ptrs to new KmerEdge and KmerEdgeHelper structs
                    let ptr_vertex = Box::into_raw(Box::new(curr_vertex));

                    // Store raw ptrs into atomic ptrs
                    a_ptr_vertex.store(ptr_vertex, Release);

                    // let elapsed = thread_now.elapsed().as_nanos();
                    // eprintln!("protein {} construction time: {} nanosecs", 
                    //     curr_protein, elapsed);
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

        let edge_indices: ArcUsize = Arc::new(AtomicUsize::new(0));
        let mut handles  = vec![];
        for _ in 0..thread_count {
            let edge_indices: ArcUsize = edge_indices.clone();
            let edges: ArcVec<ArcKmerEdge> = graph.edges.clone();
            let arc_edge_count: ArcUsize = arc_edge_count.clone();
            handles.push(spawn(move || {
                let edge_count = arc_edge_count.load(Acquire);
                loop {
                    // If we already have enough edges, break
                    let curr_edge_index = edge_indices.fetch_add(1, AcqRel);
                    if curr_edge_index >= edge_count {
                        break;
                    }

                    // No two threads will ever visit the same edge, so we are safe
                    let a_ptr_edge: ArcKmerEdge = edges[curr_edge_index].clone();
                    let ptr_edge: &KmerEdge = unsafe {& *a_ptr_edge.load(Acquire)};

                    // Check edge
                    let res = ptr_edge.double_checking();

                    if res.is_err() {
                        panic!("Oops!\n\tkmer: {:?},\tedge: {},\tvertices: {:?}",
                            ptr_edge.get_kmers(), curr_edge_index, ptr_edge.get_vertices_key());
                    }
                }
                
            }));
        }

        // Wait for threads to end
        for handle in handles {
            handle.join().unwrap();
        }

        graph_arc
    }

    // Combine edges with same two vertices
    pub fn combine_edges(&mut self, graph_weak: WeakGraph, thread_count: u32) {
        // Make ARC pointers of shared variables
        let edge_index: ArcUsize = Arc::new(AtomicUsize::new(0));
        let now = Instant::now();

        let (send, recv) = mpsc::channel::<(Vec<ArcUsize>, Vec<ArcKmerEdge>)>();
        let mut handles  = vec![]; 
        
        // Traverse through each edge and only keep unique pairs of vertices
        for _ in 0..thread_count {
            let send: mpsc::Sender<(Vec<ArcUsize>, Vec<ArcKmerEdge>)> = send.clone();
            let global_edge_keys: ArcVec<ArcUsize> = self.global_edge_keys.clone();
            let vertices: ArcVec<ArcProteinVertex> = self.vertices.clone();
            let edges: ArcVec<ArcKmerEdge> = self.edges.clone();
            let graph_weak: WeakGraph = graph_weak.clone();
            let edge_index: ArcUsize = edge_index.clone();
            handles.push(spawn(move || {
                let edge_count = edges.len();
                let mut unique_edges: Vec<ArcKmerEdge> = vec![];
                let mut unique_edge_keys: Vec<ArcUsize> = vec![];
                loop {
                    // If we already traversed through all edges, break
                    let edge_key = edge_index.fetch_add(1,  AcqRel);
                    if edge_key >= edge_count {
                        let _ = send.send((unique_edge_keys, unique_edges));
                        break;
                    }

                    // No one is updating edge, so we are safe
                    let edge_ptr: ArcKmerEdge = edges[edge_key].clone();
                    let edge_ref: &KmerEdge = unsafe {& *edge_ptr.load(Acquire)};
                    let vertices_keys: [usize; 2] = edge_ref.get_vertices_key();
                    
                    // Get edges containing first vertex
                    // No one is updating vertex, so we are safe
                    let first_vertex_ptr: ArcProteinVertex = 
                        vertices[vertices_keys[0]].clone();
                    let first_vertex: &ProteinVertex = unsafe {
                        & *first_vertex_ptr.load(Acquire)};
                    if first_vertex.get_edges_key().clone().len() == 1 {
                        unique_edge_keys.push(global_edge_keys[edge_key].clone());
                        unique_edges.push(edges[edge_key].clone());
                        if edge_key % 1000000 == 0 {
                            eprintln!("We made it to edge {}", edge_key)
                        }
                        continue;
                    }
                    let first_vertex_edges_bitarray: ArcVec<bool> = 
                        Arc::new(first_vertex.get_edges_bit_array(edges.len()));

                    // Get edges containing second vertex
                    // No one is updating vertex, so we are safe
                    let second_vertex_ptr: ArcProteinVertex = 
                        vertices[vertices_keys[1]].clone();
                    let second_vertex: &ProteinVertex = unsafe {
                        & *second_vertex_ptr.load(Acquire)};
                    let second_vertex_edges_keys: ArcVec<WeakUsize> = 
                        Arc::new(second_vertex.get_edges_key().clone());
                    if second_vertex_edges_keys.len() == 1 {
                        unique_edge_keys.push(global_edge_keys[edge_key].clone());
                        unique_edges.push(edges[edge_key].clone());
                        if edge_key % 1000000 == 0 {
                            eprintln!("We made it to edge {}", edge_key)
                        }
                        continue;
                    }
                    
                    let mut go_to_next = false;
                    let mut edges_to_merge_with: Vec<WeakUsize> = vec![];

                    for index in 0..second_vertex_edges_keys.len() {
                        let other_key = second_vertex_edges_keys[index].upgrade()
                            .unwrap().load(Acquire);
                        if first_vertex_edges_bitarray[other_key] {
                            if other_key < edge_key {
                                go_to_next = true;
                                break;
                            }
                            edges_to_merge_with.push(second_vertex_edges_keys[index].clone());
                        }
                    }

                    // If current edge share vertices with a previously
                    // visited edge, go to next edge
                    if go_to_next {
                        if edge_key % 1000000 == 0 {
                            eprintln!("We made it to edge {}", edge_key)
                        }
                        continue;
                    }

                    // Sort keys
                    edges_to_merge_with.sort_by(|a: &WeakUsize, b: &WeakUsize| 
                        a.upgrade().unwrap().load(Acquire)
                        .cmp(&b.upgrade().unwrap().load(Acquire)));

                    if edges_to_merge_with.len() == 1 {
                        unique_edge_keys.push(global_edge_keys[edge_key].clone());
                        unique_edges.push(edges[edge_key].clone());
                        if edge_key % 1000000 == 0 {
                            eprintln!("We made it to edge {}", edge_key)
                        }
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

                    unique_edge_keys.push(global_edge_keys[edge_key].clone());
                    unique_edges.push(edges[edge_key].clone());
                    if edge_key % 1000000 == 0 {
                        eprintln!("We made it to edge {}", edge_key)
                    }
                }
            }));
        }

        drop(send);

        let mut unique_edges: Vec<ArcKmerEdge> = vec![];
        let mut unique_edge_keys: Vec<ArcUsize> = vec![];

        for _ in 0..thread_count {
            let (mut unique_edge_keys_split, mut unique_edges_split) = 
                recv.recv().unwrap();
            unique_edges.append(&mut unique_edges_split);
            unique_edge_keys.append(&mut unique_edge_keys_split);
        }

        drop(recv);

        let elapsed = now.elapsed().as_nanos();
        eprintln!("edge traversal time: {} nanosecs",
            elapsed);


        // Remove weak pointers to non-unique edges
        let mut handles= vec![];
        let vertices_index: ArcUsize = Arc::new(AtomicUsize::new(0));
        let unique_edge_keys: ArcVec<ArcUsize> = Arc::new(unique_edge_keys);

        for _ in 0..thread_count {
            let vertices_index: ArcUsize = vertices_index.clone();
            let vertices: ArcVec<ArcProteinVertex> = self.vertices.clone();
            let unique_edge_keys: ArcVec<ArcUsize> = unique_edge_keys.clone();

            handles.push(spawn(move || {
                loop {
                    let curr_vertices_index = vertices_index.fetch_add(1, AcqRel);
                    if curr_vertices_index >= 10619 {
                        break;
                    }

                    // No two threads will ever traverse the same vertex, so we are safe
                    let curr_vertex: &mut ProteinVertex = unsafe {
                        &mut *vertices[curr_vertices_index].load(Acquire)};
                    curr_vertex.keep_specified_edges(unique_edge_keys.clone());
                }
            }))
        }

        for handle in handles {
            handle.join().unwrap();
        }

        let elapsed = now.elapsed().as_nanos();
        eprintln!("edge removal time: {} nanosecs",
            elapsed);


        // Update graph
        self.edges = Arc::new(unique_edges);
        self.global_edge_keys = unique_edge_keys;

        // Make ARC pointers of shared variables
        let edge_index: ArcUsize = Arc::new(AtomicUsize::new(0));

        let now = Instant::now();

        // Update keys in graph
        let mut handles= vec![];
        for _ in 0..thread_count {
            let global_edge_keys: ArcVec<ArcUsize> = self.global_edge_keys.clone();
            let edge_index: ArcUsize = edge_index.clone();
            handles.push(spawn(move || {
                let edge_count = global_edge_keys.len();
                loop {
                    let curr_index = edge_index.fetch_add(1, AcqRel);
                    if curr_index >= edge_count {
                        break;
                    }
                    
                    // No two threads will ever traverse the same edge, so we are safe
                    global_edge_keys[curr_index].store(curr_index, Release);

                    if curr_index % 1000000 == 0 {
                        eprintln!("We made it to edge {}", curr_index)
                    }
                    
                }
            }));
        }

        for handle in handles {
            handle.join().unwrap();
        }

        let elapsed = now.elapsed().as_nanos();
        eprintln!("key update time: {} nanosecs",
            elapsed);

        eprintln!("Number of edges now: {}", self.edges.len());
    }

    // Remove edges without diverging AMR labels
    pub fn remove_uninteresting_edges(&mut self, thread_count: u32) {
        // Make ARC pointers of shared variables
        let edge_index: ArcUsize = Arc::new(AtomicUsize::new(0));
        let now = Instant::now();

        // Only keep edges that represent two amr classes
        let mut handles= vec![];
        let (send, recv) = mpsc::channel::<(Vec<ArcUsize>, Vec<ArcKmerEdge>)>();
        for _ in 0..thread_count {
            let send: mpsc::Sender<(Vec<ArcUsize>, Vec<ArcKmerEdge>)> = send.clone();
            let global_edge_keys: ArcVec<ArcUsize> = self.global_edge_keys.clone();
            let proteins: ArcVec<ArcProtein> = self.protein_list.clone();
            let edges: ArcVec<ArcKmerEdge> = self.edges.clone();
            let edge_index: ArcUsize = edge_index.clone();
            handles.push(spawn(move || {
                let edge_count = edges.len();
                let mut two_amr_edges: Vec<ArcKmerEdge> = vec![];
                let mut two_amr_edge_keys: Vec<ArcUsize> = vec![];
                loop {
                    let curr_index = edge_index.fetch_add(1, AcqRel);
                    if curr_index >= edge_count {
                        let _ = send.send((two_amr_edge_keys, two_amr_edges));
                        break;
                    }
                    
                    // No two threads will ever traverse the same edge, so we are safe
                    let edge_ptr: ArcKmerEdge = edges[curr_index].clone();
                    let edge_ref: &KmerEdge = unsafe {& *edge_ptr.load(Acquire)};
                    let vertices_key: [usize; 2] = edge_ref.get_vertices_key();

                    let amr_class_0: &str = proteins[vertices_key[0]].get_amr_class();
                    let amr_class_1: &str = proteins[vertices_key[1]].get_amr_class();

                    // Can keep edge
                    if *amr_class_0 != *amr_class_1{
                        two_amr_edge_keys.push(global_edge_keys[curr_index].clone());
                        two_amr_edges.push(edges[curr_index].clone());
                    }

                    if curr_index % 1000000 == 0 {
                        eprintln!("We made it to edge {}", curr_index)
                    }
                    
                }
            }));
        }
        drop(send);

        for handle in handles {
            handle.join().unwrap();
        }

        let mut two_amr_edges: Vec<ArcKmerEdge> = vec![];
        let mut two_amr_edge_keys: Vec<ArcUsize> = vec![];

        for _ in 0..thread_count {
            let (mut two_amr_edge_keys_split, mut two_amr_edges_split) = 
                recv.recv().unwrap();
            two_amr_edges.append(&mut two_amr_edges_split);
            two_amr_edge_keys.append(&mut two_amr_edge_keys_split);
        }

        drop(recv);

        let elapsed = now.elapsed().as_nanos();
        eprintln!("edge traversal time: {} nanosecs",
            elapsed);


        // Remove weak pointers to one-amr edges
        let mut handles= vec![];
        let vertices_index: ArcUsize = Arc::new(AtomicUsize::new(0));
        let two_amr_edge_keys: ArcVec<ArcUsize> = Arc::new(two_amr_edge_keys);

        for _ in 0..thread_count {
            let vertices_index: ArcUsize = vertices_index.clone();
            let vertices: ArcVec<ArcProteinVertex> = self.vertices.clone();
            let two_amr_edge_keys: ArcVec<ArcUsize> = two_amr_edge_keys.clone();

            handles.push(spawn(move || {
                loop {
                    let curr_vertices_index = vertices_index.fetch_add(1, AcqRel);
                    if curr_vertices_index >= 10619 {
                        break;
                    }

                    // No two threads will ever traverse the same vertex, so we are safe
                    let curr_vertex: &mut ProteinVertex = unsafe {
                        &mut *vertices[curr_vertices_index].load(Acquire)};
                    curr_vertex.keep_specified_edges(two_amr_edge_keys.clone());
                }
            }))
        }

        for handle in handles {
            handle.join().unwrap();
        }

        let elapsed = now.elapsed().as_nanos();
        eprintln!("edge removal time: {} nanosecs",
            elapsed);


        // Update graph
        self.edges = Arc::new(two_amr_edges);
        self.global_edge_keys = two_amr_edge_keys;

        // Make ARC pointers of shared variables
        let edge_index: ArcUsize = Arc::new(AtomicUsize::new(0));

        let now = Instant::now();

        // Update keys in graph
        let mut handles= vec![];
        for _ in 0..thread_count {
            let global_edge_keys: ArcVec<ArcUsize> = self.global_edge_keys.clone();
            let edge_index: ArcUsize = edge_index.clone();
            handles.push(spawn(move || {
                let edge_count = global_edge_keys.len();
                loop {
                    let curr_index = edge_index.fetch_add(1, AcqRel);
                    if curr_index >= edge_count {
                        break;
                    }
                    
                    // No two threads will ever traverse the same edge, so we are safe
                    global_edge_keys[curr_index].store(curr_index, Release);

                    if curr_index % 1000000 == 0 {
                        eprintln!("We made it to edge {}", curr_index)
                    }
                    
                }
            }));
        }

        for handle in handles {
            handle.join().unwrap();
        }

        let elapsed = now.elapsed().as_nanos();
        eprintln!("key update time: {} nanosecs",
            elapsed);

        eprintln!("Number of edges now: {}", self.edges.len());
        
    }

    // pub fn get_divergent_protein_pairs(&self, thread_count: u32) {
    //     // Get bitarray representing vertex pairs we've already visited
    //     let vertex_count = self.vertices.len();
    //     let pair_count = vertex_count*(vertex_count-1)/2;

    //     let bitarray_edge_visited: ArcVec<ArcBool> = Arc::new((0..self.edges.len())
    //         .into_iter().map(|_| Arc::new(AtomicBool::new(false))).collect());
    //     // let divergent_pair: ArcVec<ArcKmerEdge> = Arc::new((0..pair_count)
    //     //     .into_iter().map(|_| Arc::new(AtomicBool::new(false))).collect());

    //     // Make ARC pointers of shared variables
    //     let edge_index: ArcUsize = Arc::new(AtomicUsize::new(0));
    //     // let mut handles= vec![];

    //     // Find protein pairs
    //     for _ in 0..thread_count {
    //         // let divergent_pair: ArcVec<ArcKmerEdge> = divergent_pair.clone();
    //         let bitarray_edge_visited: ArcVec<ArcBool> = bitarray_edge_visited.clone();
    //         let proteins: ArcVec<ArcProtein> = self.protein_list.clone();
    //         let edges: ArcVec<ArcKmerEdge> = self.edges.clone();
    //         let edge_index: ArcUsize = edge_index.clone();
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
    //     }

    //     // for handle in handles {
    //     //     handle.join().unwrap();
    //     // }
        
    // }
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
