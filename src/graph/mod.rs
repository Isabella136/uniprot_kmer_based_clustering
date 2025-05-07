mod edge;
mod vertex;

use crate::protein::Protein;
use crate::graph::edge::KmerEdge;
use crate::graph::edge::KmerEdgeHelper;
use crate::graph::vertex::ProteinVertex;

use std::fmt;
use std::rc::Rc;
use std::thread::spawn;
use std::sync::{Arc, mpsc, Weak};
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::atomic::Ordering::{Acquire, AcqRel, Release};

type ArcGraph = Arc<Graph>;
type ArcVec<T> = Arc<Vec<T>>;
type ArcProtein = Arc<Protein>;
type ArcUsize = Arc<AtomicUsize>;
type ArcKmerEdge = Arc<KmerEdge>;
type ArcProteinVertex = Arc<ProteinVertex>;
type ArcKmerEdgeHelper = Arc<KmerEdgeHelper>;

type WeakUsize = Weak<AtomicUsize>;

type OptionWeakUsize = Rc<Option<WeakUsize>>;

pub struct Graph {
    pub edges: ArcVec<ArcKmerEdge>,
    protein_list: ArcVec<ArcProtein>,
    vertices: ArcVec<ArcProteinVertex>,
    global_edge_keys: ArcVec<ArcUsize>,
    edge_update_helpers : ArcVec<ArcKmerEdgeHelper>,
}
impl Graph {
    // Create a graph and return an ARC pointer to said graph
    pub fn new(kmer_count: usize, thread_count: usize, 
            protein_list: ArcVec<ArcProtein>) -> ArcGraph {

        eprintln!("We made it to new");

        // Create Graph object skeleton
        let global_edge_keys: ArcVec<ArcUsize> = Arc::new((0..kmer_count).into_iter()
            .map(|k: usize| Arc::new(AtomicUsize::new(k))).collect());
        let mut edge_update_helpers: ArcVec<ArcKmerEdgeHelper> = Arc::new(
            Vec::with_capacity(kmer_count));
        let mut vertices: ArcVec<ArcProteinVertex> = Arc::new(
            Vec::with_capacity(protein_list.len()));
        let mut edges: ArcVec<ArcKmerEdge> = Arc::new(
            Vec::with_capacity(kmer_count));

        // Make ARC pointers of shared variables 
        let arc_kmer_count: ArcUsize = Arc::new(AtomicUsize::new(kmer_count));
        let arc_protein_count: ArcUsize = Arc::new(AtomicUsize::new(protein_list.len()));

        // Concurrently make new KmerEdge structs
        let total_edges: ArcUsize = Arc::new(AtomicUsize::new(0usize));
        let (send, recv) = mpsc::channel::<
            (ArcUsize, ArcVec<ArcKmerEdge>)>();
        let mut handles  = vec![];
        for _ in 0..thread_count {
            let total_edges: ArcUsize = total_edges.clone();
            let arc_kmer_count: ArcUsize = arc_kmer_count.clone();
            let arc_protein_count: ArcUsize = arc_protein_count.clone();
            let prefix_length: ArcUsize = Arc::new(AtomicUsize::new(0));
            let send: mpsc::Sender<(ArcUsize, ArcVec<ArcKmerEdge>)> = send.clone();
            let mut split_edges: ArcVec<ArcKmerEdge> = 
                Arc::new(Vec::with_capacity(kmer_count/(thread_count/2)));
            handles.push(spawn(move || {
                loop {

                    // If we already have enough edges, break
                    let curr_edge_count = total_edges.fetch_add(1, AcqRel);
                    let kmer_count = arc_kmer_count.load(Acquire);
                    if curr_edge_count >= kmer_count {
                        let _ = send.send((prefix_length, split_edges));
                        break;
                    }

                    // Create new edge
                    let curr_edge: ArcKmerEdge = Arc::new(KmerEdge::new(
                        prefix_length.clone(), split_edges.len(),
                        arc_protein_count.load(Acquire)));   

                    Arc::get_mut(&mut split_edges).unwrap().push(curr_edge.clone());
                }
                
            }));
        }

        drop(send);

        // Wait for threads to end
        for handle in handles {
            handle.join().unwrap();
        }

        // Aggregrate KmerEdge structs into one vector
        for _ in 0..thread_count {
            let (prefix_length, mut split_edges) = recv.recv().unwrap();
            prefix_length.store(edges.len(), Release);
            Arc::get_mut(&mut edges).unwrap().append(
                Arc::get_mut(&mut split_edges).unwrap());
        }
        
        drop(recv);
        eprintln!("We created KmerEdge structs");

        // Concurrently make new KmerEdge update helpers
        let total_edge_helpers: ArcUsize = Arc::new(AtomicUsize::new(0usize));
        let (send, recv) = mpsc::channel::<ArcVec<ArcKmerEdgeHelper>>();
        let mut handles  = vec![];
        for _ in 0..thread_count {
            let mut split_helpers: ArcVec<ArcKmerEdgeHelper> = 
                Arc::new(Vec::with_capacity(kmer_count/(thread_count/2)));
            let send: mpsc::Sender<ArcVec<ArcKmerEdgeHelper>> = send.clone();
            let total_edge_helpers: ArcUsize = total_edge_helpers.clone();
            let arc_kmer_count: ArcUsize = arc_kmer_count.clone();
            handles.push(spawn(move || {
                loop {

                    // If we already have enough helpers, break
                    let curr_helper_count = total_edge_helpers.fetch_add(1, AcqRel);
                    let kmer_count = arc_kmer_count.load(Acquire);
                    if curr_helper_count >= kmer_count {
                        let _ = send.send(split_helpers);
                        break;
                    }
                    // Create new helper
                    let curr_helper: ArcKmerEdgeHelper = Arc::new(KmerEdgeHelper::new());

                    Arc::get_mut(&mut split_helpers).unwrap().push(curr_helper.clone());
                }
            }));
        }

        drop(send);

        // Wait for threads to end
        for handle in handles {
            handle.join().unwrap();
        }

        // Aggregrate KmerEdgeHelper structs into one vector
        for _ in 0..thread_count {
            let mut split_helpers: ArcVec<ArcKmerEdgeHelper> = recv.recv().unwrap();

            // We are in single-threaded environment, so we are safe
            Arc::get_mut(&mut edge_update_helpers).unwrap().append(
                Arc::get_mut(&mut split_helpers).unwrap());
        }
        
        drop(recv);
        eprintln!("We created KmerEdgeHelper structs");

        // Construct graph
        let graph: Graph = Graph {
            vertices: vertices.clone(),
            edges,
            protein_list,
            global_edge_keys,
            edge_update_helpers,
        };

        eprintln!("Graph is made");
        
        let graph_arc: ArcGraph = Arc::new(graph);

        // Concurrently make new ProteinVertex structs
        let (send, recv) = mpsc::channel::<ArcVec<ArcProteinVertex>>();
        let total_vertices: ArcUsize = Arc::new(AtomicUsize::new(0usize));
        let mut handles  = vec![];
        for _ in 0..thread_count {
            let mut split_vertices: ArcVec<ArcProteinVertex> = 
                Arc::new(Vec::with_capacity(arc_protein_count.load(Acquire)));
            let send: mpsc::Sender<ArcVec<ArcProteinVertex>> = send.clone();
            let arc_protein_count: ArcUsize = arc_protein_count.clone();
            let total_vertices: ArcUsize = total_vertices.clone();
            let graph_arc: ArcGraph = graph_arc.clone();
            handles.push(spawn(move || {
                loop {

                    // If we already have enough vertices, break
                    let curr_vertex_count = total_vertices.fetch_add(1, AcqRel);
                    let protein_count: usize = arc_protein_count.load(Acquire);
                    if curr_vertex_count >= protein_count {
                        let _ = send.send(split_vertices);
                        break;
                    }

                    // Create new vertex
                    let curr_vertex: ArcProteinVertex = Arc::new(ProteinVertex::new(
                        Arc::downgrade(&graph_arc)));

                    // One thread per vector, so we are safe
                    Arc::get_mut(&mut split_vertices).unwrap().push(curr_vertex);
                }
            }));
        }

        drop(send);

        // Wait for threads to end
        for handle in handles {
            handle.join().unwrap();
        }

        // Aggregrate ProteinVertex structs into one vector
        for _ in 0..thread_count {
            let mut split_vertices: ArcVec<ArcProteinVertex> = recv.recv().unwrap();

            // We are in single-threaded environment, so we are safe
            unsafe {Arc::get_mut_unchecked(&mut vertices).append(
                Arc::get_mut(&mut split_vertices).unwrap())};
        }

        drop(recv);
        eprintln!("We created empty ProteinVertex structs");

        // Concurrently update ProteinVertex and KmerEdge structs
        let vertex_index: ArcUsize = Arc::new(AtomicUsize::new(0usize));
        let mut handles  = vec![];
        for _ in 0..thread_count {
            let arc_protein_count: ArcUsize = arc_protein_count.clone();
            let vertices: ArcVec<ArcProteinVertex> = vertices.clone();
            let vertex_index: ArcUsize = vertex_index.clone();
            handles.push(spawn(move || {
                loop {

                    // If we already have enough vertices, break
                    let curr_vertex_index = vertex_index.fetch_add(1, AcqRel);
                    let protein_count: usize = arc_protein_count.load(Acquire);
                    if curr_vertex_index >= protein_count {
                        break;
                    }

                    // Get current vertex
                    let mut curr_vertex: Arc<ProteinVertex> = vertices[curr_vertex_index].clone();

                    // Update key info in new vertex to be equal to its index in vector
                    // Only current thread will access this vertex, so we are safe
                    unsafe { Arc::get_mut_unchecked(&mut curr_vertex)
                        .update(curr_vertex_index);
                    }

                    // Update KmerEdge structs corresponding to edges new vertex is in 
                    // to reflect vertex's presence
                    curr_vertex.update_graph_edges();
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
                graph_arc.edge_update_helpers.clone();
            let edges: ArcVec<ArcKmerEdge> = graph_arc.edges.clone();
            let arc_kmer_count: ArcUsize = arc_kmer_count.clone();
            let edge_index: ArcUsize = edge_index.clone();
            
            handles.push(spawn(move || {
                loop {

                    // If we have gone through all edges, break
                    let curr_index = edge_index.fetch_add(1, AcqRel);
                    let kmer_count = arc_kmer_count.load(Acquire);
                    if curr_index >= kmer_count {
                        break;
                    }

                    // Retrieve current helper
                    let curr_helper: ArcKmerEdgeHelper = edge_helpers[curr_index].clone();

                    // Retrieve current edge
                    let curr_edge: &mut ArcKmerEdge = &mut edges[curr_index].clone();

                    // Only one thread can access a specific edge at a time, so we are safe
                    let curr_edge_raw: &mut KmerEdge = unsafe {
                        Arc::get_mut_unchecked(curr_edge)};

                    curr_helper.update_kmer_edge_final(curr_edge_raw);
                }
            }));
        }

        // Wait for threads to end
        for handle in handles {
            handle.join().unwrap();
        }

        // Currently in a single threaded environment, so we are safe
        //unsafe {Arc::get_mut_unchecked(&mut graph_arc).vertices = vertices};
        graph_arc
    }

    // Combine edges with lots of overlap
    pub fn combine_edges(mut graph: ArcGraph, thread_count: u32) {
        // Get current keys to edges
        let mut unordered_edge_keys: Vec<OptionWeakUsize> = (0..graph.edges.len())
            .into_iter().map(|i| Rc::new(Some(Arc::downgrade(
            &graph.global_edge_keys[i])))).collect();

        // Make a second vector of current keys, ordered by edge size
        let mut ordered_edge_keys: Vec<OptionWeakUsize> = unordered_edge_keys.clone();
        ordered_edge_keys.sort_by(|a: &OptionWeakUsize, b: &OptionWeakUsize| {
            let cloned_a: Option<WeakUsize> = a.as_ref().clone();
            let cloned_b: Option<WeakUsize> = b.as_ref().clone();
            let key_a: usize = cloned_a.unwrap().upgrade().unwrap().load(Acquire);
            let key_b: usize = cloned_b.unwrap().upgrade().unwrap().load(Acquire);
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
            let unwraped_edge_key: WeakUsize = edge_key.as_ref().clone().unwrap();
            let loaded_edge_key: usize = unwraped_edge_key.upgrade().unwrap().load(Acquire);
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
            let mut other_edges_keys: Vec<WeakUsize> = 
                first_vertex.get_edges_key().clone();

            // Make channels and thread handles
            let (send, recv) = mpsc::channel::<ArcVec<WeakUsize>>();
            let mut handles  = vec![];

            // Retrieve other edges that share a vertex with current edge
            for _ in 0..thread_count {
                let graph: ArcGraph = graph.clone();
                let send: mpsc::Sender<ArcVec<WeakUsize>> = send.clone();
                let vertices_index: ArcUsize = vertices_index.clone();
                let vertices_keys: ArcVec<usize> = vertices_keys.clone();
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
                        let curr_vertex: &ArcProteinVertex = 
                            &graph.vertices[curr_vertex_key];

                        let more_edges_keys: &Vec<WeakUsize> = &curr_vertex.get_edges_key();

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
                let graph: ArcGraph = graph.clone();
                let edge: ArcKmerEdge = edge.clone();
                let send: mpsc::Sender<ArcVec<WeakUsize>> = send.clone();
                let other_edges_index: ArcUsize = other_edges_index.clone();
                let other_edges_keys: ArcVec<WeakUsize> = other_edges_keys.clone();
                let mut edges_to_merge_with_split: ArcVec<WeakUsize> = Arc::new(vec![]);

                handles.push(spawn(move || {
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
                        let other_edge: &ArcKmerEdge = &graph.edges[loaded_other_key];

                        let other_cmp_edge: f64 = edge.compare(other_edge);
                        let edge_cmp_other: f64 = other_edge.compare(&edge);
                        if other_cmp_edge == 1.0 {
                            Arc::get_mut(&mut edges_to_merge_with_split)
                                .unwrap().push(other_key.clone());
                        }
                        else if edge_cmp_other >= 0.8 && other_cmp_edge >= 0.8 {
                            Arc::get_mut(&mut edges_to_merge_with_split)
                                .unwrap().push(other_key.clone());
                        }
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
                continue;
            }

            // Merge edges into one
            let kmer_edge_group: KmerEdge = KmerEdge::merge(
                &edges_to_merge_with, Arc::downgrade(&graph));
            kmer_edge_group.update_protein_vertices(unwraped_edge_key, thread_count);

            // Currently in a single threaded environment, so we are safe
            unsafe {
                let graph_mut = Arc::get_mut_unchecked(&mut graph);
                Arc::get_mut_unchecked(&mut graph_mut.edges)[loaded_edge_key] = 
                    Arc::new(kmer_edge_group)
            };
            
            for old_key in &*edges_to_merge_with {
                let loaded = old_key.upgrade().unwrap().load(Acquire);
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
        let vertices_index: ArcUsize = Arc::new(AtomicUsize::new(0));
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

        for handle in handles {
            handle.join().unwrap();
        }
    }

    fn remove_edges_marked_for_deletions(
            mut graph: ArcGraph, unordered_edge_keys: Vec<OptionWeakUsize>) -> Vec<bool>{

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
