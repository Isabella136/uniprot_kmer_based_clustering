#![feature(get_mut_unchecked)]

use std::sync::atomic::Ordering::Acquire;
use seq_io::parallel::parallel_fasta;
use std::sync::{Arc, Mutex, RwLock};
use seq_io::fasta::Reader;
use std::time::Instant;
use std::ops::Deref;
use std::thread;
use boomphf::*;
use std::env;

mod protein;
mod graph;
// mod tree;

use crate::protein::Protein;
use crate::graph::Graph;

type FiveMer = u32;


fn merge_sort(list: &mut Vec<(u32, u32)>, start: usize, end: usize, item: u32) {
    //println!("List: {list:#?}, Start: {start}, End: {end}, Item: {item}");
    if end - start <= 1 {
        if item < list[start].0 {
            list.insert(start, (item, 1));
        }
        else if item > list[start].0 {
            list.insert(start + 1usize, (item, 1));
        }
        else {
            list[start] = (list[start].0, list[start].1 + 1);
        }
    }
    else {
        let mid = (end+start) / 2;
        if item > list[mid].0 {
            merge_sort(list, mid, end, item);
        }
        else if item < list[mid].0 {
            merge_sort(list, start, mid, item);
        }
        else {
            list[mid] = (list[mid].0, list[mid].1 + 1);
        }
    }
}

fn main() {
    eprintln!("We start main");
    env::set_var("RUST_BACKTRACE", "1");
    // Arguments required: input, threads
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        panic!("Requires two command line arguments: input and thread");
    }
    let input: &String = &args[1];
    let threads: u32 = args[2].parse()
        .expect("threads argument should be of type int");

    let reader = Reader::from_path(input)
        .expect("input argument should refer to an existing fasta file");
    let protein_list: Arc<Mutex<Vec<Protein>>> = Arc::new(Mutex::new(Vec::new()));
    parallel_fasta(reader, threads, 64,
        |record, found|{
            *found = Protein::new(&record);
            let mut protein_list = protein_list.lock().unwrap();
            protein_list.push(found.to_owned())
        }, |_record, _found|{
            Some(())
        }).unwrap();

    eprintln!("We created Protein structs");

    // Each protein should now have a list of random 5-mers
    let five_mer_freq_list: Arc<Mutex<Vec<(FiveMer, u32)>>> = Arc::new(
        Mutex::new(Vec::new()));
    let next_protein_index = Arc::new(Mutex::new(0usize));
    let protein_list = Arc::new(protein_list.lock().unwrap().clone());
    let mut handles = vec![];

    // We'll try to combine all 5-mers
    for _ in 0..threads {
        let protein_list = protein_list.clone();
        let five_mer_freq_list = five_mer_freq_list.clone();
        let next_protein_index = next_protein_index.clone();
        handles.push(thread::spawn(move || {
            loop {
                let curr_protein_index = {
                    let mut next_protein_index = next_protein_index.lock().unwrap();
                    let curr_protein_index = *next_protein_index;
                    *next_protein_index += 1;
                    curr_protein_index
                };
                if curr_protein_index >= 10619 {
                    break;
                }
                let curr_protein = &protein_list[curr_protein_index];
                let mut five_mers = curr_protein.get_five_mers();
                five_mers.sort();
                five_mers.dedup();
                let _ = {
                    let mut five_mer_freq_list = five_mer_freq_list.lock().unwrap();
                    let list_len = five_mer_freq_list.len();
                    if list_len == 0 {
                        *five_mer_freq_list = five_mers.iter().map(|x| (*x, 1)).collect();
                    }
                    else {
                        for item in five_mers {
                            let list_len = five_mer_freq_list.len();
                            merge_sort(&mut five_mer_freq_list, 0, list_len, item);
                        }
                    }

                };
            }
        }))
    }
    for handle in handles {
        handle.join().unwrap();
    }

    eprintln!("We combined k-mers");

    // Now we can remove unique kmers and develop the minimal perfect hash functions
    let mut five_mer_unique = vec![];
    let mut five_mer_repeat = vec![];
    let five_mer_freq_list = five_mer_freq_list.lock().unwrap();
    for index in 0..five_mer_freq_list.len() {
        if five_mer_freq_list[index].1 == 1 {
            five_mer_unique.push(five_mer_freq_list[index].0);
        }
        else {
            five_mer_repeat.push(five_mer_freq_list[index].0);
        }
    }
    let five_mer_all: Vec<FiveMer> = (*five_mer_freq_list).iter().map(|(x,_y)| *x).collect();
    let five_mer_all_phf = Arc::new(Mphf::new(3.0, &five_mer_all));
    let five_mer_repeat_phf = Arc::new(Mphf::new(3.0, &five_mer_repeat));
    let five_mer_all_len = five_mer_all.len();
    let five_mer_repeat_len = Arc::new(five_mer_repeat.len());
    let mut five_mer_unique_hashmap = vec![false; five_mer_all_len];
    for five_mer in five_mer_unique {
        let curr_hash = five_mer_all_phf.hash(&five_mer) as usize;
        five_mer_unique_hashmap[curr_hash] = true;
    }
    drop(five_mer_repeat);
    drop(five_mer_freq_list);

    eprintln!("We found unique k-mers");

    let mut handles = vec![];
    let five_mer_hash_freq = Arc::new(Mutex::new(vec![0usize; *five_mer_repeat_len]));
    let five_mer_unique_hashmap = Arc::new(five_mer_unique_hashmap);
    let next_protein_index = Arc::new(Mutex::new(0usize));
    let protein_list: Arc<Vec<Arc<RwLock<Protein>>>> = Arc::new((*protein_list).iter()
        .map(|x| Arc::new(RwLock::new(x.clone()))).collect());
    for _ in 0..threads {
        let five_mer_all_phf = five_mer_all_phf.clone();
        let five_mer_repeat_phf = five_mer_repeat_phf.clone();
        let five_mer_repeat_len = five_mer_repeat_len.clone();
        let five_mer_unique_hashmap = five_mer_unique_hashmap.clone();

        let protein_list = protein_list.clone();
        let five_mer_hash_freq = five_mer_hash_freq.clone();
        let next_protein_index = next_protein_index.clone();

        handles.push(thread::spawn(move || {
            loop {
                let curr_protein_index = {
                    let mut next_protein_index = next_protein_index.lock().unwrap();
                    let curr_protein_index = *next_protein_index;
                    *next_protein_index += 1;
                    curr_protein_index
                };
                if curr_protein_index >= 10619 {
                    break;
                }
                let curr_protein = &mut protein_list[curr_protein_index]
                    .write().unwrap();
                curr_protein.remove_unique_five_mers(
                    &five_mer_all_phf, five_mer_unique_hashmap.deref());
                curr_protein.modify_hash_five_mer(
                    five_mer_repeat_len.deref(), &five_mer_repeat_phf);

                let mut five_mers: Vec<u32> = curr_protein.get_five_mers();
                five_mers.sort();
                five_mers.dedup();
                let mut five_mer_hash_freq = five_mer_hash_freq.lock().unwrap();
                for five_mer in five_mers {
                    (*five_mer_hash_freq)[five_mer_repeat_phf.hash(&five_mer) as usize] += 1;
                }
            }
        }))
    }
    for handle in handles {
        handle.join().unwrap();
    }
    drop(five_mer_all_phf);
    eprintln!("We made unique hash");

    let five_mer_hash_freq = Arc::into_inner(five_mer_hash_freq).unwrap().into_inner().unwrap();
    let protein_list: Arc<Vec<Arc<Protein>>> = Arc::new({
        let mut temp_list = vec![];
        let mut protein_list = protein_list;
        let protein_list = Arc::get_mut(&mut protein_list).unwrap();
        while protein_list.len() > 0 {
            let protein = Arc::into_inner(protein_list.remove(0)).unwrap();
            temp_list.push(Arc::new(protein.into_inner().unwrap()));
        }
        temp_list});

    eprintln!("We can make a graph");

    let now = Instant::now();
    let graph_arc = Graph::new(
        five_mer_hash_freq, threads as usize, protein_list.clone());
    let elapsed = now.elapsed().as_secs_f64();
    eprintln!("Graph construction time: {} seconds", 
        elapsed);
    let now = Instant::now();
    let mut graph_ref = unsafe {&mut *graph_arc.load(Acquire)};
    graph_ref.remove_uninteresting_edges(threads);
    graph_ref = unsafe {&mut *graph_arc.load(Acquire)};
    graph_ref.combine_edges(Arc::downgrade(&graph_arc), threads);
    graph_ref = unsafe {&mut *graph_arc.load(Acquire)};
    let elapsed = now.elapsed().as_secs_f64();
    eprintln!("Graph refinement time: {} seconds", 
        elapsed);

    graph_ref.align_and_output_pairs(threads);
    

    println!("Graph right now:\n{graph_ref:#?}");


    
}
