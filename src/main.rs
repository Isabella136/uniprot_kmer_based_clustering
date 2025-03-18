use seq_io::parallel::parallel_fasta;
use std::sync::{Arc, Mutex};
use seq_io::fasta::Reader;
use std::thread;
use boomphf::*;
use std::env;

mod protein;
mod tree;

use crate::protein::Protein;

type FiveMer = u32;
type SevenMer = u32;


fn merge_sort(list: &mut Vec<u32>, start: usize, end: usize, item: u32) {
    //println!("List: {list:#?}, Start: {start}, End: {end}, Item: {item}");
    if end - start <= 1 {
        if item < list[start] {
            list.insert(start, item);
        }
        else if item > list[start] {
            list.insert(start + 1usize, item);
        }
    }
    else {
        let mid = (end+start) / 2;
        if item > list[mid] {
            merge_sort(list, mid, end, item);
        }
        else if item < list[mid] {
            merge_sort(list, start, mid, item);
        }
    }
}

fn main() {
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
            *found = Protein::new_protein(&record);
            let mut protein_list = protein_list.lock().unwrap();
            protein_list.push(found.to_owned())
        }, |_record, _found|{
            Some(())
        }).unwrap();

    // Each protein should now have a list of random 5-mers and 7-mers
    let five_mer_list: Arc<Mutex<Vec<FiveMer>>> = Arc::new(
        Mutex::new(Vec::new()));
    let seven_mer_list: Arc<Mutex<Vec<SevenMer>>> = Arc::new(
        Mutex::new(Vec::new()));
    let next_protein_index = Arc::new(Mutex::new(0usize));
    let protein_list = Arc::new(protein_list.lock().unwrap().clone());
    let mut handles = vec![];

    // We'll try to combine all 5-mers and all 7-mers
    for _ in 0..threads {
        let protein_list = Arc::clone(&protein_list);
        let five_mer_list = Arc::clone(&five_mer_list);
        let seven_mer_list = Arc::clone(&seven_mer_list);
        let next_protein_index = Arc::clone(&next_protein_index);
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
                //println!("Current Protein: {curr_protein_index}");
                let curr_protein = &protein_list[curr_protein_index];
                let mut five_mers = curr_protein.get_five_mers();
                let mut seven_mers = curr_protein.get_seven_mers();
                five_mers.sort();
                seven_mers.sort();
                let _ = {
                    let mut five_mer_list = five_mer_list.lock().unwrap();
                    let list_len = five_mer_list.len();
                    if list_len == 0 {
                        *five_mer_list = five_mers
                    }
                    else {
                        for item in five_mers {
                            //println!("No stack overflow with 5-mer");
                            let list_len = five_mer_list.len();
                            //println!("List length: {list_len}");
                            merge_sort(&mut five_mer_list, 0, list_len, item);
                        }
                    }

                };
                let _ = {
                    let mut seven_mer_list = seven_mer_list.lock().unwrap();
                    let list_len = seven_mer_list.len();
                    if list_len == 0 {
                        *seven_mer_list = seven_mers
                    }
                    else {
                        for item in seven_mers {
                            //println!("No stack overflow with 7-mer");
                            let list_len = seven_mer_list.len();
                            //println!("List length: {list_len}");
                            merge_sort(&mut seven_mer_list, 0, list_len, item);
                        }
                    }
                };
            }
        }))
    }
    for handle in handles {
        handle.join().unwrap();
    }
    
    // Now we can develop the minimal perfect hash functions
    let five_mer_list = five_mer_list.lock().unwrap();
    let five_phf = Arc::new(Mphf::new(3.0, &five_mer_list));
    let five_mer_list_len = Arc::new(five_mer_list.len());

    let seven_mer_list = seven_mer_list.lock().unwrap();
    let seven_phf = Arc::new(Mphf::new(3.0, &seven_mer_list));
    let seven_mer_list_len = Arc::new(seven_mer_list.len());
    
    drop(five_mer_list);
    drop(seven_mer_list);

    // We can now populate hash maps for each protein
    let mut handles = vec![];
    let next_protein_index = Arc::new(Mutex::new(0usize));
    let protein_list: Arc<Vec<Mutex<Protein>>> = Arc::new((*protein_list).iter()
        .map(|x| Mutex::new(x.clone())).collect());
    for _ in 0..threads {
        let five_phf = Arc::clone(&five_phf);
        let seven_phf = Arc::clone(&seven_phf);
        let five_mer_list_len = Arc::clone(&five_mer_list_len);
        let seven_mer_list_len = Arc::clone(&seven_mer_list_len);
        let protein_list = Arc::clone(&protein_list);
        let next_protein_index = Arc::clone(&next_protein_index);

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
                    .lock().unwrap();
                curr_protein.modify_hash_five_mer(&five_mer_list_len, &five_phf);
                curr_protein.modify_hash_seven_mer(&seven_mer_list_len, &seven_phf);
                
            }
        }))
    }
    for handle in handles {
        handle.join().unwrap();
    }

    // Goal is to take all sequence and find all top ten 5-mers
    // Use those 5-mers to come up with a minimal perfect hash
    // Upper bound for number of total top 10 5-mers: 10 * 10619


    
}
