use seq_io::parallel::parallel_fasta;
use std::sync::{Arc, Mutex, RwLock};
use seq_io::fasta::Reader;
use std::thread;
use boomphf::*;
use std::time;
use std::env;

mod protein;
mod tree;

use crate::protein::Protein;
use crate::tree::Tree;

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
        let protein_list = protein_list.clone();
        let five_mer_list = five_mer_list.clone();
        let seven_mer_list = seven_mer_list.clone();
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
    let five_mer_freq = Arc::new(Mutex::new(vec![0usize; *five_mer_list_len]));
    let next_protein_index = Arc::new(Mutex::new(0usize));
    let protein_list: Arc<Vec<Arc<RwLock<Protein>>>> = Arc::new((*protein_list).iter()
        .map(|x| Arc::new(RwLock::new(x.clone()))).collect());
    for _ in 0..threads {
        let five_phf = five_phf.clone();
        let seven_phf = seven_phf.clone();
        let five_mer_freq = five_mer_freq.clone();
        let five_mer_list_len = five_mer_list_len.clone();
        let seven_mer_list_len = seven_mer_list_len.clone();
        let protein_list = protein_list.clone();
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
                curr_protein.modify_hash_five_mer(&five_mer_list_len, &five_phf);
                curr_protein.modify_hash_seven_mer(&seven_mer_list_len, &seven_phf);

                let five_mers = curr_protein.get_five_mers();
                let mut five_mer_freq = five_mer_freq.lock().unwrap();
                for five_mer in five_mers {
                    (*five_mer_freq)[five_phf.hash(&five_mer) as usize] += 1;
                }
            }
        }))
    }
    for handle in handles {
        handle.join().unwrap();
    }

    let five_mer_freq = Arc::into_inner(five_mer_freq).unwrap().into_inner().unwrap();
    let five_mer_freq_enum: Vec<(usize, &usize)> = five_mer_freq.iter().enumerate().collect();
    let max_five_mer = five_mer_freq_enum.iter().max_by(
        |x, y| x.1.cmp(y.1)).unwrap();
    eprintln!("{max_five_mer:#?}");

    let next_protein_index = Arc::new(Mutex::new(1usize));
    let protein_list: Arc<Vec<Arc<Protein>>> = Arc::new({
        let mut temp_list = vec![];
        let mut protein_list = protein_list;
        let protein_list = Arc::get_mut(&mut protein_list).unwrap();
        while protein_list.len() > 0 {
            let protein = Arc::into_inner(protein_list.remove(0)).unwrap();
            temp_list.push(Arc::new(protein.into_inner().unwrap()));
        }
        temp_list});
    let tree = Arc::new(Tree::new(5u8, protein_list[0].clone()));
    let mut handles = vec![];
    let now = time::Instant::now();

    for i in 0..threads {
        let tree = tree.clone();
        let protein_list = protein_list.clone();
        let next_protein_index = next_protein_index.clone();

        handles.push(thread::spawn(move || {
            loop {
                let now_local = time::Instant::now();

                let curr_protein_index = {
                    let mut next_protein_index = next_protein_index.lock().unwrap();
                    let curr_protein_index = *next_protein_index;
                    *next_protein_index += 1;
                    curr_protein_index
                };
                if curr_protein_index >= 10619 {
                    //println!("Done");
                    break;
                }
                let curr_protein = protein_list[curr_protein_index].clone();
                tree.add_protein(curr_protein);
                
                eprintln!("Thread {i} - Current Protein: {curr_protein_index} - Time: {} secs", now_local.elapsed().as_secs());
            }
        }))
    }
    for handle in handles {
        handle.join().unwrap();
    }
    println!("Tree building took {} seconds", now.elapsed().as_secs());

    print!("{tree:#?}");

    // Goal is to take all sequence and find all top ten 5-mers
    // Use those 5-mers to come up with a minimal perfect hash
    // Upper bound for number of total top 10 5-mers: 10 * 10619


    
}
