use std::ops::Deref;
//use std::sync::{Arc, RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::sync::{Arc, RwLock, RwLockWriteGuard};
use std::{fmt, vec};

use crate::protein::Protein;

fn unionize_bitarrays(
        mut bit_array: Vec<bool>, mut bit_array_indices: Vec<u32>, 
        array_of_bit_array_indices: Vec<Arc<Vec<u32>>>) -> (Vec<bool>, Vec<u32>) {
    for other_bit_array_indices in array_of_bit_array_indices {
        for other_bit_index in other_bit_array_indices.deref() {
            if !bit_array[*other_bit_index as usize] {
                bit_array[*other_bit_index as usize] = true;
                bit_array_indices.push(*other_bit_index);
            }
        }
    }
    return (bit_array, bit_array_indices);
}
fn intersect_bitarrays(
    mut bit_array: Vec<bool>, mut bit_array_indices: Vec<u32>, 
    array_of_bit_arrays: Vec<Arc<Vec<bool>>>) -> Vec<u32> {
    
    for index in bit_array_indices.iter() {
        for other_bit_array in array_of_bit_arrays.iter() {
            if !other_bit_array.deref()[*index as usize] {
                bit_array[*index as usize] = false;
                break;
            }
        }
    }
    
    let mut bit_array_indices_vec_index = 0usize;
    while bit_array_indices_vec_index < bit_array_indices.len() {
        while !bit_array[bit_array_indices[bit_array_indices_vec_index] as usize] {
            bit_array_indices.remove(bit_array_indices_vec_index);
            if bit_array_indices_vec_index == bit_array_indices.len() {
                break;
            }
        }
        bit_array_indices_vec_index += 1;
    }
    return bit_array_indices;
}


#[derive(Clone)]
struct Node {
    children: Vec<Arc<RwLock<Node>>>,

    // u_bitarray = union of all bit arrays
    // c_bitarray = complete intersection of all bit arrays

    u_bitarray: Arc<Vec<bool>>,
    c_bitarray: Arc<Vec<bool>>,

    u_bitarray_indices: Arc<Vec<u32>>,
    c_bitarray_indices: Arc<Vec<u32>>,

    protein: Option<Arc<Protein>>,
}
impl Node {
    fn new_leaf(protein: Arc<Protein>, kmer_size: &u8) -> Node {
        let mut node = Node {
            children: vec![],

            u_bitarray: Arc::new(Vec::new()),
            c_bitarray: Arc::new(Vec::new()),

            u_bitarray_indices: Arc::new(Vec::new()),
            c_bitarray_indices: Arc::new(Vec::new()),

            protein: Some(protein.clone()),
        };
        node.add_bit_array_from_protein(protein, kmer_size);
        // eprintln!("Node u_bitarray: {:#?}; c_bitarray: {:#?}; u_bitarray_indices: {:#?}; c_bitarray_indices: {:#?}",
        //     node.u_bitarray,
        //     node.c_bitarray,
        //     node.u_bitarray_indices,
        //     node.c_bitarray_indices);
        node
    }

    fn add_bit_array_from_protein(&mut self, protein: Arc<Protein>, kmer_size: &u8) {
        if !self.children.is_empty() {
            panic!("add_bit_array_from_protein should only be called on new leaf");
        }
        if *kmer_size == 5u8 {
            self.u_bitarray = Arc::new(protein.get_five_hash_map());
            self.c_bitarray = Arc::new(protein.get_five_hash_map());

            self.u_bitarray_indices = Arc::new(protein.get_five_hash());
            self.c_bitarray_indices = Arc::new(protein.get_five_hash());
        }
        else if *kmer_size == 7u8 {
            self.u_bitarray = Arc::new(protein.get_seven_hash_map());
            self.c_bitarray = Arc::new(protein.get_seven_hash_map());

            self.u_bitarray_indices = Arc::new(protein.get_seven_hash());
            self.c_bitarray_indices = Arc::new(protein.get_seven_hash());
        }
        else {
            panic!("{}", format!("kmer-size should be 5 or 7, not {kmer_size}"));
        }
    }

    // fn new_node_with_children(children: Vec<Arc<RwLock<Node>>>) -> Node {
    //     let array_of_u_bitarray_indices: Vec<Arc<Vec<u32>>> = children[1..]
    //         .iter().map(|x| x.read().unwrap().u_bitarray_indices
    //         .clone()).collect();
    //     let array_of_c_bitarrays: Vec<Arc<Vec<bool>>> = children[1..]
    //         .iter().map(|x| x.read().unwrap().c_bitarray
    //         .clone()).collect();

    //     let child_read_locked = children[0].read().unwrap();
    //     let u_bitarray = child_read_locked
    //         .u_bitarray.deref().clone();
    //     let u_bitarray_indices= child_read_locked
    //         .u_bitarray_indices.deref().clone();
    //     let c_bitarray = child_read_locked
    //         .c_bitarray.deref().clone();
    //     let c_bitarray_indices = child_read_locked
    //         .c_bitarray_indices.deref().clone();
    //     drop(child_read_locked);

    //     let (u_bitarray, u_bitarray_indices) = unionize_bitarrays(
    //         u_bitarray, u_bitarray_indices, 
    //         array_of_u_bitarray_indices);
    //     let c_bitarray_indices = intersect_bitarrays(
    //         c_bitarray, c_bitarray_indices, 
    //         array_of_c_bitarrays);
        
    //     let mut c_bitarray = vec![false; u_bitarray.len()];
    //     for index in c_bitarray_indices.iter() {
    //         c_bitarray[*index as usize] = true;
    //     }
        
    //     Node{
    //         children: children, 

    //         u_bitarray: Arc::new(u_bitarray),
    //         c_bitarray: Arc::new(c_bitarray),

    //         u_bitarray_indices: Arc::new(u_bitarray_indices),
    //         c_bitarray_indices: Arc::new(c_bitarray_indices),

    //         protein: None}
    // }

    fn clone_and_clean(&mut self) -> Node {
        let mut children: Vec<Arc<RwLock<Node>>> = vec![];
        while self.children.len() > 0 {
            children.push(self.children.pop().unwrap());
        }

        let u_bitarray = self.u_bitarray.clone();
        let c_bitarray = self.c_bitarray.clone();
        let u_bitarray_indices = self.u_bitarray_indices.clone();
        let c_bitarray_indices = self.c_bitarray_indices.clone();
        let protein = self.protein.clone();

        self.u_bitarray = Arc::new(Vec::new());
        self.c_bitarray = Arc::new(Vec::new());
        self.u_bitarray_indices = Arc::new(Vec::new());
        self.c_bitarray_indices = Arc::new(Vec::new());
        self.protein = None;

        Node { 
            children, 
            u_bitarray,
            c_bitarray, 
            u_bitarray_indices, 
            c_bitarray_indices, 
            protein,
        }
    }

    fn balance(mut curr: RwLockWriteGuard<'_, Node>) {
        // let mut max_u_similarity = (0usize, 0usize, 0usize);
        // let mut min_u_similarity: Option<usize> = None;

        let mut max_c_similarity = (0usize, 0usize, 0usize);
        let mut min_c_similarity: Option<usize> = None;
        for i in 1..curr.children.len() {
            let child_i_readlock = curr.children[i].read().unwrap();

            // let child_i_u_bitarray = child_i_readlock.u_bitarray.deref();
            let child_i_c_bitarray = child_i_readlock.c_bitarray.deref();

            // let child_i_u_bitarray_indices = child_i_readlock.u_bitarray_indices.deref();
            let child_i_c_bitarray_indices = child_i_readlock.c_bitarray_indices.deref();

            for j in 0..i {
                let child_j_readlock = curr.children[j].read().unwrap();
                // let u_intersect = intersect_bitarrays(
                //     child_i_u_bitarray.clone(), child_i_u_bitarray_indices.clone(),
                //     vec![child_j_readlock.u_bitarray.clone()]);
                let c_intersect = intersect_bitarrays(
                    child_i_c_bitarray.clone(), child_i_c_bitarray_indices.clone(),
                    vec![child_j_readlock.c_bitarray.clone()]);
                // let u_similarity = u_intersect.len();
                let c_similarity = c_intersect.len();
                // if u_similarity > max_u_similarity.0 {
                //     max_u_similarity = (u_similarity, i, j);
                // }
                if c_similarity > max_c_similarity.0 {
                    max_c_similarity = (c_similarity, i, j);
                }
                // if min_u_similarity.is_none() || min_u_similarity.unwrap() > u_similarity {
                //     min_u_similarity = Some(u_similarity);
                // }
                if min_c_similarity.is_none() || min_c_similarity.unwrap() > c_similarity {
                    min_c_similarity = Some(c_similarity);
                }
            }
        }
        
        if max_c_similarity.0 > min_c_similarity.unwrap() {
            eprintln!("Merging");
            if max_c_similarity.1 <= max_c_similarity.2 {
                panic!("How")
            }
            let child_one = curr.children[max_c_similarity.1].clone();
            let child_two = curr.children[max_c_similarity.2].clone();
            let child_one_children_len = child_one.read().unwrap().children.len();
            let child_two_children_len = child_two.read().unwrap().children.len();
            if child_one_children_len < child_two_children_len {
                let combined_lock = child_one.write().unwrap();
                curr.children.remove(max_c_similarity.2);
                drop(curr);
                Node::add_child(combined_lock, child_two);
            }
            else {
                let combined_lock = child_two.write().unwrap();
                curr.children.remove(max_c_similarity.1);
                drop(curr);
                Node::add_child(combined_lock, child_one);
            }
        }
        // else {
        //     eprintln!("No merging");
        // }
        // else if max_u_similarity.0 > min_u_similarity.unwrap() && curr.children.len() > 100 {
        //     if max_u_similarity.1 <= max_u_similarity.2 {
        //         panic!("How")
        //     }
        //     let child_one = curr.children[max_u_similarity.1].clone();
        //     let child_two = curr.children[max_u_similarity.2].clone();
        //     let child_one_children_len = child_one.read().unwrap().children.len();
        //     let child_two_children_len = child_two.read().unwrap().children.len();
        //     if child_one_children_len < child_two_children_len {
        //         let combined_lock = child_one.write().unwrap();
        //         curr.children.remove(max_u_similarity.2);
        //         drop(curr);
        //         Node::add_child(combined_lock, child_two);
        //     }
        //     else {
        //         let combined_lock = child_two.write().unwrap();
        //         curr.children.remove(max_u_similarity.1);
        //         drop(curr);
        //         Node::add_child(combined_lock, child_one);
        //     }
        // }
    }

    fn add_child(mut curr: RwLockWriteGuard<'_, Node>, child: Arc<RwLock<Node>>) {
        // If curr is a leaf: 
        //      Create new leaf node with curr protein,
        //      Transform curr into an inner node, and 
        //      If child is leaf: have curr adopt new leaf and child
        //      Else: have curr adopt new leaf and child's children
        if curr.children.len() == 0 {
            let child_read_lock = child.read().unwrap();
            let cloned = Arc::new(RwLock::new(curr.clone_and_clean()));
            let cloned_read_lock = cloned.read().unwrap();

            // eprintln!("Clone u_bitarray: {:#?}; c_bitarray: {:#?}; u_bitarray_indices: {:#?}; c_bitarray_indices: {:#?}",
            //     cloned_read_lock.u_bitarray,
            //     cloned_read_lock.c_bitarray,
            //     cloned_read_lock.u_bitarray_indices,
            //     cloned_read_lock.c_bitarray_indices);

            let array_of_u_bitarray_indices = vec![child_read_lock
                .u_bitarray_indices.clone()];
            let array_of_c_bitarrays = vec![child_read_lock
                .c_bitarray.clone()];

            let u_bitarray = cloned_read_lock
                .u_bitarray.deref().clone();
            let u_bitarray_indices= cloned_read_lock
                .u_bitarray_indices.deref().clone();
            let c_bitarray = cloned_read_lock
                .c_bitarray.deref().clone();
            let c_bitarray_indices = cloned_read_lock
                .c_bitarray_indices.deref().clone();

            let (u_bitarray, u_bitarray_indices) = unionize_bitarrays(
                u_bitarray.clone(), u_bitarray_indices.clone(), 
                array_of_u_bitarray_indices.clone());
            let c_bitarray_indices = intersect_bitarrays(
                c_bitarray, c_bitarray_indices, 
                array_of_c_bitarrays);
            

            let mut c_bitarray = vec![false; u_bitarray.len()];
            for index in c_bitarray_indices.iter() {
                c_bitarray[*index as usize] = true;
            }

            curr.u_bitarray = Arc::new(u_bitarray);
            curr.c_bitarray = Arc::new(c_bitarray);
            curr.u_bitarray_indices = Arc::new(u_bitarray_indices);
            curr.c_bitarray_indices = Arc::new(c_bitarray_indices);
            curr.children.push(cloned.clone());
            if child_read_lock.children.len() == 0 {
                curr.children.push(child.clone());
            }
            else {
                let mut children_of_child = child_read_lock.children.clone();
                curr.children.append(&mut children_of_child);
                drop(child_read_lock);
                drop(cloned_read_lock);
            }
        }
        // If curr is not a leaf
        else {
            let child_read_lock = child.read().unwrap();
            let curr_u_bitarray = curr.u_bitarray.deref().clone();
            let curr_u_bitarray_indices = curr.u_bitarray_indices.deref().clone();
            let intersect = intersect_bitarrays(
                    curr_u_bitarray, curr_u_bitarray_indices,
                    vec![child_read_lock.u_bitarray.clone()]);

            // eprintln!("Curr u_bitarray: {:#?}; u_bitarray_indices: {:#?}. \n Child: u_bitarray_indices: {:#?}. \n Intersect: {:#?}",
            //     curr.u_bitarray,
            //     curr.u_bitarray_indices,
            //     child_read_lock.u_bitarray_indices,
            //     intersect);

            // If curr's proteins and child have no kmer in common:
            //      Have curr adopt child
            //      Balance curr if needed
            // if intersect.is_empty() {
            let array_of_u_bitarray_indices = vec![child_read_lock
                .u_bitarray_indices.clone()];
            let array_of_c_bitarrays = vec![child_read_lock
                .c_bitarray.clone()];

            let u_bitarray = curr
                .u_bitarray.deref().clone();
            let u_bitarray_indices= curr
                .u_bitarray_indices.deref().clone();
            let c_bitarray = curr
                .c_bitarray.deref().clone();
            let c_bitarray_indices = curr
                .c_bitarray_indices.deref().clone();

            let (u_bitarray, u_bitarray_indices) = unionize_bitarrays(
                u_bitarray.clone(), u_bitarray_indices.clone(), 
                array_of_u_bitarray_indices.clone());
            let c_bitarray_indices = intersect_bitarrays(
                c_bitarray, c_bitarray_indices, 
                array_of_c_bitarrays);
            

            let mut c_bitarray = vec![false; u_bitarray.len()];
            for index in c_bitarray_indices.iter() {
                c_bitarray[*index as usize] = true;
            }

            curr.u_bitarray = Arc::new(u_bitarray);
            curr.c_bitarray = Arc::new(c_bitarray);
            curr.u_bitarray_indices = Arc::new(u_bitarray_indices);
            curr.c_bitarray_indices = Arc::new(c_bitarray_indices);
            curr.children.push(child.clone());
            drop(child_read_lock);
            //println!("Newly added protein has no kmer in common with current node");
            if !intersect.is_empty() {
                Node::balance(curr);
            }
            else {
                eprintln!("No kmers in common");
            }
            }
            // else {
            //     let mut common_kmer: Vec<Arc<RwLock<Node>>> = vec![child.clone()];
            //     let mut no_common_kmer: Vec<Arc<RwLock<Node>>> = vec![];
            //     let mut children_read_locks: Vec<RwLockReadGuard<'_, Node>> = vec![];

            //     // Divide children into two lists and keep track of read locks
            //     for curr_child in &curr.children {
            //         children_read_locks.push(curr_child.read().unwrap());
            //         eprintln!("218 - Read lock successful");
            //         let intersect = intersect_bitarrays(
            //             children_read_locks.last().unwrap().u_bitarray.clone(),
            //             child_read_lock.u_bitarray_indices.clone());
            //         if !intersect.is_empty() {
            //             common_kmer.push(curr_child.clone());
            //         }
            //         else {
            //             no_common_kmer.push(curr_child.clone());
            //         }
            //     }
            //     drop(children_read_locks);

            //     // Find bit array of kmers for both lists
            //     let common_kmer = Node::new_node_with_children(common_kmer);
            //     let no_common_kmer = {
            //         if no_common_kmer.len() > 0 {
            //             Some(Node::new_node_with_children(no_common_kmer))
            //         }
            //         else {
            //             None
            //         }
            //     };
            //     let intersect = {
            //         if no_common_kmer.is_some() {
            //             intersect_bitarrays(common_kmer.u_bitarray.clone(),
            //                 no_common_kmer.as_ref().unwrap().u_bitarray_indices.clone())
            //         }
            //         else {
            //             vec![1u32]
            //         }
            //     };
                
            //     // If aformentioned list of nodes also have kmer in common with remaining nodes:
            //     //      Have curr adopt child
            //     //      Balance curr if needed
            //     if !intersect.is_empty() {
            //         drop(common_kmer);
            //         drop(no_common_kmer);
            //         let array_of_bit_array_indices = vec![child_read_lock
            //             .u_bitarray_indices.clone()];
            //         let u_bitarray = curr.u_bitarray.deref().clone();
            //         let u_bitarray_indices= curr.u_bitarray_indices.deref().clone();
            //         let (u_bitarray, u_bitarray_indices) = unionize_bitarrays(
            //             u_bitarray, u_bitarray_indices, array_of_bit_array_indices);
            //         curr.u_bitarray = Arc::new(u_bitarray);
            //         curr.u_bitarray_indices = Arc::new(u_bitarray_indices);
            //         curr.children.push(child.clone());
            //         drop(child_read_lock);
            //         eprintln!("263 - Drop lock successful");
            //         //println!("Can't separate based on kmer");
            //         if curr.children.len() >= 5 {
            //             Node::balance(curr);
            //         }
            //     }

            //     // If only one node has kmers in common with child:
            //     //      Update curr's bit array
            //     //      Call add_child on that node
            //     else if common_kmer.children.len() == 2 {
            //         let new_parent = common_kmer.children[1].clone();
            //         drop(common_kmer);
            //         drop(no_common_kmer);
            //         let array_of_bit_array_indices = vec![child_read_lock
            //             .u_bitarray_indices.clone()];
            //         let u_bitarray = curr.u_bitarray.deref().clone();
            //         let u_bitarray_indices= curr.u_bitarray_indices.deref().clone();
            //         let (u_bitarray, u_bitarray_indices) = unionize_bitarrays(
            //             u_bitarray, u_bitarray_indices, array_of_bit_array_indices);
            //         curr.u_bitarray = Arc::new(u_bitarray);
            //         curr.u_bitarray_indices = Arc::new(u_bitarray_indices);
                    
            //         let new_parent = new_parent.write().unwrap();
            //         drop(curr);
            //         drop(child_read_lock);
            //         eprintln!("280 - Lock successful");
            //         Node::add_child(new_parent, child.clone());
            //     }

            //     // If aformentioned list of nodes have no kmer in common with remaining nodes:
            //     //      Merge nodes in common into single parent node
            //     //      Merge nodes in no_commin into another parent node
            //     //      Repeat what was done in scenario where curr had no children
            //     else {
            //         let array_of_bit_array_indices = vec![child_read_lock
            //             .u_bitarray_indices.clone()];
            //         let u_bitarray = curr.u_bitarray.deref().clone();
            //         let u_bitarray_indices= curr.u_bitarray_indices.deref().clone();
            //         let (u_bitarray, u_bitarray_indices) = unionize_bitarrays(
            //             u_bitarray, u_bitarray_indices, array_of_bit_array_indices);
            //         curr.u_bitarray = Arc::new(u_bitarray);
            //         curr.u_bitarray_indices = Arc::new(u_bitarray_indices);
            //         drop(child_read_lock);
            //         let common_node = Arc::new(RwLock::new(common_kmer));
            //         let no_common_node = Arc::new(RwLock::new(no_common_kmer.unwrap()));
            //         curr.children = vec![common_node.clone(), no_common_node];
            //         //println!("Merged based on kmer");
            //         let common_node_read_lock = common_node.read().unwrap();
            //         if common_node_read_lock.children.len() > 5 {
            //             drop(common_node_read_lock);
            //             let common_node_write_lock = common_node.write().unwrap();
            //             drop(curr);
            //             eprintln!("318 - Lock successful");
            //             Node::balance(common_node_write_lock);
            //         }
            //     }
        //    }
    //    }
    }

}
impl fmt::Debug for Node {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.protein.is_some() {
            f.debug_struct("Leaf").field("protein", &self.protein).finish()
        }
        else {
            let read_locked_nodes: Vec<Node> = self.children.iter().map(
                |x| x.read().unwrap().clone()).collect();
            f.debug_list().entries(read_locked_nodes.iter()).finish()
        }
    }
}

#[derive(Debug)]
pub struct Tree {
    root: Arc<RwLock<Node>>,
    kmer_size: u8,
}
impl Tree {
    pub fn new(kmer_size: u8, protein: Arc<Protein>) -> Tree {
        let root = Arc::new(RwLock::new(Node::new_leaf(protein, &kmer_size)));
        Tree{
            root,
            kmer_size
        }
    }
    pub fn add_protein(&self, protein: Arc<Protein>) {
        let node = Arc::new(RwLock::new(Node::new_leaf(protein, &self.kmer_size)));
        let root_clone = self.root.clone();
        let root_clone = root_clone.write().unwrap();
        Node::add_child(root_clone, node);
    }
}