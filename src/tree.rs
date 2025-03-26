use std::ops::Deref;
use std::sync::{Arc, RwLock, RwLockReadGuard, RwLockWriteGuard};
//use std::sync::{Arc, RwLock, RwLockWriteGuard};
use std::fmt;

use crate::protein::Protein;

fn unionize_bitarrays(
        mut bit_array: Vec<bool>, mut bit_index_array: Vec<u32>, 
        array_of_bit_index_array: Vec<Arc<Vec<u32>>>) -> (Vec<bool>, Vec<u32>) {
    for other_bit_index_array in array_of_bit_index_array {
        for other_bit_index in other_bit_index_array.deref() {
            if !bit_array[*other_bit_index as usize] {
                bit_array[*other_bit_index as usize] = true;
                bit_index_array.push(*other_bit_index);
            }
        }
    }
    return (bit_array, bit_index_array);
}
fn intersect_bitarrays(
    first_bit_array: Arc<Vec<bool>>, second_bit_index_array: Arc<Vec<u32>>) -> Vec<u32> {
    let mut intersect = vec![];
    for other_bit_index in second_bit_index_array.deref() {
        if first_bit_array[*other_bit_index as usize] {
            intersect.push(*other_bit_index);
        }
    }
    return intersect;
}


#[derive(Clone)]
struct Node {
    children: Vec<Arc<RwLock<Node>>>,
    bit_array: Arc<Vec<bool>>,
    bit_index_array: Arc<Vec<u32>>,
    protein: Option<Arc<Protein>>,
}
impl Node {
    fn new_leaf(protein: Arc<Protein>, kmer_size: &u8) -> Node {
        let mut node = Node {
            children: vec![],
            bit_array: Arc::new(Vec::new()),
            bit_index_array: Arc::new(Vec::new()),
            protein: Some(protein.clone()),
        };
        node.add_bit_array_from_protein(protein, kmer_size);
        node
    }

    fn add_bit_array_from_protein(&mut self, protein: Arc<Protein>, kmer_size: &u8) {
        if !self.children.is_empty() {
            panic!("add_bit_array_from_protein should only be called on new leaf");
        }
        if *kmer_size == 5u8 {
            self.bit_array = Arc::new(protein.get_five_hash_map());
            self.bit_index_array = Arc::new(protein.get_five_hash());
        }
        else if *kmer_size == 7u8 {
            self.bit_array = Arc::new(protein.get_seven_hash_map());
            self.bit_index_array = Arc::new(protein.get_seven_hash());
        }
        else {
            panic!("{}", format!("kmer-size should be 5 or 7, not {kmer_size}"));
        }
    }

    fn new_node_with_children(children: Vec<Arc<RwLock<Node>>>) -> Node {
        let array_of_bit_index_array: Vec<Arc<Vec<u32>>> = children[1..]
            .iter().map(|x| x.read().unwrap().bit_index_array.clone()).collect();
        let child_read_locked = children[0].read().unwrap();
        let bit_array = child_read_locked.bit_array.deref().clone();
        let bit_index_array= child_read_locked.bit_index_array.deref().clone();
        drop(child_read_locked);
        let (bit_array, bit_index_array) = unionize_bitarrays(
            bit_array, bit_index_array.clone(), array_of_bit_index_array.clone());
        Node{
            children: children, 
            bit_array: Arc::new(bit_array), 
            bit_index_array: Arc::new(bit_index_array), 
            protein: None}
    }

    fn clone_and_clean(&mut self) -> Node {
        let mut children_copy: Vec<Arc<RwLock<Node>>> = vec![];
        while self.children.len() > 0 {
            children_copy.push(self.children.pop().unwrap());
        }
        let bit_array_copy = self.bit_array.clone();
        let bit_index_array_clone = self.bit_index_array.clone();
        let protein_copy = self.protein.clone();
        self.bit_array = Arc::new(Vec::new());
        self.bit_index_array = Arc::new(Vec::new());
        self.protein = None;
        Node {
            children: children_copy, 
            bit_array: bit_array_copy, 
            bit_index_array: bit_index_array_clone,
            protein: protein_copy}
        
    }
    
    fn balance(mut curr: RwLockWriteGuard<'_, Node>) {
        let mut max_similarity = (0usize, 0usize, 0usize);
        let mut min_similarity: Option<usize> = None;
        for i in 1..curr.children.len() {
            for j in 0..i {
                let intersect = intersect_bitarrays(
                    curr.children[i].read().unwrap().bit_array.clone(),
                    curr.children[j].read().unwrap().bit_index_array.clone());
                let similarity = intersect.len();
                if similarity > max_similarity.0 {
                    max_similarity = (similarity, i, j);
                }
                if min_similarity.is_none() || min_similarity.unwrap() > similarity {
                    min_similarity = Some(similarity);
                }
            }
        }
        eprintln!("136 - Read locks successful");
        if max_similarity.0 > min_similarity.unwrap() {
            if max_similarity.1 <= max_similarity.2 {
                panic!("How")
            }
            let child_one = curr.children[max_similarity.1].clone();
            let child_two = curr.children[max_similarity.2].clone();
            let child_one_children_len = child_one.read().unwrap().children.len();
            let child_two_children_len = child_two.read().unwrap().children.len();
            if child_one_children_len < child_two_children_len {
                let mut combined_lock = child_one.write().unwrap();
                eprintln!("147 - Write lock successful");
                curr.children.remove(max_similarity.2);
                if child_one_children_len > 0 {
                    let array_of_bit_index_array = vec![child_two.read()
                        .unwrap().bit_index_array.clone()];
                    eprintln!("152 - Read lock successful");
                    let bit_array = combined_lock.bit_array.deref().clone();
                    let bit_index_array= combined_lock.bit_index_array.deref().clone();
                    let (bit_array, bit_index_array) = unionize_bitarrays(
                        bit_array, bit_index_array, array_of_bit_index_array);
                    combined_lock.bit_array = Arc::new(bit_array);
                    combined_lock.bit_index_array = Arc::new(bit_index_array);
                    combined_lock.children.push(child_two);
                    //self.balance();
                }
                else {
                    drop(curr);
                    Node::add_child(combined_lock, child_two);
                }
            }
            else {
                let mut combined_lock = child_two.write().unwrap();
                eprintln!("169 - Write lock successful");
                curr.children.remove(max_similarity.1);
                if child_two_children_len > 0 {
                    let array_of_bit_index_array = vec![child_one.read()
                        .unwrap().bit_index_array.clone()];
                    eprintln!("174 - Read lock successful");
                    let bit_array = combined_lock.bit_array.deref().clone();
                    let bit_index_array= combined_lock.bit_index_array.deref().clone();
                    let (bit_array, bit_index_array) = unionize_bitarrays(
                        bit_array, bit_index_array, array_of_bit_index_array);
                    combined_lock.bit_array = Arc::new(bit_array);
                    combined_lock.bit_index_array = Arc::new(bit_index_array);
                    combined_lock.children.push(child_one);
                    //self.balance();
                }
                else {
                    drop(curr);
                    Node::add_child(combined_lock, child_one);
                }
            }
        }
        // else {
        //     if curr.children.len() >= 10 {
        //         let first = curr.children.len() / 5;
        //         let second = 2 * curr.children.len() / 5;
        //         let third = 3 * curr.children.len() / 5;
        //         let fourth = 4 * curr.children.len() / 5;
        //         let child_one = Arc::new(RwLock::new(
        //             Node::new_node_with_children(curr.children[..first].to_vec())));
        //         let child_two = Arc::new(RwLock::new(
        //             Node::new_node_with_children(curr.children[first..second].to_vec())));
        //         let child_three = Arc::new(RwLock::new(
        //             Node::new_node_with_children(curr.children[second..third].to_vec())));
        //         let child_four = Arc::new(RwLock::new(
        //             Node::new_node_with_children(curr.children[third..fourth].to_vec())));
        //         let child_five = Arc::new(RwLock::new(
        //             Node::new_node_with_children(curr.children[fourth..].to_vec())));
        //         curr.children = vec![child_one, child_two, child_three, child_four, child_five];
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
            let cloned_read_locked = cloned.read().unwrap();
            eprintln!("221 - Read locks successful");
            let array_of_bit_index_array = vec![child_read_lock
                .bit_index_array.clone()];
            let bit_array = cloned_read_locked.bit_array.deref().clone();
            let bit_index_array= cloned_read_locked.bit_index_array.deref().clone();
            let (bit_array, bit_index_array) = unionize_bitarrays(
                bit_array, bit_index_array, array_of_bit_index_array);
            curr.bit_array = Arc::new(bit_array);
            curr.bit_index_array = Arc::new(bit_index_array);
            curr.children.push(cloned.clone());
            if child_read_lock.children.len() == 0 {
                curr.children.push(child.clone());
            }
            else {
                let mut children_of_child = child_read_lock.children.clone();
                curr.children.append(&mut children_of_child);
                drop(child_read_lock);
                drop(cloned_read_locked);
                eprintln!("239 - Drop lock successful");
                // if curr.children.len() >= 5 {
                //     Node::balance(curr);
                // }
            }
        }
        // If curr is not a leaf
        else {
            let child_read_lock = child.read().unwrap();
            let intersect = intersect_bitarrays(
                    curr.bit_array.clone(),
                    child_read_lock.bit_index_array.clone());

            // If curr's proteins and child have no kmer in common:
            //      Have curr adopt child
            //      Balance curr if needed
            // if intersect.is_empty() {
                let array_of_bit_index_array = vec![child_read_lock
                    .bit_index_array.clone()];
                let bit_array = curr.bit_array.deref().clone();
                let bit_index_array= curr.bit_index_array.deref().clone();
                let (bit_array, bit_index_array) = unionize_bitarrays(
                    bit_array, bit_index_array, array_of_bit_index_array);
                curr.bit_array = Arc::new(bit_array);
                curr.bit_index_array = Arc::new(bit_index_array);
                curr.children.push(child.clone());
                drop(child_read_lock);
                eprintln!("265 - Drop lock successful");
                //println!("Newly added protein has no kmer in common with current node");
                if curr.children.len() >= 5 {
                    Node::balance(curr);
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
            //             children_read_locks.last().unwrap().bit_array.clone(),
            //             child_read_lock.bit_index_array.clone());
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
            //             intersect_bitarrays(common_kmer.bit_array.clone(),
            //                 no_common_kmer.as_ref().unwrap().bit_index_array.clone())
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
            //         let array_of_bit_index_array = vec![child_read_lock
            //             .bit_index_array.clone()];
            //         let bit_array = curr.bit_array.deref().clone();
            //         let bit_index_array= curr.bit_index_array.deref().clone();
            //         let (bit_array, bit_index_array) = unionize_bitarrays(
            //             bit_array, bit_index_array, array_of_bit_index_array);
            //         curr.bit_array = Arc::new(bit_array);
            //         curr.bit_index_array = Arc::new(bit_index_array);
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
            //         let array_of_bit_index_array = vec![child_read_lock
            //             .bit_index_array.clone()];
            //         let bit_array = curr.bit_array.deref().clone();
            //         let bit_index_array= curr.bit_index_array.deref().clone();
            //         let (bit_array, bit_index_array) = unionize_bitarrays(
            //             bit_array, bit_index_array, array_of_bit_index_array);
            //         curr.bit_array = Arc::new(bit_array);
            //         curr.bit_index_array = Arc::new(bit_index_array);
                    
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
            //         let array_of_bit_index_array = vec![child_read_lock
            //             .bit_index_array.clone()];
            //         let bit_array = curr.bit_array.deref().clone();
            //         let bit_index_array= curr.bit_index_array.deref().clone();
            //         let (bit_array, bit_index_array) = unionize_bitarrays(
            //             bit_array, bit_index_array, array_of_bit_index_array);
            //         curr.bit_array = Arc::new(bit_array);
            //         curr.bit_index_array = Arc::new(bit_index_array);
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