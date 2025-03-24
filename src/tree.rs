use std::sync::{Arc, RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::fmt;

use crate::protein::Protein;

fn binary_or_bit_array(bit_arrays: Vec<Arc<Vec<bool>>>) -> Arc<Vec<bool>> {
    let length = bit_arrays.len();
    if length == 1 {
        bit_arrays[0].clone()
    }
    else if length == 2 {
        Arc::new(bit_arrays[0].iter().zip(bit_arrays[1].iter())
            .map(|(x,y)| x | y).collect())
    }
    else {
        Arc::new(binary_or_bit_array(bit_arrays[..length/2].to_vec()).iter()
            .zip(binary_or_bit_array(bit_arrays[length/2..].to_vec()).iter())
            .map(|(x,y)| x | y).collect())
    }
}
fn binary_and_bit_array(first_bit_array: Arc<Vec<bool>>, second_bit_array: Arc<Vec<bool>>) -> Arc<Vec<bool>> {
    Arc::new(first_bit_array.iter().zip(second_bit_array.iter())
        .map(|(x,y)| x & y).collect())
}

fn at_least_one_true(bit_array: Arc<Vec<bool>>) -> bool {
    for bit in bit_array.iter() {
        if *bit {
            return true
        }
    }
    false
}
fn calculate_similarity(bit_array: Arc<Vec<bool>>) -> u32 {
    let mut number_true = 0u32;
    for bit in bit_array.iter() {
        if *bit {
            number_true += 1;
        }
    }
    number_true
}

#[derive(Clone)]
struct Node {
    children: Vec<Arc<RwLock<Node>>>,
    bit_array: Arc<Vec<bool>>,
    protein: Option<Arc<Protein>>,
}
impl Node {
    fn new_leaf(protein: Arc<Protein>, kmer_size: &u8) -> Node {
        let mut node = Node {
            children: vec![],
            bit_array: Arc::new(Vec::new()),
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
            self.bit_array = Arc::new(protein.get_five_hash());
        }
        else if *kmer_size == 7u8 {
            self.bit_array = Arc::new(protein.get_seven_hash());
        }
        else {
            panic!("{}", format!("kmer-size should be 5 or 7, not {kmer_size}"));
        }
    }

    fn new_node_with_children(children: Vec<Arc<RwLock<Node>>>) -> Node {
        let array_of_bit_arrays: Vec<Arc<Vec<bool>>> = children.iter()
            .map(|x| x.read().unwrap().bit_array.clone()).collect();
        let bit_array = binary_or_bit_array(array_of_bit_arrays);
        Node{children, bit_array, protein: None}
    }

    fn clone_and_clean(&mut self) -> Node {
        let mut children_copy: Vec<Arc<RwLock<Node>>> = vec![];
        while self.children.len() > 0 {
            children_copy.push(self.children.pop().unwrap());
        }
        let bit_array_copy = self.bit_array.clone();
        let protein_copy = self.protein.clone();
        self.bit_array = Arc::new(Vec::new());
        self.protein = None;
        Node {children: children_copy, bit_array: bit_array_copy, protein: protein_copy}
        
    }
    
    fn balance(mut curr: RwLockWriteGuard<'_, Node>) {
        let mut max_similarity = (0u32, 0usize, 0usize);
        let mut min_similarity: Option<u32> = None;
        for i in 1..curr.children.len() {
            for j in 0..i {
                let similarity_bit_array = binary_and_bit_array(
                    curr.children[i].read().unwrap().bit_array.clone(),
                    curr.children[j].read().unwrap().bit_array.clone());
                let similarity = calculate_similarity(similarity_bit_array);
                if similarity > max_similarity.0 {
                    max_similarity = (similarity, i, j);
                }
                if min_similarity.is_none() || min_similarity.unwrap() > similarity {
                    min_similarity = Some(similarity);
                }
            }
        }
        //eprintln!("109 - Read locks successful");
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
                //eprintln!("120 - Write lock successful");
                curr.children.remove(max_similarity.2);
                if child_one_children_len > 0 {
                    combined_lock.bit_array = binary_or_bit_array(vec![
                        combined_lock.bit_array.clone(), 
                        child_two.read().unwrap().bit_array.clone()]);
                    //eprintln!("126 - Read lock successful");
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
                //eprintln!("136 - Write lock successful");
                curr.children.remove(max_similarity.1);
                if child_two_children_len > 0 {
                    combined_lock.bit_array = binary_or_bit_array(vec![
                        combined_lock.bit_array.clone(), 
                        child_one.read().unwrap().bit_array.clone()]);
                    //eprintln!("142 - Read lock successful");
                    combined_lock.children.push(child_one);
                    //self.balance();
                }
                else {
                    drop(curr);
                    Node::add_child(combined_lock, child_one);
                }
            }
        }
        else {
            if curr.children.len() >= 10 {
                let first = curr.children.len() / 5;
                let second = 2 * curr.children.len() / 5;
                let third = 3 * curr.children.len() / 5;
                let fourth = 4 * curr.children.len() / 5;
                let child_one = Arc::new(RwLock::new(
                    Node::new_node_with_children(curr.children[..first].to_vec())));
                let child_two = Arc::new(RwLock::new(
                    Node::new_node_with_children(curr.children[first..second].to_vec())));
                let child_three = Arc::new(RwLock::new(
                    Node::new_node_with_children(curr.children[second..third].to_vec())));
                let child_four = Arc::new(RwLock::new(
                    Node::new_node_with_children(curr.children[third..fourth].to_vec())));
                let child_five = Arc::new(RwLock::new(
                    Node::new_node_with_children(curr.children[fourth..].to_vec())));
                curr.children = vec![child_one, child_two, child_three, child_four, child_five];
            }
        }
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
            //eprintln!("183 - Read locks successful");
            curr.bit_array = binary_or_bit_array(vec![
                    child_read_lock.bit_array.clone(), 
                    cloned_read_locked.bit_array.clone()]);
            curr.children.push(cloned.clone());
            if child_read_lock.children.len() == 0 {
                curr.children.push(child.clone());
            }
            else {
                let mut children_of_child = child_read_lock.children.clone();
                curr.children.append(&mut children_of_child);
                drop(child_read_lock);
                drop(cloned_read_locked);
                //eprintln!("196 - Drop lock successful");
                if curr.children.len() >= 5 {
                    Node::balance(curr);
                }
            }
        }
        // If curr is not a leaf
        else {
            let child_read_lock = child.read().unwrap();
            // let comparison_bit_array = binary_and_bit_array(
            //     curr.bit_array.clone(), child_read_lock.bit_array.clone());

            // // If curr's proteins and child have no kmer in common:
            // //      Have curr adopt child
            // //      Balance curr if needed
            // if !at_least_one_true(comparison_bit_array.clone()) {
                curr.bit_array = binary_or_bit_array(vec![
                    child_read_lock.bit_array.clone(), curr.bit_array.clone()]);
                curr.children.push(child.clone());
                drop(child_read_lock);
                //eprintln!("340 - Drop lock successful");
                //println!("Newly added protein has no kmer in common with current node");
                if curr.children.len() >= 5 {
                    Node::balance(curr);
                }
            // }
       }
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

    // fn add_child(mut curr: RwLockWriteGuard<'_, Node>, child: Arc<RwLock<Node>>) {
    //     // If curr is a leaf: 
    //     //      Create new leaf node with curr protein,
    //     //      Transform curr into an inner node, and 
    //     //      If child is leaf: have curr adopt new leaf and child
    //     //      Else: have curr adopt new leaf and child's children
    //     if curr.children.len() == 0 {
    //         let child_read_lock = child.read().unwrap();
    //         let cloned = Arc::new(RwLock::new(curr.clone_and_clean()));
    //         let cloned_read_locked = cloned.read().unwrap();
    //         //eprintln!("183 - Read locks successful");
    //         curr.bit_array = binary_or_bit_array(vec![
    //                 child_read_lock.bit_array.clone(), 
    //                 cloned_read_locked.bit_array.clone()]);
    //         curr.children.push(cloned.clone());
    //         if child_read_lock.children.len() == 0 {
    //             curr.children.push(child.clone());
    //         }
    //         else {
    //             let mut children_of_child = child_read_lock.children.clone();
    //             curr.children.append(&mut children_of_child);
    //             drop(child_read_lock);
    //             drop(cloned_read_locked);
    //             //eprintln!("196 - Drop lock successful");
    //             Node::balance(curr);
    //         }
    //     }
    //     // If curr is not a leaf
    //     else {
    //         let child_read_lock = child.read().unwrap();
    //         let comparison_bit_array = binary_and_bit_array(
    //             curr.bit_array.clone(), child_read_lock.bit_array.clone());

    //         // If curr's proteins have at least one kmer in common with child:
    //         //      Find among curr's children a list of nodes with at least one kmer in common with child
    //         if at_least_one_true(comparison_bit_array.clone()) {
    //             // let mut common_kmer: Vec<Arc<RwLock<Node>>> = vec![];
    //             // let mut no_common_kmer: Vec<Arc<RwLock<Node>>> = vec![];
    //             // let mut children_read_locks: Vec<RwLockReadGuard<'_, Node>> = vec![];

    //             // // Divide children into two lists and keep track of read locks
    //             // for curr_child in &curr.children {
    //             //     children_read_locks.push(curr_child.read().unwrap());
    //             //     //eprintln!("218 - Read lock successful");
    //             //     let comparison_bit_array = binary_and_bit_array(
    //             //         children_read_locks.last().unwrap().bit_array.clone(),
    //             //         child_read_lock.bit_array.clone());
    //             //     if at_least_one_true(comparison_bit_array) {
    //             //         common_kmer.push(curr_child.clone());
    //             //     }
    //             //     else {
    //             //         no_common_kmer.push(curr_child.clone());
    //             //     }
    //             // }
    //             // drop(children_read_locks);

    //             // // Find bit array of kmers for both lists
    //             // let common_kmer_bit_arrays: Vec<Arc<Vec<bool>>> = common_kmer.iter()
    //             //     .map(|x| x.read().unwrap().bit_array.clone()).collect();
    //             // let common_kmer_bit_array = binary_or_bit_array(common_kmer_bit_arrays);
    //             // let no_common_kmer_bit_array = {
    //             //     if no_common_kmer.len() > 0 {
    //             //         let no_common_kmer_bit_arrays: Vec<Arc<Vec<bool>>> = no_common_kmer.iter()
    //             //             .map(|x| x.read().unwrap().bit_array.clone()).collect();
    //             //         binary_or_bit_array(no_common_kmer_bit_arrays)
    //             //     }
    //             //     else {
    //             //         common_kmer_bit_array.clone()
    //             //     }
    //             // };
    //             // let comparison_bit_array = {
    //             //     if no_common_kmer.len() > 0 {
    //             //         binary_and_bit_array(common_kmer_bit_array.clone(),
    //             //             no_common_kmer_bit_array.clone())
    //             //     }
    //             //     else {
    //             //         Arc::new(vec![true])
    //             //     }
    //             // };
                
    //             // // If aformentioned list of nodes also have kmer in common with remaining nodes:
    //             // //      Have curr adopt child
    //             // //      Balance curr if needed
    //             // if at_least_one_true(comparison_bit_array) {
    //             //     curr.bit_array = binary_or_bit_array(vec![
    //             //         child_read_lock.bit_array.clone(), curr.bit_array.clone()]);
    //             //     curr.children.push(child.clone());
    //             //     drop(child_read_lock);
    //             //     //eprintln!("263 - Drop lock successful");
    //             //     //println!("Can't separate based on kmer");
    //             //     if curr.children.len() >= 5 {
    //             //         Node::balance(curr);
    //             //     }
    //             // }

    //             // // If only one node has kmers in common with child:
    //             // //      Update curr's bit array
    //             // //      Call add_child on that node
    //             // else if common_kmer.len() == 1 {
    //             //     let new_parent = common_kmer.first().unwrap().clone();
    //             //     curr.bit_array = binary_or_bit_array(vec![
    //             //         child_read_lock.bit_array.clone(), curr.bit_array.clone()]);
    //             //     drop(curr);
    //             //     let new_parent = new_parent.write().unwrap();
    //             //     drop(child_read_lock);
    //             //     //eprintln!("280 - Lock successful");
    //             //     Node::add_child(new_parent, child.clone());
    //             // }

    //             // // If aformentioned list of nodes have no kmer in common with remaining nodes:
    //             // //      Merge nodes in common into single parent node
    //             // //      Merge nodes in no_commin into another parent node
    //             // //      Repeat what was done in scenario where curr had no children
    //             // else {
    //             //     common_kmer.push(child.clone());
    //             //     let common_kmer_bit_array = binary_or_bit_array(vec![
    //             //         child_read_lock.bit_array.clone(), common_kmer_bit_array.clone()]);
    //             //     let common_node = Arc::new(RwLock::new(
    //             //         Node {children: common_kmer, bit_array: common_kmer_bit_array, protein: None}
    //             //     ));
                        
    //             //     let no_common_node = {
    //             //         if no_common_kmer.len() == 1 {
    //             //             no_common_kmer.first().unwrap().clone()
    //             //         }
    //             //         else {
    //             //             Arc::new(RwLock::new(
    //             //                 Node {
    //             //                     children: no_common_kmer, 
    //             //                     bit_array: no_common_kmer_bit_array,
    //             //                     protein: None}))
    //             //         }
    //             //     };
    //             //     curr.bit_array = binary_or_bit_array(vec![
    //             //         child_read_lock.bit_array.clone(), curr.bit_array.clone()]);
    //             //     drop(child_read_lock);
    //             //     curr.children = vec![common_node.clone(), no_common_node.clone()];
    //             //     //println!("Merged based on kmer");
    //             //     let common_node_read_lock = common_node.read().unwrap();
    //             //     if common_node_read_lock.children.len() > 5 {
    //             //         drop(common_node_read_lock);
    //             //         let mut common_node_write_lock = common_node.write().unwrap();
    //             //         drop(curr);
    //             //         //eprintln!("318 - Lock successful");
    //             //         Node::balance(common_node_write_lock);
    //             //     }
    //                 curr.bit_array = binary_or_bit_array(vec![
    //                     child_read_lock.bit_array.clone(), curr.bit_array.clone()]);
    //                 curr.children.push(child.clone());
    //                 drop(child_read_lock);
    //                 // //println!("Can't separate based on kmer");
    //                 if curr.children.len() > 2 {
    //                     Node::balance(curr);
    //                 }
    //             //}
    //         }

    //         // If curr's proteins and child have no kmer in common:
    //         //      Have curr adopt child
    //         //      Balance curr if needed
    //         else {
    //             curr.bit_array = binary_or_bit_array(vec![
    //                 child_read_lock.bit_array.clone(), curr.bit_array.clone()]);
    //             curr.children.push(child.clone());
    //             drop(child_read_lock);
    //             //eprintln!("340 - Drop lock successful");
    //             //println!("Newly added protein has no kmer in common with current node");
    //             if curr.children.len() >= 10 {
    //                 Node::balance(curr);
    //             }
    //         }
    //     }
    //     // let array_of_bit_arrays
    // }