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

#[derive(Clone)]
struct Node {
    children: Vec<Arc<RwLock<Node>>>,
    bit_array: Arc<Vec<bool>>,
    protein: Option<Arc<RwLock<Protein>>>,
}
impl Node {
    fn new_leaf(protein: Arc<RwLock<Protein>>, kmer_size: &u8) -> Node {
        let mut node = Node {
            children: vec![],
            bit_array: Arc::new(Vec::new()),
            protein: Some(protein.clone()),
        };
        node.add_bit_array_from_protein(protein, kmer_size);
        node
    }

    fn add_bit_array_from_protein(&mut self, protein: Arc<RwLock<Protein>>, kmer_size: &u8) {
        if !self.children.is_empty() {
            panic!("add_bit_array_from_protein should only be called on new leaf");
        }
        if *kmer_size == 5u8 {
            self.bit_array = Arc::new(protein.read().unwrap().get_five_hash());
        }
        else if *kmer_size == 7u8 {
            self.bit_array = Arc::new(protein.read().unwrap().get_seven_hash());
        }
        else {
            panic!("{}", format!("kmer-size should be 5 or 7, not {kmer_size}"));
        }
    }

    // fn new_node_with_children(children: Vec<Arc<RwLock<Node>>>) -> Node {
    //     let array_of_bit_arrays: Vec<Arc<Vec<bool>>> = children.iter()
    //         .map(|x| x.read().unwrap().bit_array.clone()).collect();
    //     let bit_array = binary_or_bit_array(array_of_bit_arrays);
    //     Node{children, bit_array, protein: None}
    // }

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
    fn add_child(mut curr: RwLockWriteGuard<'_, Node>, child: Arc<RwLock<Node>>) {
        // If curr is a leaf: 
        //      Create new leaf node with curr protein,
        //      Transform curr into an inner node, and 
        //      Make new leaf and child as curr's children
        if curr.children.len() == 0 {
            let cloned = Arc::new(RwLock::new(curr.clone_and_clean()));
            let child_read_locked = child.read().unwrap();
            let cloned_read_locked = cloned.read().unwrap();
            curr.bit_array = binary_or_bit_array(vec![
                child_read_locked.bit_array.clone(), 
                cloned_read_locked.bit_array.clone()]);
            curr.children.push(child.clone());
            curr.children.push(cloned.clone());
        }
        // If curr is not a leaf
        else {
            let child_read_locked = child.read().unwrap();
            let comparison_bit_array = binary_and_bit_array(
                curr.bit_array.clone(), child_read_locked.bit_array.clone());

            // If curr's proteins have at least one kmer in common with child:
            //      Find among curr's children a list of nodes with at least one kmer in common with child
            if at_least_one_true(comparison_bit_array.clone()) {
                let mut common_kmer: Vec<Arc<RwLock<Node>>> = vec![];
                let mut no_common_kmer: Vec<Arc<RwLock<Node>>> = vec![];
                let mut children_read_locks: Vec<RwLockReadGuard<'_, Node>> = vec![];

                // Divide children into two lists and keep track of read locks
                for curr_child in &curr.children {
                    children_read_locks.push(curr_child.read().unwrap());
                    let comparison_bit_array = binary_and_bit_array(
                        children_read_locks.last().unwrap().bit_array.clone(),
                        child_read_locked.bit_array.clone());
                    if at_least_one_true(comparison_bit_array) {
                        common_kmer.push(curr_child.clone());
                    }
                    else {
                        no_common_kmer.push(curr_child.clone());
                    }
                }
                // Find bit array of kmers for both lists
                let common_kmer_bit_arrays: Vec<Arc<Vec<bool>>> = common_kmer.iter()
                    .map(|x| x.read().unwrap().bit_array.clone()).collect();
                let common_kmer_bit_array = binary_or_bit_array(common_kmer_bit_arrays);
                let no_common_kmer_bit_array = {
                    if no_common_kmer.len() > 0 {
                        let no_common_kmer_bit_arrays: Vec<Arc<Vec<bool>>> = no_common_kmer.iter()
                            .map(|x| x.read().unwrap().bit_array.clone()).collect();
                        binary_or_bit_array(no_common_kmer_bit_arrays)
                    }
                    else {
                        common_kmer_bit_array.clone()
                    }
                };
                let comparison_bit_array = {
                    if no_common_kmer.len() > 0 {
                        binary_and_bit_array(common_kmer_bit_array.clone(),
                            no_common_kmer_bit_array.clone())
                    }
                    else {
                        Arc::new(vec![true])
                    }
                };
                
                // If aformentioned list of nodes also have kmer in common with remaining nodes:
                //      Just make curr the parent
                if at_least_one_true(comparison_bit_array) {
                    drop(children_read_locks);
                    curr.bit_array = binary_or_bit_array(vec![
                        child_read_locked.bit_array.clone(), curr.bit_array.clone()]);
                    curr.children.push(child.clone());
                }
                // If only one node has kmers in common with child:
                //      Update curr's bit array
                //      Call add_child on that node
                else if common_kmer.len() == 1 {
                    let new_parent = common_kmer.first().unwrap().clone();
                    drop(children_read_locks);
                    curr.bit_array = binary_or_bit_array(vec![
                        child_read_locked.bit_array.clone(), curr.bit_array.clone()]);
                    drop(curr);
                    let new_parent = new_parent.write().unwrap();
                    drop(child_read_locked);
                    Node::add_child(new_parent, child.clone());
                }

                // If aformentioned list of nodes have no kmer in common with remaining nodes:
                //      Merge nodes in common into single parent node
                //      Merge nodes in no_commin into another parent node
                //      Repeat what was done in scenario where curr had no children
                else {
                    common_kmer.push(child.clone());
                    let common_kmer_bit_array = binary_or_bit_array(vec![
                        child_read_locked.bit_array.clone(), common_kmer_bit_array.clone()]);
                    let common_node = Arc::new(RwLock::new(
                        Node {children: common_kmer, bit_array: common_kmer_bit_array, protein: None}
                    ));
                        
                    let no_common_node = {
                        if no_common_kmer.len() == 1 {
                            no_common_kmer.first().unwrap().clone()
                        }
                        else {
                            Arc::new(RwLock::new(
                                Node {
                                    children: no_common_kmer, 
                                    bit_array: no_common_kmer_bit_array,
                                    protein: None}))
                        }
                    };
                    drop(children_read_locks);
                    curr.bit_array = binary_or_bit_array(vec![
                        child_read_locked.bit_array.clone(), curr.bit_array.clone()]);
                    curr.children = vec![common_node.clone(), no_common_node.clone()];
                }
            }

            // If curr's proteins and child have no kmer in common
            else {
                let child_read_locked = child.read().unwrap();
                curr.bit_array = binary_or_bit_array(vec![
                    child_read_locked.bit_array.clone(), curr.bit_array.clone()]);
                curr.children.push(child.clone());
            }
        }
        // let array_of_bit_arrays
    }
}
impl fmt::Debug for Node {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.protein.is_some() {
            f.debug_struct("Leaf").field("protein", &self.protein).finish()
        }
        else {
            f.debug_struct("Node").field("Children", &self.children).finish()
        }
    }
}

#[derive(Debug)]
pub struct Tree {
    root: Arc<RwLock<Node>>,
    kmer_size: u8,
}
impl Tree {
    pub fn new(kmer_size: u8, protein: Arc<RwLock<Protein>>) -> Tree {
        let root = Arc::new(RwLock::new(Node::new_leaf(protein, &kmer_size)));
        Tree{
            root,
            kmer_size
        }
    }
    pub fn add_protein(&self, protein: Arc<RwLock<Protein>>) {
        let node = Arc::new(RwLock::new(Node::new_leaf(protein, &self.kmer_size)));
        let root_clone = self.root.clone();
        let root_clone = root_clone.write().unwrap();
        Node::add_child(root_clone, node);
    }
}