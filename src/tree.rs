use std::ops::Deref;

use crate::protein::Protein;

fn binary_or_bit_array(bit_arrays: Vec<Box<Vec<bool>>>) -> Box<Vec<bool>> {
    let length = bit_arrays.len();
    if length == 1 {
        bit_arrays[0].clone()
    }
    else if length == 2 {
        Box::new(bit_arrays[0].iter().zip(bit_arrays[1].iter())
            .map(|(x,y)| x | y).collect())
    }
    else {
        Box::new(binary_or_bit_array(bit_arrays[..length/2].to_vec()).iter()
            .zip(binary_or_bit_array(bit_arrays[length/2..].to_vec()).iter())
            .map(|(x,y)| x | y).collect())
    }
}
fn binary_and_bit_array(first_bit_array: &Vec<bool>, second_bit_array: &Vec<bool>) -> Box<Vec<bool>> {
    Box::new(first_bit_array.iter()        .zip(second_bit_array.iter())
        .map(|(x,y)| x | y).collect())
}

fn at_least_one_true(bit_array: Box<Vec<bool>>) -> bool {
    for bit in bit_array.iter() {
        if *bit {
            return true
        }
    }
    false
}

#[derive(Clone)]
struct Node {
    children: Vec<Box<Node>>,
    bit_array: Vec<bool>,
}
impl Node {
    fn new_leaf(protein: &Protein, kmer_size: &u8) -> Node {
        let mut node = Node {
            children: vec![],
            bit_array: Vec::new(),
        };
        node.add_bit_array_from_protein(protein, kmer_size);
        node
    }

    fn add_bit_array_from_protein(&mut self, protein: &Protein, kmer_size: &u8) {
        if !self.children.is_empty() {
            panic!("add_bit_array_from_protein should only be called on new leaf");
        }
        if *kmer_size == 5u8 {
            self.bit_array = protein.get_five_hash();
        }
        else if *kmer_size == 7u8 {
            self.bit_array = protein.get_seven_hash();
        }
        else {
            panic!("{}", format!("kmer-size should be 5 or 7, not {kmer_size}"));
        }
    }

    fn new_node_with_children(children: Vec<Box<Node>>) -> Node {
        let array_of_bit_arrays: Vec<Box<Vec<bool>>> = children.iter()
            .map(|x| Box::new(x.bit_array.clone())).collect();
        let bit_array = binary_or_bit_array(array_of_bit_arrays);
        Node{children, bit_array: bit_array.deref().to_vec()}
    }

    fn clone_and_clean(&mut self) -> Node {
        let mut children_copy: Vec<Box<Node>> = vec![];
        while self.children.len() > 0 {
            children_copy.push(self.children.pop().unwrap());
        }
        let bit_array_copy = self.bit_array.clone();
        self.bit_array = Vec::new();
        Node {children: children_copy, bit_array: bit_array_copy}
        
    }

    fn add_child(&mut self, child: Box<Node>) {
        if self.children.len() == 0 {
            let cloned = Box::new(self.clone_and_clean());
            self.bit_array = binary_or_bit_array(vec![
                Box::new(child.bit_array.clone()), 
                Box::new(cloned.bit_array.clone())]).deref().to_vec();
            self.children.push(child);
            self.children.push(cloned);
        }
        else {
            let comparison_bit_array = binary_and_bit_array(
                &self.bit_array, &child.bit_array);
            if at_least_one_true(comparison_bit_array) {
                let mut common_kmer: Vec<&Box<Node>> = vec![];
                let mut no_common_kmer: Vec<&Box<Node>> = vec![];
                for curr_child in &self.children {
                    let comparison_bit_array = binary_and_bit_array(
                    &curr_child.bit_array, &child.bit_array);
                    if at_least_one_true(comparison_bit_array) {
                        common_kmer.push(curr_child);
                    }
                    else {
                        no_common_kmer.push(curr_child);
                    }
                }
                let common_kmer_bit_arrays: Vec<Box<Vec<bool>>> = common_kmer.iter()
                    .map(|x| Box::new(x.bit_array.clone())).collect();
                let common_kmer_bit_array = binary_or_bit_array(common_kmer_bit_arrays);
                let no_common_kmer_bit_arrays: Vec<Box<Vec<bool>>> = no_common_kmer.iter()
                    .map(|x| Box::new(x.bit_array.clone())).collect();
                let no_common_kmer_bit_array = binary_or_bit_array(no_common_kmer_bit_arrays);
                let comparison_bit_array = binary_and_bit_array(
                    &*common_kmer_bit_array, &*no_common_kmer_bit_array);
                if at_least_one_true(comparison_bit_array) {
                    self.bit_array = binary_or_bit_array(vec![
                        Box::new(child.bit_array.clone()), 
                        Box::new(self.bit_array.clone())]).deref().to_vec();
                    self.children.push(child);
                }
                else {
                    // Merge nodes in common into single parent node
                    // Merge nodes in no_commin into another parent node
                    // Repeat what was done in scenario where self had no children
                }
            }
            else {
                self.bit_array = binary_or_bit_array(vec![
                    Box::new(child.bit_array.clone()), 
                    Box::new(self.bit_array.clone())]).deref().to_vec();
                self.children.push(child);
            }
        }
        // let array_of_bit_arrays
    }
}

struct Tree {
    root: Box<Node>,
    kmer_size: u8,
}
impl Tree {
    fn new(kmer_size: u8, protein: Protein) -> Tree {
        let root = Box::new(Node::new_leaf(&protein, &kmer_size));
        Tree{
            root,
            kmer_size
        }
    }
    fn add_protein(&mut self, protein: Protein) {
        let node = Box::new(Node::new_leaf(&protein, &self.kmer_size));
        self.root.add_child(node);
    }
}