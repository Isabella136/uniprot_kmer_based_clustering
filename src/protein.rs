use seq_io::fasta::{Record, RefRecord};
// use rand::seq::index::sample;
use boomphf::*;
use std::fmt;

type BitwiseAminoAcid = u8;
type FiveMer = u32;
type SevenMer = u32;

const AMINO_ACID_LIST: [char; 21] = [
    'C', 'S', 'T', 'A', 'G',
    'P', 'D', 'E', 'Q', 'N',
    'H', 'R', 'K', 'M', 'I', 
    'L', 'V', 'W', 'Y', 'F', '*'];

fn bitwise_power(mut a: u32, mut n: u32) -> u32 {
    let mut ans: u32 = 1;
    while n > 0 {
        let odd: u32 = n & 1u32;        // is 1 if n is odd, else is zero
        if odd > 0u32 {                 // if odd
            ans *= a;                   // can multiply
        }
        if n > 1 {
            a = a * a;                  // raise a by 2
        }
        n = n >> 1;                     // floor divide n by 2
    }
    ans
}
fn create_five_mer(amino_acids: [BitwiseAminoAcid; 5]) -> FiveMer {
    let powers_of_twenty_one = [
        1, 21, 21*21, bitwise_power(21, 3), bitwise_power(21, 4)];
    let mut five_mer: u32= 0u32;
    for i in 0..5 {
        five_mer += Into::<u32>::into(amino_acids[i]) * powers_of_twenty_one[4-i];
    }
    five_mer
}
fn five_mer_back_to_amino_acid(five_mer: FiveMer) -> String {
    let powers_of_twenty_one = [
        1, 21, 21*21, bitwise_power(21, 3), bitwise_power(21, 4)];
    let mut amino_acids = String::new();
    let mut five_mer_clone = five_mer.clone();
    for i in 0..5 {
        amino_acids.push(AMINO_ACID_LIST[(five_mer_clone / powers_of_twenty_one[4-i]) as usize]);
        five_mer_clone %= powers_of_twenty_one[4-i];
    }
    amino_acids
}
fn create_seven_mer(amino_acids: [BitwiseAminoAcid; 7]) -> SevenMer {
    let powers_of_twenty_one = [
        1, 21, 21*21, bitwise_power(21, 3), bitwise_power(21, 4),
        bitwise_power(21, 5), bitwise_power(21, 6)];
    let mut seven_mer: u32= 0u32;
    for i in 0..7 {
        seven_mer += Into::<u32>::into(amino_acids[i]) * powers_of_twenty_one[6-i];
    }
    seven_mer
}
fn amino_acid_to_bits(amino_acid: &char) -> BitwiseAminoAcid {
    let index_usize = AMINO_ACID_LIST.iter().position(|x| x==amino_acid)
        .unwrap_or(20usize);
    let index_u8: u8 = index_usize.try_into().unwrap();
    index_u8
}


#[derive(Default)]
#[derive(Clone)]
pub struct Protein {
    id: String,
    five_mers: Vec<FiveMer>,
    seven_mers: Vec<SevenMer>,
    hash_five_mers: Vec<u32>,
    hash_seven_mers: Vec<u32>,
    hashmap_five_mers: Vec<bool>,
    hashmap_seven_mers: Vec<bool>,
}
impl fmt::Debug for Protein {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let five_mers: Vec<String> =  self.five_mers.iter()
                .map(|x| five_mer_back_to_amino_acid(*x)).collect();
        f.debug_struct("Protein")
            .field("id", &self.id)
            .field("5-mers", &five_mers)
            .finish()
    }
}
impl Protein {
    pub fn new_protein(record: &RefRecord) -> Protein {
        // let mut rng = rand::rng();
        let id = record.id().unwrap().to_string();
        let seq = record.seq();
        let seq_len = seq.len();
        let five_mer_indices: Vec<usize> = (0..seq_len-4).collect();
        // let five_mer_indices: Vec<usize> = sample(
        //     &mut rng, seq_len-4, (seq_len-4)/10).into_vec();
        let seven_mer_indices: Vec<usize> = (0..seq_len-6).collect();
        // let seven_mer_indices: Vec<usize> = sample(
        //     &mut rng, seq_len-6, (seq_len-6)/4).into_vec();
        let five_mers: Vec<FiveMer> = five_mer_indices.iter().map(
            |start| {
                let end = start + 5;
                let mut amino_acids = [0u8; 5];
                for i in *start..end {
                    amino_acids[i-start] = amino_acid_to_bits(&(seq[i] as char));
                }
                create_five_mer(amino_acids)
            }).collect();
        let seven_mers: Vec<SevenMer> = seven_mer_indices.iter().map(
            |start| {
                let end = start + 7;
                let mut amino_acids = [0u8; 7];
                for i in *start..end {
                    amino_acids[i-start] = amino_acid_to_bits(&(seq[i] as char));
                }
                create_seven_mer(amino_acids)
            }).collect();
        let hash_five_mers: Vec<u32> = vec![];
        let hash_seven_mers: Vec<u32> = vec![];
        let hashmap_five_mers: Vec<bool> = vec![];
        let hashmap_seven_mers: Vec<bool> = vec![];
        Protein {
            id,
            five_mers,
            seven_mers,
            hash_five_mers,
            hash_seven_mers,
            hashmap_five_mers,
            hashmap_seven_mers,
        }
    }
    pub fn get_five_mers(&self) -> Vec<FiveMer> {
        return self.five_mers.clone();
    }
    pub fn get_seven_mers(&self) -> Vec<SevenMer> {
        return self.seven_mers.clone();
    }
    pub fn get_five_hash_map(&self) -> Vec<bool> {
        return self.hashmap_five_mers.clone();
    }
    pub fn get_five_hash(&self) -> Vec<u32> {
        return self.hash_five_mers.clone();
    }
    pub fn get_seven_hash_map(&self) -> Vec<bool> {
        return self.hashmap_seven_mers.clone();
    }
    pub fn get_seven_hash(&self) -> Vec<u32> {
        return self.hash_seven_mers.clone();
    }
    pub fn remove_unique_five_mers(&mut self, phf: &Mphf<u32>, unique_kmers: &Vec<bool>) {
        let mut index = 0usize;
        while index < self.five_mers.len() {
            while unique_kmers[phf.hash(&self.five_mers[index]) as usize] {
                self.five_mers.remove(index);
                if index == self.five_mers.len() {
                    break;
                }
            }
            index += 1;
        }
    }
    pub fn remove_unique_seven_mers(&mut self, phf: &Mphf<u32>, unique_kmers: &Vec<bool>) {
        let mut index = 0usize;
        while index < self.seven_mers.len() {
            while unique_kmers[phf.hash(&self.seven_mers[index]) as usize] {
                self.seven_mers.remove(index);
                if index == self.seven_mers.len() {
                    break;
                }
            }
            index += 1;
        }
    }
    pub fn modify_hash_five_mer(&mut self, len: &usize, phf: &Mphf<u32>) {
        let mut hash_map = vec![false; *len];
        for five_mer in &self.five_mers {
            let curr_hash = phf.hash(five_mer) as usize;
            hash_map[curr_hash] = true;
            self.hash_five_mers.push(curr_hash as u32);
        }
        self.hashmap_five_mers = hash_map;
    }
    pub fn modify_hash_seven_mer(&mut self, len: &usize, phf: &Mphf<u32>) {
        let mut hash_map = vec![false; *len];
        for seven_mer in &self.seven_mers {
            let curr_hash = phf.hash(seven_mer) as usize;
            hash_map[curr_hash] = true;
            self.hash_seven_mers.push(curr_hash as u32);
        }
        self.hashmap_seven_mers = hash_map;
    }
}
