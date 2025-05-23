use seq_io::fasta::{Record, RefRecord};
use rand::seq::index::sample;
use boomphf::*;
use std::fmt;

type BitwiseAminoAcid = u8;
type FiveMer = u32;

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
    seq: String,
    five_mers: Vec<FiveMer>,
    hash_five_mers: Vec<u32>,
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
    // Make new protein struct, but select a tenth of all 5-mers
    pub fn new_with_rand_fivemers(record: &RefRecord) -> Protein {
        // Read from fasta record
        let id = record.id().unwrap().to_string();
        let seq = str::from_utf8(record.seq()).unwrap().to_string();

        // Randomly select a tenth of all 5-mers
        let mut rng = rand::rng();
        let seq_len = seq.len();
        let five_mer_indices: Vec<usize> = sample(
            &mut rng, seq_len-4, (seq_len-4)/10).into_vec();
        let five_mers: Vec<FiveMer> = five_mer_indices.iter().map(
            |start| {
                let end = start + 5;
                let mut amino_acids = [0u8; 5];
                for i in *start..end {amino_acids[i-start] = 
                    amino_acid_to_bits(&(record.seq()[i] as char));}
                create_five_mer(amino_acids)
            }).collect();

        // Empty vector to fill out later
        let hash_five_mers: Vec<u32> = vec![];
        Protein {
            id,
            seq,
            five_mers,
            hash_five_mers,
        }
    }

    // Make new protein struct
    pub fn new(record: &RefRecord) -> Protein {
        // Read from fasta record
        let id = record.id().unwrap().to_string();
        let seq = str::from_utf8(record.seq()).unwrap().to_string();

        // Get all 5-mers in protein
        let seq_len = seq.len();
        let five_mer_indices: Vec<usize> = (0..seq_len-4).collect();
        let five_mers: Vec<FiveMer> = five_mer_indices.iter().map(
            |start| {
                let end = start + 5;
                let mut amino_acids = [0u8; 5];
                for i in *start..end {amino_acids[i-start] = 
                    amino_acid_to_bits(&(record.seq()[i] as char));}
                create_five_mer(amino_acids)
            }).collect();

        // Empty vector to fill out later
        let hash_five_mers: Vec<u32> = vec![];
        Protein {
            id,
            seq,
            five_mers,
            hash_five_mers,
        }
    }

    // Get amr class of protein
    pub fn get_amr_class(&self) -> &str {
        let protein_attr: Vec<&str> = self.id.split_terminator('|').collect();
        protein_attr[3]
    }

    // Get 5-mers in protein struct
    pub fn get_five_mers(&self) -> Vec<FiveMer> {
        return self.five_mers.clone();
    }

    // Get hash keys for the 5-mers in protein struct
    pub fn get_five_hash(&self) -> Vec<u32> {
        return self.hash_five_mers.clone();
    }

    // Remove 5-mers that only appear in one protein
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

    // Populate hash keys vector
    pub fn modify_hash_five_mer(&mut self, len: &usize, phf: &Mphf<u32>) {
        let mut hash_map = vec![false; *len];
        for five_mer in &self.five_mers {
            let curr_hash = phf.hash(five_mer) as usize;
            if !hash_map[curr_hash] {
                hash_map[curr_hash] = true;
                self.hash_five_mers.push(curr_hash as u32);
            }
        }
    }

    // Get protein id and sequence
    pub fn get_id_and_seq(&self) -> (&String, &String) {
        return (&self.id, &self.seq);
    }
}
