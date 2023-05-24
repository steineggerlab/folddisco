// pub fn interact(structure: &core::CompactStructure) {
//     for i in 0..structure.num_residues {
//         for j in i+1..structure.num_residues {

//         }
//     }
// }


pub enum HashType {
    SimpleHash,
    TriadHash,
    PPFHash,
    TRRosettaHash,
    TRRosettaReducedHash,
    None,
}

pub trait GeometricHash {
    // fn hash(&self, structure: &core::CompactStructure) -> Vec<u64>;
    fn hash_type(&self) -> HashType;
}