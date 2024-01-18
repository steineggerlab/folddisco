// pub fn interact(structure: &core::CompactStructure) {
//     for i in 0..structure.num_residues {
//         for j in i+1..structure.num_residues {

//         }
//     }
// }

pub enum HashType {
    PDBMotif,
    FoldDiscoDefault,
    Other,
}

pub trait GeometricHash {
    fn perfect_hash(feature: Vec<f32>) -> Self;
    fn reverse_hash(&self) -> Vec<f32>;
    fn hash_type(&self) -> HashType;

}
