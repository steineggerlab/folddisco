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
