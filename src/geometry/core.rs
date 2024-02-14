// File: core.rs
// Author: Hyunbin Kim (khb7840@gmail.com)
// Description: Core geometric hash traits and types

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
