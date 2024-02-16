// File: core.rs
// Author: Hyunbin Kim (khb7840@gmail.com)
// Description: Core geometric hash enum and types

use std::fmt;
use crate::HashableSync;

#[derive(Clone, Copy, Eq, PartialEq)]
pub enum HashType {
    PDBMotif,
    PDBMotifSinCos,
    FoldDiscoDefault,
    Other,
}

#[derive(Clone, Copy, Hash, Eq, PartialEq, Ord, PartialOrd)]
pub enum GeometricHash {
    PDBMotif(super::pdb_motif::HashValue),
    PDBMotifSinCos(super::pdb_motif_sincos::HashValue),
    FoldDiscoDefault(super::default::HashValue),
}

impl HashableSync for GeometricHash {}

impl GeometricHash {
    pub fn perfect_hash(feature: Vec<f32>, hash_type: HashType) -> Self {
        match hash_type {
            HashType::PDBMotif => GeometricHash::PDBMotif(
                super::pdb_motif::HashValue::perfect_hash(feature)
            ),
            HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                super::pdb_motif_sincos::HashValue::perfect_hash(feature)
            ),
            HashType::FoldDiscoDefault => GeometricHash::FoldDiscoDefault(
                super::default::HashValue::perfect_hash(feature)
            ),
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn reverse_hash(&self) -> Vec<f32> {
        match self {
            GeometricHash::PDBMotif(hash) => hash.reverse_hash(),
            GeometricHash::PDBMotifSinCos(hash) => hash.reverse_hash(),
            GeometricHash::FoldDiscoDefault(hash) => hash.reverse_hash(),
            _ => panic!("Invalid hash type"),
        }
    }

    pub fn hash_type(&self) -> HashType {
        match self {
            GeometricHash::PDBMotif(hash) => hash.hash_type(),
            GeometricHash::PDBMotifSinCos(hash) => hash.hash_type(),
            GeometricHash::FoldDiscoDefault(hash) => hash.hash_type(),
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn from_u32(hashvalue: u32, hash_type: HashType) -> Self {
        match hash_type {
            HashType::PDBMotif => GeometricHash::PDBMotif(
                super::pdb_motif::HashValue::from_u32(hashvalue)
            ),
            HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                super::pdb_motif_sincos::HashValue::from_u32(hashvalue)
            ),
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn from_u64(hashvalue: u64, hash_type: HashType) -> Self {
        match hash_type {
            HashType::PDBMotif => GeometricHash::PDBMotif(
                super::pdb_motif::HashValue::from_u64(hashvalue)
            ),
            HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                super::pdb_motif_sincos::HashValue::from_u64(hashvalue)
            ),
            HashType::FoldDiscoDefault => GeometricHash::FoldDiscoDefault(
                super::default::HashValue::from_u64(hashvalue)
            ),
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn as_u32(&self) -> u32 {
        match self {
            GeometricHash::PDBMotif(hash) => hash.as_u32(),
            GeometricHash::PDBMotifSinCos(hash) => hash.as_u32(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn as_u64(&self) -> u64 {
        match self {
            GeometricHash::PDBMotif(hash) => hash.as_u64(),
            GeometricHash::PDBMotifSinCos(hash) => hash.as_u64(),
            GeometricHash::FoldDiscoDefault(hash) => hash.as_u64(),
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn downcast_pdb_motif(&self) -> super::pdb_motif::HashValue {
        match self {
            GeometricHash::PDBMotif(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_pdb_motif_sincos(&self) -> super::pdb_motif_sincos::HashValue {
        match self {
            GeometricHash::PDBMotifSinCos(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_fold_disco_default(&self) -> super::default::HashValue {
        match self {
            GeometricHash::FoldDiscoDefault(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    
}

impl fmt::Debug for GeometricHash {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GeometricHash::PDBMotif(hash) => {
                write!(f, "PDBMotif({:?})", hash)
            },
            GeometricHash::PDBMotifSinCos(hash) => {
                write!(f, "PDBMotifSinCos({:?})", hash)
            },
            GeometricHash::FoldDiscoDefault(hash) => {
                write!(f, "FoldDiscoDefault({:?})", hash)
            },
        }
    }
}

impl fmt::Display for GeometricHash {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GeometricHash::PDBMotif(hash) => {
                write!(f, "PDBMotif\t{:?}", hash)
            },
            GeometricHash::PDBMotifSinCos(hash) => {
                write!(f, "PDBMotifSinCos\t{:?}", hash)
            },
            GeometricHash::FoldDiscoDefault(hash) => {
                write!(f, "FoldDiscoDefault\t{:?}", hash)
            },
        }
    }
}