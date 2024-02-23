// File: core.rs
// Author: Hyunbin Kim (khb7840@gmail.com)
// Description: Core geometric hash enum and types

use std::{fmt, io::{BufRead, Write}};
use crate::HashableSync;

#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub enum HashType {
    PDBMotif,
    PDBMotifSinCos,
    TrRosetta,
    FoldDiscoDefault,
    Other,
}

impl HashType {
    pub fn get_with_index(index: usize) -> Self {
        match index {
            0 => HashType::PDBMotif,
            1 => HashType::PDBMotifSinCos,
            2 => HashType::TrRosetta,
            3 => HashType::FoldDiscoDefault,
            _ => HashType::Other,
        }
    }
    pub fn get_with_str(hash_type: &str) -> Self {
        match hash_type {
            "0" | "PDBMotif" => HashType::PDBMotif,
            "1" | "PDBMotifSinCos" | "pdb" => HashType::PDBMotifSinCos,
            "2" | "TrRosetta" | "trrosetta" => HashType::TrRosetta,
            "3" | "FoldDiscoDefault" | "default" => HashType::FoldDiscoDefault,
            _ => HashType::Other,
        }
    }
    pub fn encoding_type(&self) -> usize {
        match self {
            HashType::PDBMotif => 32usize,
            HashType::PDBMotifSinCos => 32usize,
            HashType::TrRosetta => 64usize,
            HashType::FoldDiscoDefault => 64usize,
            HashType::Other => 0usize,
        }
    }
    
    pub fn save_to_file(&self, path: &str) {
        let mut file = std::fs::File::create(path).unwrap();
        file.write_all(format!("{:?}", self).as_bytes()).unwrap();
        println!("{:?}", self);
    }

    pub fn load_from_file(path: &str) -> Self {
        let file = std::fs::File::open(path).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut hash_type = HashType::FoldDiscoDefault;
        for line in reader.lines() {
            let line = line.unwrap();
            hash_type = match line.as_str() {
                "PDBMotif" => HashType::PDBMotif,
                "PDBMotifSinCos" => HashType::PDBMotifSinCos,
                "TrRosetta" => HashType::TrRosetta,
                "FoldDiscoDefault" => HashType::FoldDiscoDefault,
                _ => HashType::Other,
            };
        }
        hash_type
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash_type() {
        let path = "test_hash_type.txt";
        let hash_type = HashType::PDBMotif;
        hash_type.save_to_file(path);
        let hash_type = HashType::PDBMotifSinCos;
        hash_type.save_to_file(path);
        let hash_type = HashType::TrRosetta;
        hash_type.save_to_file(path);
        let hash_type = HashType::FoldDiscoDefault;
        hash_type.save_to_file(path);
    }
}

#[derive(Clone, Copy, Hash, Eq, PartialEq, Ord, PartialOrd)]
pub enum GeometricHash {
    PDBMotif(super::pdb_motif::HashValue),
    PDBMotifSinCos(super::pdb_motif_sincos::HashValue),
    TrRosetta(super::trrosetta::HashValue),
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
            HashType::TrRosetta => GeometricHash::TrRosetta(
                super::trrosetta::HashValue::perfect_hash(feature)
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
            GeometricHash::TrRosetta(hash) => hash.reverse_hash(),
            GeometricHash::FoldDiscoDefault(hash) => hash.reverse_hash(),
            _ => panic!("Invalid hash type"),
        }
    }

    pub fn hash_type(&self) -> HashType {
        match self {
            GeometricHash::PDBMotif(hash) => hash.hash_type(),
            GeometricHash::PDBMotifSinCos(hash) => hash.hash_type(),
            GeometricHash::TrRosetta(hash) => hash.hash_type(),
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
            HashType::TrRosetta => GeometricHash::TrRosetta(
                super::trrosetta::HashValue::from_u64(hashvalue)
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
            GeometricHash::TrRosetta(hash) => hash.as_u64(),
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
    pub fn downcast_trrosetta(&self) -> super::trrosetta::HashValue {
        match self {
            GeometricHash::TrRosetta(hash) => hash.clone(),
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
            GeometricHash::TrRosetta(hash) => {
                write!(f, "TrRosetta({:?})", hash)
            },
            GeometricHash::FoldDiscoDefault(hash) => {
                write!(f, "FoldDiscoDefault({:?})", hash)
            },
            _ => panic!("Invalid hash type"),
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
            GeometricHash::TrRosetta(hash) => {
                write!(f, "TrRosetta\t{:?}", hash)
            },
            GeometricHash::FoldDiscoDefault(hash) => {
                write!(f, "FoldDiscoDefault\t{:?}", hash)
            },
            _ => panic!("Invalid hash type"),
        }
    }
}