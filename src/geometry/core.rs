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
    Default32bit,
    PointPairFeature,
    // append new hash type here
    Other,
}

impl HashType {
    pub fn get_with_index(index: usize) -> Self {
        match index {
            0 => HashType::PDBMotif,
            1 => HashType::PDBMotifSinCos,
            2 => HashType::TrRosetta,
            3 => HashType::FoldDiscoDefault,
            4 => HashType::Default32bit,
            5 => HashType::PointPairFeature,
            // append new hash type here
            _ => HashType::Other,
        }
    }
    pub fn get_with_str(hash_type: &str) -> Self {
        match hash_type {
            "0" | "PDBMotif" | "pyscomotif" => HashType::PDBMotif,
            "1" | "PDBMotifSinCos" | "pdb" => HashType::PDBMotifSinCos,
            "2" | "TrRosetta" | "trrosetta" => HashType::TrRosetta,
            "3" | "FoldDiscoDefault" | "default" => HashType::FoldDiscoDefault,
            "4" | "Default32bit" | "default32" => HashType::Default32bit,
            "5" | "PointPairFeature" | "ppf" => HashType::PointPairFeature,
            // append new hash type here
            _ => HashType::Other,
        }
    }
    pub fn encoding_type(&self) -> usize {
        match self {
            HashType::PDBMotif => 32usize,
            HashType::PDBMotifSinCos => 32usize,
            HashType::TrRosetta => 64usize,
            HashType::FoldDiscoDefault => 64usize,
            HashType::Default32bit => 32usize,
            HashType::PointPairFeature => 32usize,
            // append new hash type here
            HashType::Other => 0usize,
        }
    }
    
    pub fn save_to_file(&self, path: &str) {
        let mut file = std::fs::File::create(path).unwrap();
        file.write_all(format!("{:?}", self).as_bytes()).unwrap();
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
                "Default32bit" => HashType::Default32bit,
                "PointPairFeature" => HashType::PointPairFeature,
                // append new hash type here
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
        let path = "data/test.type";
        let hash_type_vec = vec![
            HashType::PDBMotif,
            HashType::PDBMotifSinCos,
            HashType::TrRosetta,
            HashType::FoldDiscoDefault,
            HashType::Default32bit,
            HashType::PointPairFeature,
            // append new hash type here
        ];
        for hash_type in hash_type_vec {
            hash_type.save_to_file(path);
            let loaded_hash_type = HashType::load_from_file(path);
            assert_eq!(hash_type, loaded_hash_type);
        }
    }
}

#[derive(Clone, Copy, Hash, Eq, PartialEq, Ord, PartialOrd)]
pub enum GeometricHash {
    PDBMotif(super::pdb_motif::HashValue),
    PDBMotifSinCos(super::pdb_motif_sincos::HashValue),
    TrRosetta(super::trrosetta::HashValue),
    FoldDiscoDefault(super::default::HashValue),
    Default32bit(super::default_32bit::HashValue),
    PointPairFeature(super::ppf::HashValue),
    // append new hash type here
}

impl HashableSync for GeometricHash {}

impl GeometricHash {
    pub fn perfect_hash_default(feature: Vec<f32>, hash_type: HashType) -> Self {
        match hash_type {
            HashType::PDBMotif => GeometricHash::PDBMotif(
                super::pdb_motif::HashValue::perfect_hash_default(feature)
            ),
            HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                super::pdb_motif_sincos::HashValue::perfect_hash_default(feature)
            ),
            HashType::TrRosetta => GeometricHash::TrRosetta(
                super::trrosetta::HashValue::perfect_hash_default(feature)
            ),
            HashType::FoldDiscoDefault => GeometricHash::FoldDiscoDefault(
                super::default::HashValue::perfect_hash_default(feature)
            ),
            HashType::Default32bit => GeometricHash::Default32bit(
                super::default_32bit::HashValue::perfect_hash_default(feature)
            ),
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue::perfect_hash_default(feature)
            ),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn perfect_hash(
        feature: Vec<f32>, hash_type: HashType, nbin_dist: usize, nbin_angle: usize
    ) -> Self {
            match hash_type {
                HashType::PDBMotif => GeometricHash::PDBMotif(
                    super::pdb_motif::HashValue::perfect_hash(
                        feature, nbin_dist, nbin_angle
                    )
                ),
                HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                    super::pdb_motif_sincos::HashValue::perfect_hash(
                        feature, nbin_dist, nbin_angle
                    )
                ),
                HashType::TrRosetta => GeometricHash::TrRosetta(
                    super::trrosetta::HashValue::perfect_hash(
                        feature, nbin_dist, nbin_angle
                    )
                ),
                HashType::FoldDiscoDefault => GeometricHash::FoldDiscoDefault(
                    super::default::HashValue::perfect_hash(
                        feature, nbin_dist, nbin_angle
                    )
                ),
                HashType::Default32bit => GeometricHash::Default32bit(
                    super::default_32bit::HashValue::perfect_hash(
                        feature, nbin_dist, nbin_angle
                    )
                ),
                HashType::PointPairFeature => GeometricHash::PointPairFeature(
                    super::ppf::HashValue::perfect_hash(
                        feature, nbin_dist, nbin_angle
                    )
                ),
                // append new hash type here
                _ => panic!("Invalid hash type"),
            }
    }
    
    
    pub fn reverse_hash_default(&self) -> Vec<f32> {
        match self {
            GeometricHash::PDBMotif(hash) => hash.reverse_hash_default(),
            GeometricHash::PDBMotifSinCos(hash) => hash.reverse_hash_default(),
            GeometricHash::TrRosetta(hash) => hash.reverse_hash_default(),
            GeometricHash::FoldDiscoDefault(hash) => hash.reverse_hash_default(),
            GeometricHash::Default32bit(hash) => hash.reverse_hash_default(),
            GeometricHash::PointPairFeature(hash) => hash.reverse_hash_default(),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }

    pub fn reverse_hash(&self, nbin_dist: usize, nbin_angle: usize) -> Vec<f32> {
        match self {
            GeometricHash::PDBMotif(hash) => hash.reverse_hash(nbin_dist, nbin_angle),
            GeometricHash::PDBMotifSinCos(hash) => hash.reverse_hash(nbin_dist, nbin_angle),
            GeometricHash::TrRosetta(hash) => hash.reverse_hash(nbin_dist, nbin_angle),
            GeometricHash::FoldDiscoDefault(hash) => hash.reverse_hash(nbin_dist, nbin_angle),
            GeometricHash::Default32bit(hash) => hash.reverse_hash(nbin_dist, nbin_angle),
            GeometricHash::PointPairFeature(hash) => hash.reverse_hash(nbin_dist, nbin_angle),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn hash_type(&self) -> HashType {
        match self {
            GeometricHash::PDBMotif(hash) => hash.hash_type(),
            GeometricHash::PDBMotifSinCos(hash) => hash.hash_type(),
            GeometricHash::TrRosetta(hash) => hash.hash_type(),
            GeometricHash::FoldDiscoDefault(hash) => hash.hash_type(),
            GeometricHash::Default32bit(hash) => hash.hash_type(),
            GeometricHash::PointPairFeature(hash) => hash.hash_type(),
            // append new hash type here
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
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue::from_u32(hashvalue)
            ),
            HashType::Default32bit => GeometricHash::Default32bit(
                super::default_32bit::HashValue::from_u32(hashvalue)
            ),
            // append new hash type here if it is encoded as u32
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
            HashType::Default32bit => GeometricHash::Default32bit(
                super::default_32bit::HashValue::from_u64(hashvalue)
            ),
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue::from_u64(hashvalue)
            ),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn as_u32(&self) -> u32 {
        match self {
            GeometricHash::PDBMotif(hash) => hash.as_u32(),
            GeometricHash::PDBMotifSinCos(hash) => hash.as_u32(),
            GeometricHash::PointPairFeature(hash) => hash.as_u32(),
            GeometricHash::Default32bit(hash) => hash.as_u32(),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn as_u64(&self) -> u64 {
        match self {
            GeometricHash::PDBMotif(hash) => hash.as_u64(),
            GeometricHash::PDBMotifSinCos(hash) => hash.as_u64(),
            GeometricHash::TrRosetta(hash) => hash.as_u64(),
            GeometricHash::FoldDiscoDefault(hash) => hash.as_u64(),
            GeometricHash::Default32bit(hash) => hash.as_u64(),
            GeometricHash::PointPairFeature(hash) => hash.as_u64(),
            // append new hash type here
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
    pub fn downcast_default_32bit(&self) -> super::default_32bit::HashValue {
        match self {
            GeometricHash::Default32bit(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_point_pair_feature(&self) -> super::ppf::HashValue {
        match self {
            GeometricHash::PointPairFeature(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    // append the downcast method for new hash type here

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
            GeometricHash::Default32bit(hash) => {
                write!(f, "Default32bit({:?})", hash)
            },
            GeometricHash::PointPairFeature(hash) => {
                write!(f, "PointPairFeature({:?})", hash)
            },
            // append new hash type here
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
            GeometricHash::Default32bit(hash) => {
                write!(f, "Default32bit\t{:?}", hash)
            },
            GeometricHash::PointPairFeature(hash) => {
                write!(f, "PointPairFeature\t{:?}", hash)
            },
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }
}