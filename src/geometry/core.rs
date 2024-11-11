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
    PDBTrRosetta,
    PointPairFeature,
    TertiaryInteraction,
    // append new hash type here
    Other,
}

impl HashType {
    
    pub fn get_with_index(index: usize) -> Self {
        match index {
            0 => HashType::PDBMotif,
            1 => HashType::PDBMotifSinCos,
            2 => HashType::TrRosetta,
            3 => HashType::PDBTrRosetta,
            4 => HashType::PointPairFeature,
            5 => HashType::TertiaryInteraction,
            // append new hash type here
            _ => HashType::Other,
        }
    }
    
    pub fn get_with_str(hash_type: &str) -> Self {
        match hash_type {
            "0" | "PDBMotif" | "pyscomotif" | "orig_pdb" => HashType::PDBMotif,
            "1" | "PDBMotifSinCos" | "pdb" => HashType::PDBMotifSinCos,
            "2" | "TrRosetta" | "trrosetta" | "tr" => HashType::TrRosetta,
            "3" | "PDBTrRosetta" | "pdbtr" | "default" | "folddisco" => HashType::PDBTrRosetta,
            "4" | "PointPairFeature" | "ppf" => HashType::PointPairFeature,
            "5" | "TertiaryInteraction" | "tertiary" | "3di" => HashType::TertiaryInteraction,
            // append new hash type here
            _ => HashType::Other,
        }
    }
    
    pub fn to_string(&self) -> String {
        match self {
            HashType::PDBMotif => "PDBMotif".to_string(),
            HashType::PDBMotifSinCos => "PDBMotifSinCos".to_string(),
            HashType::TrRosetta => "TrRosetta".to_string(),
            HashType::PDBTrRosetta => "PDBTrRosetta".to_string(),
            HashType::PointPairFeature => "PointPairFeature".to_string(),
            HashType::TertiaryInteraction => "TertiaryInteraction".to_string(),
            // append new hash type here
            HashType::Other => "Other".to_string(),
        }
    }

    pub fn encoding_type(&self) -> usize {
        // Unified to u32 encoding
        32usize
    }

    pub fn encoding_bits(&self) -> usize {
        match self {
            HashType::PDBMotif => 25usize,
            HashType::PDBMotifSinCos => 26usize,
            HashType::TrRosetta => 32usize,
            HashType::PDBTrRosetta => 30usize,
            HashType::PointPairFeature => 32usize,
            HashType::TertiaryInteraction => 29usize,
            // append new hash type here
            HashType::Other => 32usize,
        }
    }
    
    pub fn save_to_file(&self, path: &str) {
        let mut file = std::fs::File::create(path).unwrap();
        file.write_all(format!("{:?}", self).as_bytes()).unwrap();
    }

    pub fn load_from_file(path: &str) -> Self {
        let file = std::fs::File::open(path).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut hash_type = HashType::PDBTrRosetta;
        for line in reader.lines() {
            let line = line.unwrap();
            hash_type = match line.as_str() {
                "PDBMotif" => HashType::PDBMotif,
                "PDBMotifSinCos" => HashType::PDBMotifSinCos,
                "TrRosetta" => HashType::TrRosetta,
                "PDBTrRosetta" => HashType::PDBTrRosetta,
                "PointPairFeature" => HashType::PointPairFeature,
                "TertiaryInteraction" => HashType::TertiaryInteraction,
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
            HashType::PDBTrRosetta,
            HashType::PointPairFeature,
            HashType::TertiaryInteraction,
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
    PDBTrRosetta(super::pdb_tr::HashValue),
    PointPairFeature(super::ppf::HashValue),
    TertiaryInteraction(super::tertiary_interaction::HashValue),
    // append new hash type here
}

impl HashableSync for GeometricHash {}

impl GeometricHash {
    pub fn perfect_hash_default_as_u32(feature: &Vec<f32>, hash_type: HashType) -> u32 {
        match hash_type {
            HashType::PDBMotif => super::pdb_motif::HashValue::perfect_hash_default(feature),
            HashType::PDBMotifSinCos => super::pdb_motif_sincos::HashValue::perfect_hash_default(feature),
            HashType::TrRosetta => super::trrosetta::HashValue::perfect_hash_default(feature),
            HashType::PDBTrRosetta => super::pdb_tr::HashValue::perfect_hash_default(feature),
            HashType::PointPairFeature => super::ppf::HashValue::perfect_hash_default(feature),
            HashType::TertiaryInteraction => super::tertiary_interaction::HashValue::perfect_hash_default(feature),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }

    pub fn perfect_hash_as_u32(
        feature: &Vec<f32>, hash_type: HashType, nbin_dist: usize, nbin_angle: usize
    ) -> u32 {
        match hash_type {
            HashType::PDBMotif => super::pdb_motif::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::PDBMotifSinCos => super::pdb_motif_sincos::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::TrRosetta => super::trrosetta::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::PDBTrRosetta => super::pdb_tr::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::PointPairFeature => super::ppf::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            HashType::TertiaryInteraction => super::tertiary_interaction::HashValue::perfect_hash(
                feature, nbin_dist, nbin_angle
            ),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn perfect_hash_default(feature: &Vec<f32>, hash_type: HashType) -> Self {
        match hash_type {
            HashType::PDBMotif => GeometricHash::PDBMotif(
                super::pdb_motif::HashValue(
                    super::pdb_motif::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                super::pdb_motif_sincos::HashValue(
                    super::pdb_motif_sincos::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::TrRosetta => GeometricHash::TrRosetta(
                super::trrosetta::HashValue(
                    super::trrosetta::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::PDBTrRosetta => GeometricHash::PDBTrRosetta(
                super::pdb_tr::HashValue(
                    super::pdb_tr::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue(
                    super::ppf::HashValue::perfect_hash_default(feature)
                )
            ),
            HashType::TertiaryInteraction => GeometricHash::TertiaryInteraction(
                super::tertiary_interaction::HashValue(
                    super::tertiary_interaction::HashValue::perfect_hash_default(feature)
                )
            ),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }

    pub fn perfect_hash(
        feature: &Vec<f32>, hash_type: HashType, nbin_dist: usize, nbin_angle: usize
    ) -> Self {
        match hash_type {
            HashType::PDBMotif => GeometricHash::PDBMotif(
                super::pdb_motif::HashValue(
                    super::pdb_motif::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::PDBMotifSinCos => GeometricHash::PDBMotifSinCos(
                super::pdb_motif_sincos::HashValue(
                    super::pdb_motif_sincos::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::TrRosetta => GeometricHash::TrRosetta(
                super::trrosetta::HashValue(
                    super::trrosetta::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::PDBTrRosetta => GeometricHash::PDBTrRosetta(
                super::pdb_tr::HashValue(
                    super::pdb_tr::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue(
                    super::ppf::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            HashType::TertiaryInteraction => GeometricHash::TertiaryInteraction(
                super::tertiary_interaction::HashValue(
                    super::tertiary_interaction::HashValue::perfect_hash(feature, nbin_dist, nbin_angle)
                )
            ),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }

    pub fn reverse_hash_default(&self, output: &mut Vec<f32>) {
        match self {
            GeometricHash::PDBMotif(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::PDBMotifSinCos(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::TrRosetta(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::PDBTrRosetta(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            }
            GeometricHash::PointPairFeature(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::TertiaryInteraction(hash) => {
                let reversed = hash.reverse_hash_default();
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            // append new hash type here
            // _ => panic!("Invalid hash type"),
        }
    }


    pub fn reverse_hash(&self, nbin_dist: usize, nbin_angle: usize, output: &mut Vec<f32>) {
        match self {
            GeometricHash::PDBMotif(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::PDBMotifSinCos(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::TrRosetta(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::PDBTrRosetta(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::PointPairFeature(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            GeometricHash::TertiaryInteraction(hash) => {
                let reversed = hash.reverse_hash(nbin_dist, nbin_angle);
                for i in 0..reversed.len() {
                    output[i] = reversed[i];
                }
            },
            // append new hash type here
            // _ => panic!("Invalid hash type"),
        }
    }


    pub fn hash_type(&self) -> HashType {
        match self {
            GeometricHash::PDBMotif(hash) => hash.hash_type(),
            GeometricHash::PDBMotifSinCos(hash) => hash.hash_type(),
            GeometricHash::TrRosetta(hash) => hash.hash_type(),
            GeometricHash::PointPairFeature(hash) => hash.hash_type(),
            GeometricHash::PDBTrRosetta(hash) => hash.hash_type(),
            GeometricHash::TertiaryInteraction(hash) => hash.hash_type(),
            // append new hash type here
            // _ => panic!("Invalid hash type"),
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
            HashType::TrRosetta => GeometricHash::TrRosetta(
                super::trrosetta::HashValue::from_u32(hashvalue)
            ),
            HashType::PDBTrRosetta => GeometricHash::PDBTrRosetta(
                super::pdb_tr::HashValue::from_u32(hashvalue)
            ),
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue::from_u32(hashvalue)
            ),
            HashType::TertiaryInteraction => GeometricHash::TertiaryInteraction(
                super::tertiary_interaction::HashValue::from_u32(hashvalue)
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
            HashType::PDBTrRosetta => GeometricHash::PDBTrRosetta(
                super::pdb_tr::HashValue::from_u64(hashvalue)
            ),
            HashType::PointPairFeature => GeometricHash::PointPairFeature(
                super::ppf::HashValue::from_u64(hashvalue)
            ),
            HashType::TertiaryInteraction => GeometricHash::TertiaryInteraction(
                super::tertiary_interaction::HashValue::from_u64(hashvalue)
            ),
            // append new hash type here
            _ => panic!("Invalid hash type"),
        }
    }
    
    pub fn as_u32(&self) -> u32 {
        match self {
            GeometricHash::PDBMotif(hash) => hash.as_u32(),
            GeometricHash::PDBMotifSinCos(hash) => hash.as_u32(),
            GeometricHash::TrRosetta(hash) => hash.as_u32(),
            GeometricHash::PDBTrRosetta(hash) => hash.as_u32(),
            GeometricHash::PointPairFeature(hash) => hash.as_u32(),
            GeometricHash::TertiaryInteraction(hash) => hash.as_u32(),
            // append new hash type here
        }
    }
    pub fn as_u64(&self) -> u64 {
        match self {
            GeometricHash::PDBMotif(hash) => hash.as_u64(),
            GeometricHash::PDBMotifSinCos(hash) => hash.as_u64(),
            GeometricHash::TrRosetta(hash) => hash.as_u64(),
            GeometricHash::PDBTrRosetta(hash) => hash.as_u64(),
            GeometricHash::PointPairFeature(hash) => hash.as_u64(),
            GeometricHash::TertiaryInteraction(hash) => hash.as_u64(),
            // append new hash type here
        }
    }
    
    pub fn is_symmetric(&self) -> bool {
        match self {
            GeometricHash::PDBMotif(hash) => hash.is_symmetric(),
            GeometricHash::PDBMotifSinCos(hash) => hash.is_symmetric(),
            GeometricHash::TrRosetta(hash) => hash.is_symmetric(),
            GeometricHash::PDBTrRosetta(hash) => hash.is_symmetric(),
            GeometricHash::PointPairFeature(hash) => hash.is_symmetric(),
            GeometricHash::TertiaryInteraction(hash) => hash.is_symmetric(),
            // append new hash type here
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
    pub fn downcast_default_32bit(&self) -> super::trrosetta::HashValue {
        match self {
            GeometricHash::TrRosetta(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_point_pair_feature(&self) -> super::ppf::HashValue {
        match self {
            GeometricHash::PointPairFeature(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_pdb_tr(&self) -> super::pdb_tr::HashValue {
        match self {
            GeometricHash::PDBTrRosetta(hash) => hash.clone(),
            _ => panic!("Invalid hash type"),
        }
    }
    pub fn downcast_tertiary_interaction(&self) -> super::tertiary_interaction::HashValue {
        match self {
            GeometricHash::TertiaryInteraction(hash) => hash.clone(),
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
            GeometricHash::PDBTrRosetta(hash) => {
                write!(f, "PDBTrRosetta({:?})", hash)
            },
            GeometricHash::PointPairFeature(hash) => {
                write!(f, "PointPairFeature({:?})", hash)
            },
            GeometricHash::TertiaryInteraction(hash) => {
                write!(f, "TertiaryInteraction({:?})", hash)
            },
            // append new hash type here
            // _ => panic!("Invalid hash type"),
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
            GeometricHash::PDBTrRosetta(hash) => {
                write!(f, "PDBTrRosetta\t{:?}", hash)
            },
            GeometricHash::PointPairFeature(hash) => {
                write!(f, "PointPairFeature\t{:?}", hash)
            },
            GeometricHash::TertiaryInteraction(hash) => {
                write!(f, "TertiaryInteraction\t{:?}", hash)
            },
            // append new hash type here
            // _ => panic!("Invalid hash type"),
        }
    }
}