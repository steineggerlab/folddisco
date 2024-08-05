//!

use std::fmt;
pub mod parser;
pub mod pdb;

#[cfg(feature = "foldcomp")]
pub mod fcz;

#[derive(Debug, Clone, PartialEq)]
pub enum StructureFileFormat {
    PDB,
    PDBGZ,
    CIF, // IMPORTANT: TODO: Implement CIF parser
    FCZ,
    FCZDB,
    MMTF, // TODO: Implement MMTF parser
    Unknown,
}

impl StructureFileFormat {
    pub fn to_string(&self) -> String {
        match *self {
            StructureFileFormat::PDB => "PDB".to_string(),
            StructureFileFormat::PDBGZ => "PDBGZ".to_string(),
            StructureFileFormat::CIF => "CIF".to_string(),
            StructureFileFormat::FCZ => "FCZ".to_string(),
            StructureFileFormat::FCZDB => "FCZDB".to_string(),
            StructureFileFormat::MMTF => "MMTF".to_string(),
            StructureFileFormat::Unknown => "Unknown".to_string(),
        }
    }
    
    pub fn get_with_string(s: &str) -> StructureFileFormat {
        match s {
            "0" | "PDB" | "pdb" => StructureFileFormat::PDB,
            "1" | "PDBGZ" | "pdbgz" => StructureFileFormat::PDBGZ,
            "2" | "CIF" | "cif" => StructureFileFormat::CIF,
            "3" | "FCZ" | "fcz" => StructureFileFormat::FCZ,
            "4" | "FCZDB" | "fczdb" => StructureFileFormat::FCZDB,
            "5" | "MMTF" | "mmtf" => StructureFileFormat::MMTF,
            _ => StructureFileFormat::Unknown,
        }
    }
}


impl fmt::Display for StructureFileFormat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            StructureFileFormat::PDB => write!(f, "PDB"),
            StructureFileFormat::PDBGZ => write!(f, "PDBGZ"),
            StructureFileFormat::CIF => write!(f, "CIF"),
            StructureFileFormat::FCZ => write!(f, "FCZ"),
            StructureFileFormat::FCZDB => write!(f, "FCZDB"),
            StructureFileFormat::MMTF => write!(f, "MMTF"),
            StructureFileFormat::Unknown => write!(f, "Unknown"),
        }
    }
}

impl From<&str> for StructureFileFormat {
    fn from(s: &str) -> Self {
        match s {
            "PDB" => StructureFileFormat::PDB,
            "PDBGZ" => StructureFileFormat::PDBGZ,
            "CIF" => StructureFileFormat::CIF,
            "FCZ" => StructureFileFormat::FCZ,
            "FCZDB" => StructureFileFormat::FCZDB,
            "MMTF" => StructureFileFormat::MMTF,
            _ => StructureFileFormat::Unknown,
        }
    }
}

impl From<String> for StructureFileFormat {
    fn from(s: String) -> Self {
        match s.as_str() {
            "PDB" => StructureFileFormat::PDB,
            "PDBGZ" => StructureFileFormat::PDBGZ,
            "CIF" => StructureFileFormat::CIF,
            "FCZ" => StructureFileFormat::FCZ,
            "FCZDB" => StructureFileFormat::FCZDB,
            "MMTF" => StructureFileFormat::MMTF,
            _ => StructureFileFormat::Unknown,
        }
    }
}

impl From<&String> for StructureFileFormat {
    fn from(s: &String) -> Self {
        match s.as_str() {
            "PDB" => StructureFileFormat::PDB,
            "PDBGZ" => StructureFileFormat::PDBGZ,
            "CIF" => StructureFileFormat::CIF,
            "FCZ" => StructureFileFormat::FCZ,
            "FCZDB" => StructureFileFormat::FCZDB,
            "MMTF" => StructureFileFormat::MMTF,
            _ => StructureFileFormat::Unknown,
        }
    }
}
