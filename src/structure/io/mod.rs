//!
pub mod parser;
pub mod pdb;

#[cfg(feature = "foldcomp")]
pub mod fcz;

#[derive(Debug)]
pub enum StructureFileFormat {
    PDB,
    PDBGZ,
    CIF,
    FCZ,
    MMTF,
    Unknown,
}

// pub trait Parse {

// }
