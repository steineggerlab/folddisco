//!
pub mod parser;
pub mod pdb;

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
