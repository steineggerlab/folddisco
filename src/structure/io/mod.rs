//!
pub mod pdb;
pub mod parser;

#[derive(Debug)]
pub enum StructureFileFormat {
    PDB,
    CIF,
    FCZ,
    MMTF,
    Unknown,
}

// pub trait Parse {

// }
