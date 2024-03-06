//!
pub mod parser;
pub mod pdb;

#[derive(Debug)]
pub enum StructureFileFormat {
    PDB,
    PDB_GZ,
    CIF,
    FCZ,
    MMTF,
    Unknown,
}

// pub trait Parse {

// }
