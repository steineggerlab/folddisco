//!

#[derive(Debug)]
pub enum StructureFileFormat {
    PDB,
    CIF,
    FCZ,
    MMTF,
    Unknown,
}

pub mod pdb;