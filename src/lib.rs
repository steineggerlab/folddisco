//! # About project
//!
//! Motifsearch is a tool for finding discontinuous motifs in protein structures.

pub mod cli;
pub mod controller;
pub mod geometry;
pub mod index;
pub mod structure;
pub mod utils;

/* re-export: pub use */
pub use structure::io::pdb::Reader as PDBReader;

pub struct MotifSearch {
    pub pdb_files: Vec<String>,
    // pub pdb_table: HashMap<index, PDB>,
    pub controller: controller::Controller,
}

pub fn run() {}
