//! # About project
//!
//! Motifsearch is a tool for finding discontinuous motifs in protein structures.

//! # Things to keep in mind when developing
//! * No dependencies between huge modules
//!   * geometry ←→ structure ←→ index
//! * Make errors visible
//!   * No unwrap, use Option/Result/expect/match
//! * Write tests
//!   * Unit tests at source files
//!   * Integration tests at `tests/`

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
