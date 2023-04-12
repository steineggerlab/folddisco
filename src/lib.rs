//! # About project
//!
//! Motifsearch is a tool for finding discontinuous motifs in protein structures.

pub mod cli;
pub mod controller;
pub mod geometry;
pub mod index;
pub mod structure;
pub mod test;
pub mod utils;

pub struct MotifSearch {
    pub pdb_files: Vec<String>,
    // pub pdb_table: HashMap<index, PDB>,
    pub controller: controller::Controller,
}

pub fn run() {}
