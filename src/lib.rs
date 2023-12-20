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

pub mod prelude {
    pub use crate::controller::Controller;
    pub use crate::geometry::trrosetta_subfamily::{HashCollection, HashValue};
    pub use crate::index::builder::IndexBuilder;
    pub use crate::index::query_multiple_with_neighbors;
    pub use crate::index::{IndexTablePrinter, query_single, query_multiple};
    pub use crate::index::alloc::IndexAllocator;
    pub use crate::PDBReader;
    pub use crate::index::index_table::IndexTable;
    pub use crate::index::io::save_offset_map;
    pub use crate::utils::loader::{load_path, get_all_combination};
    pub use crate::utils::benchmark::{Metrics, calculate_metrics, compare_target_answer};
    pub use crate::utils::log::{INFO, FAIL, WARN, DONE, log_msg, print_log_msg};
    pub use crate::measure_time;
}