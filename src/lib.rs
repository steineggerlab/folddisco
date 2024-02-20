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

use std::hash::Hash;
use std::fmt::Debug;
pub mod cli;
pub mod controller;
pub mod geometry;
pub mod index;
pub mod structure;
pub mod utils;

/* re-export: pub use */
pub use structure::io::pdb::Reader as PDBReader;

// Declare a new trait that supports required traits
pub trait HashableSync: Clone + Copy + Hash + Sync + Send + Eq + PartialEq + Ord + Debug + 'static {}

impl HashableSync for usize {}
impl HashableSync for u64 {}
impl HashableSync for u32 {}
impl HashableSync for u16 {}
impl HashableSync for u8 {}
impl HashableSync for isize {}
impl HashableSync for i64 {}
impl HashableSync for i32 {}
impl HashableSync for i16 {}
impl HashableSync for i8 {}
impl HashableSync for char {}


pub mod prelude {
    pub use crate::PDBReader;
    pub use crate::HashableSync;
    pub use crate::measure_time;

    pub use crate::controller::FoldDisco;
    
    pub use crate::geometry::core::{GeometricHash, HashType};
    
    pub use crate::index::builder::IndexBuilder;
    pub use crate::index::query_multiple_with_neighbors;
    pub use crate::index::{IndexTablePrinter, query_single, query_multiple};
    pub use crate::index::alloc::IndexAllocator;
    pub use crate::index::index_table::IndexTable;

    pub use crate::utils::loader::{load_path, get_all_combination};
    pub use crate::utils::benchmark::{Metrics, calculate_metrics, compare_target_answer};
    pub use crate::utils::log::{INFO, FAIL, WARN, DONE, log_msg, print_log_msg};
}