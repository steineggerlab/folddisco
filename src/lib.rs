//! # Folddisco
//!
//! Folddisco is a tool for finding discontinuous motifs in protein structures.

use std::hash::{Hash, Hasher};
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
pub trait HashableSync: Clone + Copy + Hash + Sync + Send + Eq + PartialEq + Ord + Debug + 'static {
    fn hash_u32(&self) -> u32 {
        use rustc_hash::FxHasher;
        let mut hasher = FxHasher::default();
        self.hash(&mut hasher);
        hasher.finish() as u32
    }
}

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
    pub use crate::controller::io::{read_offset_map, save_offset_map, write_usize_vector};
    pub use crate::controller::query::{make_query_map, parse_query_string};

    pub use crate::geometry::core::{GeometricHash, HashType};
    
    pub use crate::index::lookup::{save_lookup_to_file, load_lookup_from_file};
    pub use crate::index::alloc::IndexBuilder;
    pub use crate::index::alloc::convert_sorted_pairs_to_offset_and_values_vec;
    
    pub use crate::utils::loader::load_path;
    pub use crate::utils::benchmark::{Metrics, compare_target_answer_set};
    pub use crate::utils::log::{INFO, FAIL, WARN, DONE, log_msg, print_log_msg};
}