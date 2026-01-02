// A utility module defining traits used across the project.

use std::hash::{Hash, Hasher};
use std::fmt::Debug;
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

