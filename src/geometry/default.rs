// Default hasher from PDB's implementation
// 5bit aa1, 5bit aa2, 

use std::fmt;
use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u32);

