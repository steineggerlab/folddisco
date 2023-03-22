
use std::hash::Hash;
use std::cmp::{Ord, Eq};
use crate::index::IndexTable;

pub struct IndexBuilder{
}

pub trait Indexable: Hash + Ord + Eq {}

impl IndexBuilder {
    pub fn concat<T: Indexable>(&self, hash_vec: Vec<T>) -> IndexTable {
        todo!()
    }
    pub fn sort<T: Indexable>(&self, index_table: &mut Vec<T>) {
        todo!()
    }
    pub fn diff<T: Indexable>(&self, index_table: &mut Vec<T>) {
        todo!()
    }
}