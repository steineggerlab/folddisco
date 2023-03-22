use std::collections::HashMap;
use crate::index::IndexTable;

pub struct IndexAllocator {
    pub lookup_table: Vec<usize>,
    pub starting_index: usize,
    pub size: usize,
    pub offset_table: HashMap<usize, usize>,
    pub big_allocation: Box<[usize]>,
}

impl IndexAllocator {
    pub fn allocate(&self, size: usize) -> Result<usize, &str> {
        // Ok(pointer)
        // Err("Not enough space")
        todo!()
    }
    pub fn fill(&self, index_table: &IndexTable) {
        todo!()
    }
    pub fn build_offset_table(&self) {
        todo!()
    }
}