
pub struct IndexAllocator {
    pub lookup_table: Vec<usize>,
    pub starting_index: usize,
    pub size: usize,
    pub offset_table: Vec<usize>,
    pub big_allocation: Box<[usize]>,
}

impl IndexAllocator {
    pub fn allocate(&self, size: usize) -> Result<usize, &'static str> {
        // Ok(pointer)
        // Err("Not enough space")
        todo!()
    }
    pub fn fill(&self, pointer: usize, size: usize, value: u8) {
        todo!()
    }
    pub fn build_offset_table(&self) {
        todo!()
    }
}