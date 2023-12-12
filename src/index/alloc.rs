use std::{sync::{Arc, atomic::{AtomicUsize, Ordering}, Mutex}, io::Write, fs::OpenOptions};
use std::thread;
use std::time::Instant;
use std::env;
use std::process;
use memmap2::{MmapOptions, MmapMut};
use crate::index::IndexTable;
use std::collections::HashMap;

pub struct IndexAllocator {
    pub size: usize,
    pub offset_table: HashMap<usize, usize>,
    pub big_allocation: Arc<Vec<usize>>,
}

impl IndexAllocator {
    pub fn new(size: usize) -> Self {
        Self {
            size,
            offset_table: HashMap::new(),
            big_allocation: Arc::new(vec![0; size]),
        }
    }

    pub fn allocate(&mut self, size: usize) -> Result<_,()> {
        // Allocate a big chunk of memory
        self.big_allocation = Arc::new(vec![0; size]);
        Ok(())
    }

    pub fn fill(&mut self, values: Vec<usize>, offset: usize, num_threads: usize) {
        // Fill values into big_allocation
        let mut handles = vec![];
        let data = Arc::new(values);
        for i in 0..num_threads {
            let data = Arc::clone(&data);
            let handle = thread::spawn(move || {
                let offset = i * (data.len() / num_threads); // Calculate the offset for this thread
                // If the data size is not evenly divisible by the number of threads then
                // the last thread will have to do more work
                let data_size = if i == num_threads - 1 {
                    data.len() - offset
                } else {
                    data.len() / num_threads
                };
                for i in 0..data_size {
                    // Save value to the data vector at the correct offset
                    // Value should be the offset + the current value at the offset
                    let index = offset + i;
                    self.big_allocation[index] = data[index];
                }
            });
            handles.push(handle);
        }
    }
    pub fn build_offset_table(&self) {
        todo!()
    }
    
    pub fn save(&self, path: &str) {
        // Unsafe.
        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(path)
            .unwrap();

        file.set_len((self.size * std::mem::size_of::<usize>()) as u64).unwrap();

        let mut mmap = unsafe { MmapOptions::new().map_mut(&file).unwrap() };
        for (i, value) in self.big_allocation.iter().enumerate() {
            let bytes = unsafe { std::slice::from_raw_parts((value as usize).to_ne_bytes().as_ptr(), std::mem::size_of::<usize>()) };
            mmap[i * std::mem::size_of::<usize>()..(i+1) * std::mem::size_of::<usize>()].copy_from_slice(bytes);
        }

        mmap.flush().expect("Failed to write");
    }
}
