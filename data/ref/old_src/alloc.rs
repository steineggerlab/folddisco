use std::{sync::{Arc, atomic::{AtomicUsize, Ordering}, Mutex}, thread, io::Write};
use std::fs::OpenOptions;
use std::io::Error;
use memmap2::{MmapOptions, Mmap};

pub struct IndexAllocator {
    pub num_threads: usize,
    pub data_size: AtomicUsize,
    pub allocation: Arc<Vec<AtomicUsize>>,
    pub cursor: Arc<AtomicUsize>,
}

impl IndexAllocator {
    pub fn new(num_threads: usize, data_size: usize) -> Self {
        let zero_vec = Arc::new(
            (0..data_size).map(|_| AtomicUsize::new(0)).collect::<Vec<_>>()
        );

        IndexAllocator {
            num_threads,
            data_size: AtomicUsize::new(data_size),
            allocation: zero_vec,
            cursor: Arc::new(AtomicUsize::new(0)),
        }
    }

    // This method is expected to be called from multiple threads
    pub fn fill_and_update(&self, ext_data: Arc<Vec<AtomicUsize>>) {
        let len = ext_data.len();
        let start = self.cursor.fetch_add(len, Ordering::SeqCst);
        let size = self.data_size.load(Ordering::Relaxed);
        for (i, value) in ext_data.iter().enumerate() {
            self.allocation[start + i].store(
                value.load(Ordering::Relaxed), Ordering::Relaxed
            );
        }
    }
    
    // This method is expected to be called from multiple threads
    pub fn fill_usize_vec(&self, ext_data: Vec<usize>) {
        let len = ext_data.len();
        let start = self.cursor.fetch_add(len, Ordering::SeqCst);
        let size = self.data_size.load(Ordering::Relaxed);
        for (i, value) in ext_data.iter().enumerate() {
            self.allocation[start + i].store(*value, Ordering::Relaxed);
        }
    }
    
    // This method uses unsafe code
    pub fn save_to_file(&self, filename: &str) -> Result<(), Error> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(filename)?;

        file.set_len((self.data_size.load(Ordering::Relaxed) * std::mem::size_of::<usize>()) as u64)?;

        let vec_to_write = unsafe {
            std::slice::from_raw_parts(
                self.allocation.as_ptr() as *const u8,
                self.allocation.len() * std::mem::size_of::<usize>()
            )
        };
        
        let mut writer = std::io::BufWriter::new(&file);
        writer.write_all(vec_to_write)?;

        Ok(())
    }
}


// TODO: NEED TO CHECK IF THIS IMPLEMENTATION IS RIGHT
pub fn load_from_file(filename: &str) -> Result<(Mmap, &'static [usize]), Error> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(filename)?;

    let mmap = unsafe { MmapOptions::new().map(&file)? };

    let vec = unsafe {
        std::slice::from_raw_parts(
            mmap.as_ptr() as *const usize,
            mmap.len() / std::mem::size_of::<usize>()
        )
    };
    Ok((mmap, vec))
}

pub fn resize_allocation(vec: &mut [AtomicUsize], new_size: usize) {
    // Re allocate
    
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Read;
    use std::time::{Instant, Duration};
    use crate::measure_time;

    #[test]
    fn test_fill_and_update_multithread() {
        // Create an instance of IndexAllocator and allocate memory
        let allocator = Arc::new(IndexAllocator::new(8, 10000000));

        // Create a vector of vectors to be filled into allocation
        let values = Arc::new((0..1000).map(|x|
            Arc::new((0..10000).map(|y|
                AtomicUsize::new(x + y)
            ).collect::<Vec<_>>())
        ).collect::<Vec<_>>());
        
        // Spawn multiple threads and call the fill_and_update method in each thread
        let mut handles = vec![];
        let num_threads = 4;
        // Use 4 threads to iterate through the values vector
        for i in 0..num_threads {
            let allocator = Arc::clone(&allocator);
            let values = Arc::clone(&values);
            let handle = thread::spawn(move || {
                let offset = i * (values.len() / num_threads); // Calculate the offset for this thread
                // If the data size is not evenly divisible by the number of threads then
                // the last thread will have to do more work
                let data_size = if i == num_threads - 1 {
                    values.len() - offset
                } else {
                    values.len() / num_threads
                };
                for j in 0..data_size {
                    // Save value to the data vector at the correct offset
                    // Value should be the offset + the current value at the offset
                    let index = offset + j;
                    allocator.fill_and_update(Arc::clone(&values[index]));
                }
            });
            handles.push(handle);
        }
        // Join all threads to ensure all updates are completed
        for handle in handles {
            handle.join().unwrap();
        }

        // Write the allocation to a file
        allocator.save_to_file("output.bin").unwrap();
        
    }
    
    // #[test]
    // pub fn test_resize() {
    //     let mut allocator = IndexAllocator::new(8, 10000000);
    //     resize_allocation(&mut allocator, 20000000);
    //     assert_eq!(allocator.data_size.load(Ordering::Relaxed), 20000000);
    //     println!("Data size: {}", allocator.data_size.load(Ordering::Relaxed));
    //     println!("Allocation size: {}", allocator.allocation.len());
    // }
    
}

