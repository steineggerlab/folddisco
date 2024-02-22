// File: nested_vec.rs
// Created: 2023-12-28 16:44:31
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

use rayon::iter::ParallelIterator;
// external crates
use rustc_hash::{FxHashMap, FxHashSet};
use dashmap::DashMap;

use std::collections::BTreeMap;
// 
use std::sync::{Arc, Mutex};
use std::cell::UnsafeCell;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::hash::Hash;
// Measure time
use crate::measure_time;
use crate::structure::atom::Atom;
use crate::HashableSync;

const DEFAULT_NUM_THREADS: usize = 4;
const DEFAULT_ALLOCATION_SIZE: usize = 1024;

pub struct HugeAllocation<T: HashableSync> {
    pub allocation: UnsafeCell<Vec<T>>,
}
unsafe impl<T: HashableSync> Sync for HugeAllocation<T> {}

impl<T: HashableSync> HugeAllocation<T> {
    pub fn new(size: usize) -> HugeAllocation<T> {
        HugeAllocation {
            allocation: UnsafeCell::new(vec![unsafe { std::mem::zeroed() }; size]),
        }
    }
}

pub struct IndexBuilder<K: HashableSync, V: HashableSync> {
    pub offset: DashMap<V, (AtomicUsize, AtomicUsize)>, // Value, (offset in allocation, size)
    pub allocation: Arc<HugeAllocation<K>>,
    pub ids: Arc<Vec<K>>,
    pub data: Arc<Vec<Vec<V>>>,
    // Alternatively, sort with dashmap and save it. 
    // This will be tested with allocating approach
    pub data_dashmap: Arc<DashMap<V, Vec<K>>>,
    // Configurations
    pub num_threads: usize,
    pub allocation_size: usize,
    pub offset_path: String,
    pub data_path: String,
}

// IndexBuilder is accessible from multiple threads
unsafe impl<K: HashableSync, V: HashableSync> Sync for IndexBuilder<K, V> {}

impl<K: HashableSync, V: HashableSync> IndexBuilder<K, V> {
    pub fn empty() -> IndexBuilder<K, V> {
        IndexBuilder {
            offset: DashMap::new(),
            allocation: Arc::new(HugeAllocation::new(0)),
            ids: Arc::new(Vec::new()),
            data: Arc::new(Vec::new()),
            data_dashmap: Arc::new(DashMap::new()),
            num_threads: 0,
            allocation_size: 0,
            offset_path: String::from(""),
            data_path: String::from(""),
        }
    }
    
    // Constructor
    pub fn new(
        ids: &Vec<K>, data: &Vec<Vec<V>>,
        num_threads: usize, allocation_size: usize,
        offset_path: String, data_path: String,
    ) -> IndexBuilder<K, V> {
        // Before creating check if given ids and data have same length
        assert_eq!(ids.len(), data.len());
        // If num_threads is 0, set it to default. one liner
        let num_threads = if num_threads == 0 { DEFAULT_NUM_THREADS } else { num_threads };
        // If allocation_size is 0, set it to default
        let allocation_size = if allocation_size == 0 { DEFAULT_ALLOCATION_SIZE } else { allocation_size };
        // Create IndexBuilder
        IndexBuilder {
            offset: DashMap::new(),
            allocation: Arc::new(HugeAllocation::new(allocation_size)),
            ids: Arc::new(ids.to_owned()),
            data: Arc::new(data.to_owned()),
            data_dashmap: Arc::new(DashMap::new()),
            num_threads,
            allocation_size,
            offset_path,
            data_path,
        }
    }
    
    // Estimate the size of allocation
    pub fn estimate_size(&self) -> usize {
        let mut total_size = 0;
        for i in 0..self.data.len() {
            total_size += self.data[i].len();
        }
        total_size
    }

    pub fn fill_offset_map(&mut self) {
        let mut observed: FxHashMap<V, usize> = FxHashMap::default();
        for i in 0..self.data.len() {
            let data_inner = &self.data[i];
            for j in 0..data_inner.len() {
                let val = data_inner[j].clone();
                // If val is not in observed, add it
                if !observed.contains_key(&val) {
                    observed.insert(val.clone(), 1);
                } else {
                    // If val is in observed, increment the value
                    let val_count = observed.get_mut(&val).unwrap();
                    *val_count += 1;
                }
            }
        }
        // Sort observed by V
        let mut observed_vec = observed.into_iter().collect::<Vec<(V, usize)>>();
        observed_vec.sort_by(|a, b| a.0.cmp(&b.0));
        // Fill offset map
        let mut offset = 0;
        for i in 0..observed_vec.len() {
            let val = observed_vec[i].0.clone();
            let val_count = observed_vec[i].1;
            self.offset.insert(val, (AtomicUsize::new(offset), AtomicUsize::new(val_count)));
            offset += val_count;
        }
    }
    
    pub fn estimate_size_multi(&self) -> usize {
        let total_size = Arc::new(AtomicUsize::new(0));
        let ext_data_index = Arc::new(AtomicUsize::new(0));
        let mut handles = vec![];
        for i in 0..self.num_threads {
            let ext_data = self.data.clone();
            let total_size = total_size.clone();
            let ext_data_index = ext_data_index.clone();
            let handle = thread::spawn(move || {
                while ext_data_index.load(Ordering::Relaxed) < ext_data.len() {
                    let ext_data_index = ext_data_index.fetch_add(1, Ordering::Relaxed);
                    if ext_data_index >= ext_data.len() {
                        break;
                    }
                    let ext_data_inner = &ext_data[ext_data_index];
                    total_size.fetch_add(ext_data_inner.len(), Ordering::Relaxed);
                }
            });
            handles.push(handle);
        }
        for handle in handles {
            handle.join().unwrap();
        }
        total_size.load(Ordering::Relaxed)
    }
    
    
    pub fn estimate_and_set_size(&mut self) {
        let total_size = self.estimate_size();
        self.allocation_size = total_size;
    }
    
    
    pub fn fill_with_dashmap(&mut self) {
        let mut handles = vec![];
        let data_index = Arc::new(AtomicUsize::new(0));
        // 
        for _ in 0..self.num_threads {
            // Clone Arcs to move into threads
            let data = self.data.clone();   
            let ids = self.ids.clone();
            let data_index = data_index.clone();
            let data_dashmap = self.data_dashmap.clone();
            let handle = thread::spawn(move || {
                while data_index.load(Ordering::Relaxed) < data.len() {
                    let data_index = data_index.fetch_add(1, Ordering::Relaxed);
                    if data_index >= data.len() {
                        break;
                    }
                    let data_inner = &data[data_index];
                    let curr_id = ids[data_index];
                    for j in 0..data_inner.len() {
                        // If data_dashmap contains data_inner[j], append ids[data_index] to the value
                        if data_dashmap.contains_key(&data_inner[j]) {
                            let mut value = data_dashmap.get_mut(&data_inner[j]).unwrap();
                            if !value.contains(&curr_id) {
                                value.push(curr_id);
                            }
                        } else {
                            // If data_dashmap does not contain data_inner[j], insert ids[data_index] to the value
                            data_dashmap.insert(data_inner[j].clone(), vec![curr_id]);
                        }
                    }
                    
                }
            });
            handles.push(handle);
        }
        for handle in handles {
            handle.join().unwrap();
        }
    }

    pub fn fill_and_return_dashmap(&mut self) -> DashMap<V, Vec<K>>  {
        let output = Arc::new(DashMap::<V, Vec<K>>::new());
        let mut handles = vec![];
        let data_index = Arc::new(AtomicUsize::new(0));
        for _ in 0..self.num_threads {
            // Clone Arcs to move into threads
            let data = self.data.clone();   
            let ids = self.ids.clone();
            let data_index = data_index.clone();
            let output_ref_inner = output.clone();
            let handle = thread::spawn(move || {
                while data_index.load(Ordering::Acquire) < data.len() {
                    if data_index.load(Ordering::Acquire) >= data.len() {
                        break;
                    }
                    let i = data_index.fetch_add(1, Ordering::Release);
                    let data_inner = &data[i];
                    let curr_id = ids[i];
                    data_inner.iter().for_each(|curr_hash| {
                        let value = output_ref_inner.get_mut(curr_hash);
                        match value {
                            Some(mut v) => { v.push(curr_id); },
                            None => { output_ref_inner.insert(curr_hash.clone(), vec![curr_id]); }
                        }
                    });
                }
            });
            handles.push(handle);
        }
        for handle in handles {
            handle.join().unwrap();
        }
        Arc::try_unwrap(output).unwrap()
    }
    
    
    // Allocate memory with allocation size
    pub fn allocate(&mut self) {
        self.allocation = Arc::new(HugeAllocation::new(self.allocation_size));
    }
    
    pub fn fill_multiple(&mut self) {
        let mut handles = vec![];
        let data_index = Arc::new(AtomicUsize::new(0));
        // 
        for _ in 0..self.num_threads {
            let data = self.data.clone();   
            let ids = self.ids.clone();
            let data_index = data_index.clone();
            let allocation = self.allocation.clone();
            let handle = thread::spawn(move || {
                while data_index.load(Ordering::Relaxed) < data.len() {
                    let data_index = data_index.fetch_add(1, Ordering::Relaxed);
                    if data_index >= data.len() {
                        break;
                    }
                    let data_inner = &data[data_index];
                    for j in 0..data_inner.len() {
                        // Get offset from offsets map for value of j
                        let val = data_inner[j];
                    }
                    // Update offset map
                    // self.offset.insert(data_inner[j].clone(), offset);
                }
            });
            handles.push(handle);
        }
        for handle in handles {
            handle.join().unwrap();
        }
        
    }

    pub fn convert_hashmap_to_offset_and_values(
        &self,
        orig_map: DashMap<V, Vec<K>>
    ) -> (DashMap<V, (usize, usize)>, Vec<K>) {
        // OffsetMap - key: hash, value: (offset, length)
        // Vec - all values concatenated
        println!("Orig map length: {}", orig_map.len());
        let mut offset_map = DashMap::new();
        let mut vec: Arc<Mutex<Vec<K>>> = Arc::new(Mutex::new(Vec::with_capacity(orig_map.len())));
        let mut offset = AtomicUsize::new(0);
        
        // Parallel 
        orig_map.par_iter_mut().for_each(|entry| {
            
            let hash = entry.key().to_owned();
            let ids = entry.value();
            offset_map.insert(
                hash,
                (offset.fetch_add(ids.len(), Ordering::Relaxed), ids.len())
            );
            let vec_inner = vec.clone();
            vec_inner.lock().unwrap().extend_from_slice(ids);
        });
        
        // orig_map.iter().for_each(|entry| {
        //     let hash = entry.key().to_owned();
        //     let ids = entry.value();
        //     offset_map.insert(
        //         hash,
        //         (offset.fetch_add(ids.len(), Ordering::Relaxed), ids.len())
        //     );
        //     vec.extend_from_slice(ids);
        // });
        (offset_map, Arc::try_unwrap(vec).unwrap().into_inner().unwrap())
    }


    // fn fill_inner
    // Write dashmap to file
    // pub fn 

    
    
    
    
}

#[cfg(test)]
mod tests {
    use super::*;
    
    fn create_test_data(num_key: usize, num_value: usize) -> Vec<Vec<usize>> {
        let mut data = Vec::new();
        for i in 0..num_key {
            let mut inner = Vec::new();
            for j in 0..num_value {
                inner.push(i * num_value + j);
            }
            data.push(inner);
        }
        println!("Created test data with {} keys and {} values", num_key, num_value);
        println!("data.len(): {}", data.len());
        data
    }
    
    #[test]
    fn test_estimate_size() {
        let data = create_test_data(40000, 10000);
        let ids = (0..40000).collect::<Vec<usize>>();
        let index_builder = IndexBuilder::new(
            &ids, &data, 6, 1000, String::from(""), String::from("")
        );
        measure_time!(assert_eq!(index_builder.estimate_size(), 40000 * 10000));
    }
    
    #[test]
    fn test_fill_dashmap() {
        let data = create_test_data(40000, 10000);
        let ids = (0..40000).collect::<Vec<usize>>();
        let mut index_builder = IndexBuilder::new(
            &ids, &data, 6, 1000, String::from(""), String::from("")
        );
        measure_time!(index_builder.fill_with_dashmap());
        assert_eq!(index_builder.data_dashmap.len(), 40000);
    }

    #[test]
    fn test_fill_offset_map() {
        // Currently, TOO SLOW
        let data = create_test_data(40000, 10000);
        let ids = (0..40000).collect::<Vec<usize>>();
        let mut index_builder = IndexBuilder::new(
            &ids, &data, 6, 1000, String::from(""), String::from("")
        );
        index_builder.fill_offset_map();
        assert_eq!(index_builder.offset.len(), 40000);
        for i in 0..10 {
            println!("{:?}", index_builder.offset.get(&i));
        }
    }
    
}










// pub struct HugeAllocation {
//     pub allocation: UnsafeCell<Vec<usize>>,
// }

// unsafe impl Sync for HugeAllocation {}

// pub fn run(num_threads: usize, ext_data: Arc<Vec<Vec<usize>>>) -> Arc<HugeAllocation> {
//     println!("Creating {} threads", num_threads);
//     // Iterate through ext_data and get the total size
//     let total_size = Arc::new(AtomicUsize::new(0));
//     let ext_data_index = Arc::new(AtomicUsize::new(0));
//     // Spawn threads to find out the size to allocate
//     let start = Instant::now();
//     let mut handles = vec![];
//     for i in 0..num_threads {
//         let ext_data = ext_data.clone();
//         let total_size = total_size.clone();
//         let ext_data_index = ext_data_index.clone();
//         let handle = thread::spawn(move || {
//             // While there is data to check the size, keep checking
//             while ext_data_index.load(Ordering::Relaxed) < ext_data.len() {
//                 let ext_data_index = ext_data_index.fetch_add(1, Ordering::Relaxed);
//                 if ext_data_index >= ext_data.len() {
//                     break;
//                 }
//                 let ext_data_inner = &ext_data[ext_data_index];
//                 total_size.fetch_add(ext_data_inner.len(), Ordering::Relaxed);
//             }
//         });
//         handles.push(handle);
//     }
//     for handle in handles {
//         handle.join().unwrap();
//     }
//     let estimation_time = start.elapsed();
//     println!("Estimation time: {:?}", estimation_time);
//     println!(
//         "Allocating {} gigabytes", 
//         total_size.clone().load(Ordering::Relaxed) as f32 * 8.0 / 1024.0 / 1024.0 / 1024.0
//     );

//     // Allocate the memory
//     let start = Instant::now();
//     let data = Arc::new(HugeAllocation {
//         allocation: UnsafeCell::new(vec![0; total_size.load(Ordering::Relaxed)]),
//     });
//     let allocation_time = start.elapsed();
//     println!("Allocation time: {:?}", allocation_time);
//     // Spawn threads to copy the data
//     let start = Instant::now();
//     let mut handles = vec![];
//     let mut ext_data_index = Arc::new(AtomicUsize::new(0));

//     let expected_num_value = ext_data[0].len();
//     let size_per_value = total_size.load(Ordering::Relaxed) / expected_num_value;
//     let offset_vec = (0..expected_num_value).map(
//         |x| AtomicUsize::new(x * size_per_value)
//     ).collect::<Vec<AtomicUsize>>();
//     let offset_vec = Arc::new(offset_vec);
    
//     for i in 0..num_threads {
//         let data_clone = Arc::clone(&data);
//         let ext_data = ext_data.clone();
//         let ext_data_index = ext_data_index.clone();
//         let offset_vec = offset_vec.clone();
//         let handle = thread::spawn(move || {
//             while ext_data_index.load(Ordering::Relaxed) < ext_data.len() {
//                 let ext_data_index = ext_data_index.fetch_add(1, Ordering::Relaxed);
//                 if ext_data_index >= ext_data.len() {
//                     break;
//                 }
//                 let ext_data_inner = &ext_data[ext_data_index];
//                 // let offset_in_allocation = offset_in_allocation.fetch_add(ext_data_inner.len(), Ordering::Relaxed);
//                 let data = data_clone.allocation.get();
//                 for j in 0..ext_data_inner.len() {
//                     // Get offset from offsets map for value of j
//                     let val = ext_data_inner[j];
//                     let offset_in_allocation = offset_vec[val].fetch_add(1, Ordering::Relaxed);
//                     unsafe {
//                         (*data)[offset_in_allocation] = ext_data_index;
//                     }
//                 }
//             }
//         });
//         handles.push(handle);
//     }
//     for handle in handles {
//         handle.join().unwrap();
//     }
//     let computation_time = start.elapsed();
//     println!("Filling time: {:?}", computation_time);
    
//     for i in 0..10 {
//         println!("{:?}", offset_vec[i].load(Ordering::Relaxed));
//     }
    
//     return data;
// }