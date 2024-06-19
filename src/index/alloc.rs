// File: nested_vec.rs
// Created: 2023-12-28 16:44:31
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

// external crates
use rustc_hash::FxHashMap;
use dashmap::DashMap;


// 
use std::sync::{atomic, Arc};
use std::cell::UnsafeCell;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

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
    pub fn fill(&self, index: usize, data: T) {
        let allocation_ref = unsafe { &mut *self.allocation.get() };
        allocation_ref[index] = data;
    }
}

impl HugeAllocation<usize> {
    pub fn add(&self, index: usize, data: usize) {
        let allocation_ref = unsafe { &mut *self.allocation.get() };
        allocation_ref[index] += data;
    }
}

// Test for HugeAllocation
#[cfg(test)]
mod tests {
    use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

    use crate::measure_time;

    use super::*;
    #[test]
    fn test_huge_allocation() {
        // Making an allocation for 2^30 elements
        let size = 1 << 30;
        let allocation = HugeAllocation::<usize>::new(size);
        let allocation_ref = unsafe { &*allocation.allocation.get() };
        assert_eq!(allocation_ref.len(), size);
    }
    
    #[test]
    fn test_parallel_filling() {
        let size = 1 << 10;
        let allocation = measure_time!(HugeAllocation::<usize>::new(size));
        let allocation_ref = unsafe { &*allocation.allocation.get() };
        // Use rayon
        measure_time!({(0..size).into_par_iter().for_each(|i| {allocation.fill(i, i);});});
        for i in 0..size {
            assert_eq!(allocation_ref[i], i);
        }
    }
    
    #[test]
    fn test_counting_and_filling() {
        // one 1, two 2s, three 3s, four 4s, five 5s, until 32
        let size = 1024;
        let data: Vec<usize> = (1..size+1).flat_map(|i| vec![i; i]).collect();
        // Count each element and fill the vector
        let count_vector = HugeAllocation::<usize>::new(size);
        data.par_iter().for_each(|&i| {
            count_vector.add(i-1, 1);
        });
        // Print 
        
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
        ids: Vec<K>, data: Vec<Vec<V>>,
        num_threads: usize, allocation_size: usize,
        offset_path: String, data_path: String,
    ) -> IndexBuilder<K, V> {

        // If num_threads is 0, set it to default. one liner
        let num_threads = if num_threads == 0 { DEFAULT_NUM_THREADS } else { num_threads };
        // If allocation_size is 0, set it to default
        let allocation_size = if allocation_size == 0 { DEFAULT_ALLOCATION_SIZE } else { allocation_size };
        // Create IndexBuilder
        IndexBuilder {
            offset: DashMap::new(),
            allocation: Arc::new(HugeAllocation::new(allocation_size)),
            ids: Arc::new(ids),
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
        for _i in 0..self.num_threads {
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
                        let mut entry = output_ref_inner.entry(curr_hash.clone()).or_insert(Vec::new());
                        if !entry.contains(&curr_id) {
                            entry.push(curr_id);
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
            let _ids = self.ids.clone();
            let data_index = data_index.clone();
            let _allocation = self.allocation.clone();
            let handle = thread::spawn(move || {
                while data_index.load(Ordering::Relaxed) < data.len() {
                    let data_index = data_index.fetch_add(1, Ordering::Relaxed);
                    if data_index >= data.len() {
                        break;
                    }
                    let data_inner = &data[data_index];
                    for j in 0..data_inner.len() {
                        // Get offset from offsets map for value of j
                        let _val = data_inner[j];
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
        let offset_map = DashMap::new();
        let mut vec: Vec<K> = Vec::new();
        let offset = AtomicUsize::new(0);
        
        orig_map.iter().for_each(|entry| {
            let hash = entry.key().to_owned();
            let ids = entry.value();
            offset_map.insert(
                hash,
                (offset.fetch_add(ids.len(), Ordering::Relaxed), ids.len())
            );
            vec.extend_from_slice(ids);
        });
        (offset_map, vec)
    }

    pub fn convert_sorted_pairs_to_offset_and_values(
        &self,
        sorted_pairs: Vec<(V, K)>
    ) -> (DashMap<V, (usize, usize)>, Vec<K>) {
        // OffsetMap - key: hash, value: (offset, length)
        // Vec - all values concatenated
        let offset_map = DashMap::new();
        let mut vec: Vec<K> = Vec::new();
        let offset = AtomicUsize::new(0);
        
        sorted_pairs.iter().for_each(|pair| {
            // If offset_map does not contain the key, insert it
            if !offset_map.contains_key(&pair.0) {
                offset_map.insert(pair.0, (offset.load(Ordering::Relaxed), 1));
            } else {
                // If offset_map contains the key, increment the offset, size and push the value to vec
                let mut entry = offset_map.get_mut(&pair.0).unwrap();
                entry.1 += 1;
            }
            vec.push(pair.1);
            offset.fetch_add(1, Ordering::Relaxed);
        });
        
        (offset_map, vec)
    }

    // fn fill_inner
    // Write dashmap to file
    // pub fn 

    
    
    
    
}


pub fn convert_sorted_pairs_to_offset_and_values<V: HashableSync, K:HashableSync> (
    sorted_pairs: Vec<(V, K)>
) -> (DashMap<V, (usize, usize)>, Vec<K>) {
    // OffsetMap - key: hash, value: (offset, length)
    // Vec - all values concatenated
    let offset_map = DashMap::new();
    let mut vec: Vec<K> = Vec::new();
    let offset = AtomicUsize::new(0);
    
    sorted_pairs.iter().for_each(|pair| {
        // If offset_map does not contain the key, insert it
        if !offset_map.contains_key(&pair.0) {
            offset_map.insert(pair.0, (offset.load(Ordering::Relaxed), 1));
        } else {
            // If offset_map contains the key, increment the offset, size and push the value to vec
            let mut entry = offset_map.get_mut(&pair.0).unwrap();
            entry.1 += 1;
        }
        vec.push(pair.1);
        offset.fetch_add(1, Ordering::Relaxed);
    });
    
    (offset_map, vec)
}

pub fn estimate_hash_size<V: HashableSync, K:HashableSync> (
    sorted_pairs: &Vec<(V, K)>
) -> (usize, usize) {
    let mut total_hashes = 0usize;
    let mut total_values = 0usize;
    let mut curr = &sorted_pairs[0].0;
    for i in 0..sorted_pairs.len() {
        if curr != &sorted_pairs[i].0 {
            total_hashes += 1;
            curr = &sorted_pairs[i].0;
        }
        total_values += 1;
    }
    (total_hashes, total_values)
}

pub fn convert_sorted_pairs_to_offset_and_values_vec<V: HashableSync, K:HashableSync> (
    sorted_pairs: Vec<(V, K)>
) -> (Vec<(V, usize, usize)>, Vec<K>) {
    let (total_hashes, total_values) = estimate_hash_size(&sorted_pairs);
    let mut offset_list: Vec<(V, usize, usize)> = Vec::with_capacity(total_hashes);
    let mut vec: Vec<K> = Vec::with_capacity(total_values);

    if let Some((first_hash, _)) = sorted_pairs.first() {
        let mut current_hash = first_hash;
        let mut current_offset = 0;
        let mut current_count = 0;

        for (index, pair) in sorted_pairs.iter().enumerate() {
            if pair.0 == *current_hash {
                current_count += 1;
            } else {
                offset_list.push((*current_hash, current_offset, current_count));
                current_hash = &pair.0;
                current_offset = index;
                current_count = 1;
            }
            vec.push(pair.1);
        }

        offset_list.push((*current_hash, current_offset, current_count));
    }
    // Shrink offset_list and vec
    offset_list.shrink_to_fit();
    vec.shrink_to_fit();
    (offset_list, vec)
}
