use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use std::mem;
use std::path::Path;
use std::slice;
use memmap2::{Mmap, MmapMut};

use crate::HashableSync;

const INITIAL_CAPACITY: usize = 100;
const OCCUPIED_MASK: usize = 1 << (mem::size_of::<usize>() * 8 - 1); // Highest bit mask

#[derive(Debug)]
struct SimpleHashMap<K, V> {
    buckets: Vec<usize>, // Stores the hash indices with occupancy information
    keys: Vec<K>,        // Stores keys separately
    values: Vec<V>,      // Stores values separately
    size: usize,
    capacity: usize,
}

impl<K, V> SimpleHashMap<K, V>
where
    K: HashableSync,
    V: Clone,
{
    fn new(capacity: usize) -> Self {
        SimpleHashMap {
            buckets: vec![0; capacity],
            keys: Vec::with_capacity(capacity),
            values: Vec::with_capacity(capacity),
            size: 0,
            capacity,
        }
    }

    fn new_from_std_hashmap(map: std::collections::HashMap<K, V>, capacity: usize) -> Self {
        let mut simple_hash_map = SimpleHashMap::new(capacity);

        for (key, value) in map {
            simple_hash_map.insert(key, value);
        }

        simple_hash_map
    }

    fn hash(&self, hash: u32) -> usize {
        (hash % self.capacity as u32) as usize
    }

    fn is_occupied(entry: usize) -> bool {
        entry & OCCUPIED_MASK != 0
    }

    fn mark_occupied(index: usize) -> usize {
        index | OCCUPIED_MASK
    }

    fn get_index(entry: usize) -> usize {
        entry & !OCCUPIED_MASK
    }

    fn insert(&mut self, key: K, value: V) {
        let hash = self.hash(key.hash_u32()); // Assuming `hash_u32` is your perfect hash function
        let mut index = hash;

        loop {
            if !Self::is_occupied(self.buckets[index]) {
                let key_index = self.keys.len();
                self.keys.push(key);
                self.values.push(value);
                self.buckets[index] = Self::mark_occupied(key_index);
                self.size += 1;
                return;
            } else if self.keys[Self::get_index(self.buckets[index])] == key {
                let key_index = Self::get_index(self.buckets[index]);
                self.values[key_index] = value;
                return;
            } else {
                index = (index + 1) % self.capacity;
            }
        }
    }

    pub fn get(&self, key: &K) -> Option<&V> {
        let hash = self.hash(key.hash_u32()); // Assuming `hash_u32` is your perfect hash function
        let mut index = hash;

        loop {
            if !Self::is_occupied(self.buckets[index]) {
                return None;
            } else if self.keys[Self::get_index(self.buckets[index])] == *key {
                let key_index = Self::get_index(self.buckets[index]);
                return Some(&self.values[key_index]);
            } else {
                index = (index + 1) % self.capacity;
            }
        }
    }

    pub fn dump_to_disk(&self, path: &Path) -> std::io::Result<()> {
        let mut file = OpenOptions::new().write(true).create(true).open(path)?;

        // Serialize and write buckets
        let buckets_size = self.buckets.len() * mem::size_of::<usize>();
        let buckets_bytes = unsafe {
            slice::from_raw_parts(self.buckets.as_ptr() as *const u8, buckets_size)
        };
        file.write_all(buckets_bytes)?;

        // Serialize and write keys
        let keys_size = self.keys.len() * mem::size_of::<K>();
        let keys_bytes = unsafe {
            slice::from_raw_parts(self.keys.as_ptr() as *const u8, keys_size)
        };
        println!("Keys size: {:?}", keys_size);
        println!("Keys: {:?}", self.keys);
        // Append keys to the file
        file.write_all(keys_bytes)?;
        
        // Serialize and write values
        let values_size = self.values.len() * mem::size_of::<V>();
        let values_bytes = unsafe {
            slice::from_raw_parts(self.values.as_ptr() as *const u8, values_size)
        };
        println!("Values size: {:?}", values_size);
        file.write_all(values_bytes)?;

        Ok(())
    }

    pub fn load_from_disk(path: &Path, keys_size: usize) -> std::io::Result<Self> {
        // Open as read-only
        // Use open options.
        println!("Loading from disk");
        let mut file = OpenOptions::new().read(true).open(path)?;
        println!("File opened");
        let mmap = unsafe { Mmap::map(&file)? };
        println!("Mmap created");
        let buckets_size = INITIAL_CAPACITY * mem::size_of::<usize>();
        let buckets = unsafe {
            slice::from_raw_parts(mmap.as_ptr() as *const usize, INITIAL_CAPACITY).to_vec()
        };
        println!("Buckets loaded");
        let mut offset = buckets_size;
        // let keys_size = mmap.len() - buckets_size - mem::size_of::<V>() * INITIAL_CAPACITY;
        let keys_count = keys_size / mem::size_of::<K>();
        println!("Keys size: {:?}", keys_size);
        let keys = unsafe {
            slice::from_raw_parts(mmap.as_ptr().add(offset) as *const K, keys_count).to_vec()
        };
        offset += keys_size * mem::size_of::<K>();
        println!("Keys loaded");
        let values = unsafe {
            slice::from_raw_parts(mmap.as_ptr().add(offset) as *const V, keys_count).to_vec()
        };
        println!("Values loaded");
        Ok(SimpleHashMap {
            buckets,
            keys,
            values,
            size: keys_count,
            capacity: INITIAL_CAPACITY,
        })
    }
}

// trait PerfectHash {
//     fn hash_u32(&self) -> u32;
// }

// // Implement this trait for your key type
// impl PerfectHash for YourKeyType {
//     fn hash_u32(&self) -> u32 {
//         // Your perfect hash function logic here
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap as StdHashMap;
    use std::os::unix::fs::PermissionsExt;
    use std::path::PathBuf;

    #[test]
    fn test_get() {
        let mut std_map = StdHashMap::new();
        std_map.insert(1u32, 100usize);
        std_map.insert(2u32, 200usize);
        let map = SimpleHashMap::new_from_std_hashmap(std_map, INITIAL_CAPACITY);
        println!("{:?}", map);

        assert_eq!(map.get(&1u32), Some(&100usize));
        assert_eq!(map.get(&2u32), Some(&200usize));
    }

    #[test]
    fn test_dump_and_load() {
        let mut std_map = StdHashMap::new();
        std_map.insert(1u32, 100usize);
        std_map.insert(2u32, 200usize);

        let map = SimpleHashMap::new_from_std_hashmap(std_map, INITIAL_CAPACITY);
        let path = PathBuf::from("hashmap.dat");

        map.dump_to_disk(&path).expect("Failed to dump to disk");
        // Change the permissions of the file to allow read and write access

        let loaded_map = SimpleHashMap::<u32, usize>::load_from_disk(&path, 2).expect("Failed to load from disk");
        println!("{:?}", loaded_map);
        assert_eq!(loaded_map.get(&1u32), Some(&100usize));
        assert_eq!(loaded_map.get(&2u32), Some(&200usize));
        assert_eq!(loaded_map.get(&3u32), None);
        std::fs::remove_file(path).expect("Failed to remove test file");
    }
}
