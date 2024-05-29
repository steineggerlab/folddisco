use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Write};
use std::mem;
use std::path::Path;
use std::slice;
use std::sync::atomic::{AtomicUsize, Ordering};
use memmap2::Mmap;

use crate::index::alloc::estimate_hash_size;
use crate::measure_time;
use crate::prelude::GeometricHash;

// const INITIAL_CAPACITY: usize = 16;

#[derive(Debug)]
struct BitVec {
    bits: Vec<u8>,
    len: usize,
}

impl BitVec {
    fn new(size: usize) -> Self {
        let byte_size = (size + 7) / 8; // Round up to the nearest byte
        BitVec {
            bits: vec![0; byte_size],
            len: size,
        }
    }

    fn set(&mut self, index: usize, value: bool) {
        let byte_index = index / 8;
        let bit_index = index % 8;
        if value {
            self.bits[byte_index] |= 1 << bit_index;
        } else {
            self.bits[byte_index] &= !(1 << bit_index);
        }
    }

    fn get(&self, index: usize) -> bool {
        let byte_index = index / 8;
        let bit_index = index % 8;
        self.bits[byte_index] & (1 << bit_index) != 0
    }
}

#[derive(Debug)]
pub struct SimpleHashMap {
    buckets: Vec<u32>,             // Stores the hash indices with occupancy information
    occupancy: BitVec,               // Stores the occupancy information
    keys: Vec<u32>,        // Stores keys separately
    values: Vec<(usize, usize)>,     // Stores values separately
    size: usize,
    capacity: usize,
}

impl SimpleHashMap {
    fn new(capacity: usize) -> Self {
        SimpleHashMap {
            buckets: vec![0; capacity],
            occupancy: BitVec::new(capacity),
            keys: Vec::with_capacity(capacity),
            values: Vec::with_capacity(capacity),
            size: 0,
            capacity,
        }
    }

    fn new_from_std_hashmap(mut map: std::collections::HashMap<GeometricHash, (usize, usize)>, capacity: usize) -> Self {
        let mut simple_hash_map = SimpleHashMap::new(capacity);

        for (key, value) in map.drain() {
            simple_hash_map.insert(key, value);
        }

        simple_hash_map
    }

    pub fn new_from_vec(vec: Vec<(GeometricHash, usize, usize)>) -> Self {
        let capacity = vec.len();
        let mut simple_hash_map = SimpleHashMap::new(capacity);
        for (key, value1, value2) in vec {
            simple_hash_map.insert(key, (value1, value2));
        }
        simple_hash_map
    }
    
    fn new_from_dashmap(map: dashmap::DashMap<GeometricHash, (usize, usize)>, capacity: usize) -> Self {
        let mut simple_hash_map = SimpleHashMap::new(capacity);

        map.iter().for_each(|entry| {
            simple_hash_map.insert(entry.key().clone(), entry.value().clone());
        });
        
        simple_hash_map
    }

    fn hash(&self, hash: u32) -> usize {
        (hash % self.capacity as u32) as usize
    }

    fn insert(&mut self, key: GeometricHash, value: (usize, usize)) {
        let key_u32 = key.as_u32();
        let hash = self.hash(key_u32);
        let mut index = hash;

        loop {
            if !self.occupancy.get(index) {
                let key_index = self.keys.len();
                self.keys.push(key_u32);
                self.values.push(value);
                self.buckets[index] = key_index as u32;
                self.occupancy.set(index, true);
                self.size += 1;
                return;
            } else if self.keys[self.buckets[index] as usize] == key_u32 {
                let key_index = self.buckets[index] as usize;
                self.values[key_index] = value;
                return;
            } else {
                index = (index + 1) % self.capacity;
            }
        }
    }

    pub fn get(&self, key: &GeometricHash) -> Option<&(usize, usize)> {
        let key_u32 = key.as_u32();
        let hash = self.hash(key_u32); // Assuming `hash_u32` is your perfect hash function
        let mut index = hash;
        let mut count = 0;
        loop {
            if !self.occupancy.get(index) {
                return None;
            } else if self.keys[self.buckets[index] as usize] == key_u32 {
                let key_index = self.buckets[index] as usize;
                return Some(&self.values[key_index]);
            } else {
                index = (index + 1) % self.capacity;
            }
            count += 1;
            if count == self.capacity {
                return None;
            }
        }
    }
    
    fn estimate_file_size(&self) -> usize {
        let mut size = 0;
        size += mem::size_of::<usize>() * 2; // size and capacity
        size += self.buckets.len() * mem::size_of::<u32>();
        size += (self.capacity + 7) / 8; // occupancy
        size += self.keys.len() * mem::size_of::<u32>();
        size += self.values.len() * mem::size_of::<(usize, usize)>();
        size
    }
    
    
    pub fn dump_to_disk(&self, path: &Path) -> std::io::Result<()> {
        // If file exists, make a new file
        let file = OpenOptions::new().read(true).write(true).create(true).open(path)?;
        file.set_len(self.estimate_file_size() as u64)?;
        let mut writer = BufWriter::new(file);
        // Serialize and write metadata
        writer.write_all(&self.size.to_le_bytes())?;
        writer.write_all(&self.capacity.to_le_bytes())?;

        println!("[DEBUG] Size and capacity written");
        // Serialize and write buckets. 
        let buckets_size = self.buckets.len() * mem::size_of::<u32>();
        // Byte ordering should be preserved. little endian
        let buckets_bytes = unsafe {
            slice::from_raw_parts(self.buckets.as_ptr() as *const u8, buckets_size)
        };
        let buckets_bytes = unsafe {
            slice::from_raw_parts(self.buckets.as_ptr() as *const u8, buckets_size)
        };
        writer.write_all(buckets_bytes)?;
        println!("[DEBUG] Buckets written");
        // Serialize and write occupancy
        let occupancy_size = self.occupancy.bits.len();
        let occupancy_bytes = unsafe {
            slice::from_raw_parts(self.occupancy.bits.as_ptr() as *const u8, occupancy_size)
        };
        writer.write_all(occupancy_bytes)?;
        println!("[DEBUG] Occupancy written");
        // Serialize and write keys
        let keys_size = self.keys.len() * mem::size_of::<u32>();
        let keys_bytes = unsafe {
            slice::from_raw_parts(self.keys.as_ptr() as *const u8, keys_size)
        };
        // Append keys to the file
        writer.write_all(keys_bytes)?;
        println!("[DEBUG] Keys written");
        // Serialize and write values
        let values_size = self.values.len() * mem::size_of::<(usize, usize)>();
        let values_bytes = unsafe {
            slice::from_raw_parts(self.values.as_ptr() as *const u8, values_size)
        };
        writer.write_all(values_bytes)?;
        println!("[DEBUG] Values written");
        Ok(())
    }

    pub fn load_from_disk(path: &Path) -> std::io::Result<Self> {
        // Open as read-only
        let file = OpenOptions::new().read(true).open(path)?;
        println!("[DEBUG] File opened");
        let mmap = unsafe { Mmap::map(&file)? };
        println!("[DEBUG] Mmap length: {}", mmap.len());
        println!("[DEBUG] Mapped file");
        let mut offset = 0usize;
        // Deserialize and read metadata
        let mut size_bytes = [0u8; mem::size_of::<usize>()];
        size_bytes.copy_from_slice(&mmap[offset..offset + mem::size_of::<usize>()]);
        let size = usize::from_le_bytes(size_bytes);
        offset += mem::size_of::<usize>();
        println!("[DEBUG] Size: {}", size);
        let mut capacity_bytes = [0u8; mem::size_of::<usize>()];
        capacity_bytes.copy_from_slice(&mmap[offset..offset + mem::size_of::<usize>()]);
        let capacity = usize::from_le_bytes(capacity_bytes);
        offset += mem::size_of::<usize>();
        println!("[DEBUG] Capacity: {}", capacity);
        let buckets_size = capacity * mem::size_of::<u32>();
        let buckets = unsafe {
            slice::from_raw_parts(mmap.as_ptr().add(offset) as *const u32, capacity).to_vec()
        };
        offset += buckets_size;
        println!("[DEBUG] Current offset: {}", offset);
        println!("[DEBUG] Buckets read");
        // Deserialize and read occupancy
        let occupancy_size = (capacity + 7) / 8;
        let bits = unsafe {
            slice::from_raw_parts(mmap.as_ptr().add(offset) as *const u8, occupancy_size).to_vec()
        };
        let occupancy = BitVec { bits, len: capacity };
        offset += occupancy_size;
        println!("[DEBUG] Current offset: {}", offset);
        println!("[DEBUG] Occupancy read");
        let keys_count = size;
        let keys_size = keys_count * mem::size_of::<u32>();
        let keys = unsafe {
            slice::from_raw_parts(mmap.as_ptr().add(offset) as *const u32, keys_count).to_vec()
        };
        offset += keys_size;
        println!("[DEBUG] Current offset: {}", offset);
        println!("[DEBUG] Keys read");
        let values_size = keys_count * mem::size_of::<(usize, usize)>();
        let values = unsafe {
            slice::from_raw_parts(mmap.as_ptr().add(offset) as *const (usize, usize), keys_count).to_vec()
        };
        println!("[DEBUG] Current offset: {}", offset);
        println!("[DEBUG] Values read");
        Ok(SimpleHashMap {
            buckets,
            occupancy,
            keys,
            values,
            size,
            capacity,
        })
    }
}


pub fn convert_sorted_hash_pairs_to_simplemap(
    sorted_pairs: Vec<(GeometricHash, usize)>
) -> (SimpleHashMap, Vec<usize>) {
    // OffsetMap - key: hash, value: (offset, length)
    // Vec - all values concatenated
    let (total_hashes, total_values) = estimate_hash_size(&sorted_pairs);
    println!("Total hashes: {}, Total values: {}", total_hashes, total_values);
    let mut offset_map = SimpleHashMap::new(total_hashes * 3);
    let mut vec: Vec<usize> = Vec::with_capacity(total_values);
    let offset = AtomicUsize::new(0);
    measure_time!(
    // sorted_pairs.iter().for_each(|pair| {
    //     // If offset_map does not contain the key, insert it
    //     if !offset_map.contains_key(&pair.0) {
    //         offset_map.insert(pair.0, (offset.load(Ordering::Relaxed), 1));
    //     } else {
    //         // If offset_map contains the key, increment the offset, size and push the value to vec
    //         let mut entry = offset_map.get_mut(&pair.0).unwrap();
    //         entry.1 += 1;
    //     }
    //     vec.push(pair.1);
    //     offset.fetch_add(1, Ordering::Relaxed);
    // })
    
    if let Some((first_hash, _)) = sorted_pairs.first() {
        let mut current_hash = first_hash;
        let mut current_offset = 0;
        let mut current_count = 0;

        for (index, pair) in sorted_pairs.iter().enumerate() {
            if pair.0 == *current_hash {
                current_count += 1;
            } else {
                // offset_list.push((*current_hash, current_offset, current_count));
                offset_map.insert(*current_hash, (current_offset, current_count));
                current_hash = &pair.0;
                current_offset = index;
                current_count = 1;
            }
            vec.push(pair.1);
        }
        // offset_list.push((*current_hash, current_offset, current_count));
        offset_map.insert(*current_hash, (current_offset, current_count));
    }
    );
    // Convert the offset_map to SimpleHashMap
    // let offset_map = measure_time!(SimpleHashMap::new_from_std_hashmap(offset_map, total_values));
    (offset_map, vec)
}


#[cfg(test)]
mod tests {
    use dashmap::DashMap;

    use crate::measure_time;
    use crate::prelude::{read_offset_map, save_offset_map, GeometricHash};

    use super::*;
    use std::collections::HashMap as StdHashMap;
    use std::path::PathBuf;

    #[test]
    fn test_dump_and_load() {
        let mut std_map = StdHashMap::new();
        std_map.insert(GeometricHash::from_u32(2u32, crate::prelude::HashType::PDBTrRosetta), (200usize, 200usize));
        std_map.insert(GeometricHash::from_u32(1u32, crate::prelude::HashType::PDBTrRosetta), (100usize, 100usize));


        let map = SimpleHashMap::new_from_std_hashmap(std_map, 16usize);
        println!("MAP: {:?}", map);
        let path = PathBuf::from("hashmap.dat");

        map.dump_to_disk(&path).expect("Failed to dump to disk");
        // Change the permissions of the file to allow read and write access

        let loaded_map = SimpleHashMap::load_from_disk(&path).expect("Failed to load from disk");
        println!("LOADED: {:?}", loaded_map);
        assert_eq!(loaded_map.get(&GeometricHash::from_u32(1u32, crate::prelude::HashType::PDBTrRosetta)), Some(&(100usize, 100usize)));
        assert_eq!(loaded_map.get(&GeometricHash::from_u32(2u32, crate::prelude::HashType::PDBTrRosetta)), Some(&(200usize, 200usize)));
        assert_eq!(loaded_map.get(&GeometricHash::from_u32(3u32, crate::prelude::HashType::PDBTrRosetta)), None);
        std::fs::remove_file(path).expect("Failed to remove test file");
    }
    
    #[test]
    fn test_bigger() {
        let test_size = 200000usize;
        let mut std_map = DashMap::new();
        for i in 0..test_size {
            std_map.insert(GeometricHash::from_u32(i as u32, crate::prelude::HashType::PDBTrRosetta), (i, i));
        }
        let dash_map = std_map.clone();
        let map = SimpleHashMap::new_from_dashmap(std_map, test_size);
        let path = PathBuf::from("hashmap.dat");

        map.dump_to_disk(&path).expect("Failed to dump to disk");
        // Change the permissions of the file to allow read and write access
        measure_time!({
            let loaded_map = SimpleHashMap::load_from_disk(&path).expect("Failed to load from disk");
        });
        let loaded_map = SimpleHashMap::load_from_disk(&path).expect("Failed to load from disk");
        for i in 0..20 {
            assert_eq!(loaded_map.get(&GeometricHash::from_u32(i as u32, crate::prelude::HashType::PDBTrRosetta)), Some(&(i, i)));
        }
        save_offset_map("hashmap.offset", &dash_map).unwrap();
        measure_time!({
            let offset_map = read_offset_map("hashmap.offset", crate::prelude::HashType::PDBTrRosetta).unwrap();
        });
        // Delete the file
        std::fs::remove_file(path).expect("Failed to remove test file");
        std::fs::remove_file("hashmap.offset").expect("Failed to remove test file");
    }

    #[test]
    fn test_conversion_from_vector() {
        println!("Test conversion from vector");
        let vec = vec![
            (GeometricHash::from_u32(1999u32, crate::prelude::HashType::PDBTrRosetta), 100usize, 100usize),
            (GeometricHash::from_u32(22345234u32, crate::prelude::HashType::PDBTrRosetta), 200usize, 200usize),
        ];
        let map = SimpleHashMap::new_from_vec(vec);
        println!("{:?}", map);
        assert_eq!(map.get(&GeometricHash::from_u32(1999u32, crate::prelude::HashType::PDBTrRosetta)), Some(&(100usize, 100usize)));
        assert_eq!(map.get(&GeometricHash::from_u32(22345234u32, crate::prelude::HashType::PDBTrRosetta)), Some(&(200usize, 200usize)));
        assert_eq!(map.get(&GeometricHash::from_u32(3u32, crate::prelude::HashType::PDBTrRosetta)), None);
    }
}
