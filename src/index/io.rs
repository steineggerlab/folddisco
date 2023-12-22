use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write, Error};
use std::time::Instant;
use std::fs::OpenOptions;
use memmap2::{Mmap, MmapMut};
use std::collections::HashMap;
use dashmap::DashMap;
use std::mem::{self, size_of};
// Using FxHashMap
use rustc_hash::FxHashMap;

pub fn write_u64_vector(path: &str, vec: &Vec<u64>) -> Result<(), Error> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(path)?;
    let total_size: u64 = 8 * vec.len() as u64;
    file.set_len(total_size as u64)?;
    // Write as whole
    let mut writer = BufWriter::new(file);
    let vec_bytes = unsafe { std::slice::from_raw_parts(vec.as_ptr() as *const u8, total_size as usize) };
    writer.write_all(vec_bytes)?;
    Ok(())
}

pub fn read_u64_vector(path: &str) -> Result<&'static [u64], Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let vec = unsafe { std::slice::from_raw_parts(mmap.as_ptr() as *const u64, mmap.len() / size_of::<u64>()) };
    let boxed = Box::new(mmap);
    let _ = Box::leak(boxed); // Leak the memory to prevent it from being dropped
    Ok(vec)
}

pub fn read_u64_vector_with_mmap(path: &str)-> Result<(Mmap, &'static [u64]), Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let vec = unsafe { std::slice::from_raw_parts(mmap.as_ptr() as *const u64, mmap.len() / size_of::<u64>()) };
    Ok((mmap, vec))
}

// pub fn get_hashmap_size(map: &DashMap<u64, Vec<u64>>) -> usize {
//     let mut size = 0;
//     for (key, value) in map {
//         size += mem::size_of::<u64>() * value.len();
//     }
//     size
// }

// This function consumes the original hashmap
pub fn convert_hashmap_to_offset_and_values(orig_map: DashMap<u64, Vec<u64>>) -> (FxHashMap<u64, (usize, usize)>, Vec<u64>) {
    // OffsetMap - key: hash, value: (offset, length)
    // Vec - all values concatenated
    let mut offset_map = FxHashMap::default();
    let mut vec: Vec<u64> = Vec::new();
    let mut offset: usize = 0;
    for (key, value) in orig_map {
        offset_map.insert(key, (offset, value.len()));
        offset += value.len();
        vec.extend(value);
    }
    (offset_map, vec)
}

pub fn get_values_with_offset(vec: &[u64], offset: usize, length: usize) -> &[u64] {
    &vec[offset..offset + length]
}

pub fn save_offset_map(path: &str, offset_map: &FxHashMap<u64, (usize, usize)>) -> Result<(), Error> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(path)?;
    let total_size: u64 = 24 * offset_map.len() as u64;
    file.set_len(total_size as u64)?;
    // Write as whole
    let mut writer = BufWriter::new(file);
    for (key, value) in offset_map {
        writer.write_all(&key.to_le_bytes())?;
        writer.write_all(&value.0.to_le_bytes())?;
        writer.write_all(&value.1.to_le_bytes())?;
    }
    Ok(())
}

pub fn read_offset_map(path: &str) -> Result<FxHashMap<u64, (usize, usize)>, Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let mut offset_map = FxHashMap::default();
    let mut offset = 0;
    while offset < mmap.len() {
        let key = u64::from_le_bytes(mmap[offset..offset + 8].try_into().unwrap());
        let value = (usize::from_le_bytes(mmap[offset + 8..offset + 16].try_into().unwrap()), usize::from_le_bytes(mmap[offset + 16..offset + 24].try_into().unwrap()));
        offset_map.insert(key, value);
        offset += 24;
    }
    Ok(offset_map)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_write_u64_vector() {
        let mut vec = Vec::new();
        for i in 0..1000000 {
            vec.push(i);
        }
        let start = Instant::now();
        write_u64_vector("test_write_u64_vector", &vec).unwrap();
        println!("write_u64_vector: {:?}", start.elapsed());
    }
    #[test]
    fn test_read_u64_vector() {
        let start = Instant::now();
        let vec = read_u64_vector("test_write_u64_vector").unwrap();
        println!("read_u64_vector: {:?}", start.elapsed());
        assert_eq!(vec.len(), 1000000);
    }
    #[test]
    fn test_read_u64_vector_with_mmap() {
        let start = Instant::now();
        let (mmap, vec) = read_u64_vector_with_mmap("test_write_u64_vector").unwrap();
        println!("read_u64_vector_with_mmap: {:?}", start.elapsed());
        assert_eq!(vec.len(), 1000000);
        assert_eq!(mmap.len(), 1000000 * size_of::<u64>());
    }

    #[test]
    fn test_offset_is_working() {
        let mut orig_map = HashMap::new();
        orig_map.insert(1, vec![1, 2, 3, 4, 5]);
        orig_map.insert(2, vec![6, 7, 8]);
        orig_map.insert(3, vec![9, 10, 11, 12]);
        let start = Instant::now();
        let (offset_map, vec) = convert_hashmap_to_offset_and_values(orig_map);
        println!("convert_hashmap_to_offset_and_values: {:?}", start.elapsed());
        let start = Instant::now();
        assert_eq!(get_values_with_offset(&vec, 0, 5), &[1, 2, 3, 4, 5]);
        println!("get_values_with_offset: {:?}", start.elapsed());
        assert_eq!(get_values_with_offset(&vec, 5, 3), &[6, 7, 8]);
        assert_eq!(get_values_with_offset(&vec, 8, 4), &[9, 10, 11, 12]);
    }

    fn huge_hashmap(size: usize) -> HashMap<u64, Vec<u64>> {
        let mut orig_map = HashMap::new();
        for i in 1..size {
            orig_map.insert(i as u64, vec![i as u64; i]);
        }
        orig_map
    }
    
    #[test]
    fn test_offset_io_with_time() {
        // Make a huge hashmap and measure the running time
        let orig_map = huge_hashmap(10000 * 3 + 1);
        let size = get_hashmap_size(&orig_map);
        println!("Original map size - {}, total {} keys", get_hashmap_size(&orig_map), orig_map.len());
        // Get how much memory it takes
        println!("Total {:.3} GB", size as f64 / 1024.0 / 1024.0 / 1024.0);
        let start = Instant::now();
        let (offset_map, vec) = convert_hashmap_to_offset_and_values(orig_map);
        println!("convert_hashmap_to_offset_and_values: {:?}", start.elapsed());
        // Save offset map
        let start = Instant::now();
        save_offset_map("test_map.offset", &offset_map).unwrap();
        println!("save_offset_map: {:?}", start.elapsed());
        // Load offset map
        let start = Instant::now();
        let offset_map = read_offset_map("test_map.offset").unwrap();
        println!("load_offset_map: {:?}", start.elapsed());
        // Save vector
        let start = Instant::now();
        write_u64_vector("test_map.bin", &vec).unwrap();
        println!("save_u64_vector: {:?}", start.elapsed());

        // Load vector
        let start = Instant::now();
        let vec = read_u64_vector("test_map.bin").unwrap();
        println!("load_u64_vector: {:?}", start.elapsed());

        // Test if offset is working
        let start = Instant::now();
        for i in 1..5 {
            let offset = offset_map.get(&(i as u64)).unwrap();
            let values = get_values_with_offset(&vec, offset.0, offset.1);
            match i {
                1 => assert_eq!(values, &[1]),
                2 => assert_eq!(values, &[2, 2]),
                3 => assert_eq!(values, &[3, 3, 3]),
                4 => assert_eq!(values, &[4, 4, 4, 4]),
                _ => panic!("Should not happen"),
            }
        }
        println!("get_values_with_offset (4 times): {:?}", start.elapsed());

    }
    
    
}