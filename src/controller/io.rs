// File: io.rs
// Created: 2024-02-20 20:06:01
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use dashmap::DashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write, Error};
use memmap2::Mmap;
use crate::prelude::{GeometricHash, HashType};
use std::mem::size_of;


pub fn save_offset_map(
    path: &str, offset_map: &DashMap<GeometricHash, (usize, usize)>
) -> Result<(), Error> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(path)?;

    // Get hash type
    let hash_type = offset_map.iter().next().unwrap().key().hash_type();
    let total_size = match hash_type.encoding_type() {
        32 => 20 * offset_map.len() as u64,
        64 => 24 * offset_map.len() as u64,
        _ => { panic!("Invalid hash type"); }
    };
    file.set_len(total_size as u64)?;
    // Write as whole
    let mut writer = BufWriter::new(file);
    // Iterate through dashmap
    offset_map.iter().for_each(|entry| {
        let key = entry.key();
        let value = entry.value();
        match key.hash_type().encoding_type() {
            32 => { writer.write_all(&key.as_u32().to_le_bytes()).unwrap(); },
            64 => { writer.write_all(&key.as_u64().to_le_bytes()).unwrap(); },
            _ => { panic!("Invalid hash type"); }
        }
        writer.write_all(&value.0.to_le_bytes()).unwrap();
        writer.write_all(&value.1.to_le_bytes()).unwrap();
    });
    Ok(())
}

pub fn read_offset_map(path: &str, hash_type: HashType) -> Result<DashMap<GeometricHash, (usize, usize)>, Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let mut offset_map = DashMap::default();
    let mut offset = 0;
    while offset < mmap.len() {
        let (key, w): (GeometricHash, usize) = match hash_type.encoding_type() {
            32 => {
                (GeometricHash::from_u32(
                    u32::from_le_bytes(mmap[offset..offset + 4].try_into().unwrap()),
                    hash_type
                ), 4_usize)
            },
            64 => {
                (GeometricHash::from_u64(
                    u64::from_le_bytes(mmap[offset..offset + 8].try_into().unwrap()),
                    hash_type
                ), 8_usize)
            },
            _ => { panic!("Invalid hash type"); }
        };
        let value = (
            usize::from_le_bytes(mmap[offset + w..offset + w + 8].try_into().unwrap()),
            usize::from_le_bytes(mmap[offset + w + 8..offset + w + 16].try_into().unwrap())
        );
        offset_map.insert(key, value);
        offset = offset + w + 16;
    }
    Ok(offset_map)
}


pub fn write_usize_vector(path: &str, vec: &Vec<usize>) -> Result<(), Error> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(path)?;
    let total_size: u64 = 8 * vec.len() as u64;
    file.set_len(total_size as u64)?;
    // Write as whole
    let mut writer = BufWriter::new(file);
    let vec_bytes = unsafe { 
        std::slice::from_raw_parts(vec.as_ptr() as *const u8, 
        total_size as usize) 
    };
    writer.write_all(vec_bytes)?;
    Ok(())
}

pub fn read_usize_vector(path: &str)-> Result<(Mmap, &'static [u64]), Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let vec = unsafe {
        std::slice::from_raw_parts(mmap.as_ptr() as *const u64,
        mmap.len() / size_of::<u64>())
    };
    Ok((mmap, vec))
}

pub fn get_values_with_offset(vec: &[u64], offset: usize, length: usize) -> &[u64] {
    &vec[offset..offset + length]
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_offset_map_io() {
        let mut offset_map = DashMap::new();
        offset_map.insert(GeometricHash::from_u64(0, HashType::PDBMotif), (0, 10));
        offset_map.insert(GeometricHash::from_u64(1, HashType::PDBMotif), (10, 10));
        offset_map.insert(GeometricHash::from_u64(2, HashType::PDBMotif), (20, 10));
        println!("{:?}", offset_map);

        save_offset_map("test_offset_map_io.offset", &offset_map).unwrap();
        let offset_map = read_offset_map("test_offset_map_io.offset", HashType::PDBMotif).unwrap();
        offset_map.iter().for_each(|entry| {
            let key = entry.key();
            let value = entry.value();
            println!("{:?} -> {:?}", key, value);
        });
    }
    #[test]
    fn test_usize_vector_io() {
        let vec = vec![1, 2, 3, 4, 5];
        write_usize_vector("test_usize_vector_io.value", &vec).unwrap();
        let (mmap, vec) = read_usize_vector("test_usize_vector_io.value").unwrap();
        assert_eq!(vec, &[1, 2, 3, 4, 5]);
    }
    
}