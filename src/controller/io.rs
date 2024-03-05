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

pub fn save_offset_vec(
    path: &str, offset_map: &Vec<(GeometricHash, usize, usize)>
) -> Result<(), Error> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(path)?;

    // Get hash type
    let hash_type = offset_map[0].0.hash_type();
    let total_size = match hash_type.encoding_type() {
        32 => 20 * offset_map.len() as u64,
        64 => 24 * offset_map.len() as u64,
        _ => { panic!("Invalid hash type"); }
    };
    file.set_len(total_size as u64)?;
    // Write as whole
    let mut writer = BufWriter::new(file);
    // Iterate through vector
    offset_map.iter().for_each(|(key, offset, length)| {
        match key.hash_type().encoding_type() {
            32 => { writer.write_all(&key.as_u32().to_le_bytes()).unwrap(); },
            64 => { writer.write_all(&key.as_u64().to_le_bytes()).unwrap(); },
            _ => { panic!("Invalid hash type"); }
        }
        writer.write_all(&offset.to_le_bytes()).unwrap();
        writer.write_all(&length.to_le_bytes()).unwrap();
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
    // let total_size: u64 = 8 * vec.len() as u64;
    let total_size: u64 = (size_of::<usize>() * vec.len()) as u64;
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
pub fn get_values_with_offset_u16(vec: &[u16], offset: usize, length: usize) -> &[u16] {
    &vec[offset..offset + length]
}
pub fn get_values_with_offset_u32(vec: &[u32], offset: usize, length: usize) -> &[u32] {
    &vec[offset..offset + length]
}
pub fn get_values_with_offset_u8(vec: &[u8], offset: usize, length: usize) -> &[u8] {
    &vec[offset..offset + length]
}

pub fn write_usize_vector_in_bits(path: &str, vec: &Vec<usize>, num_bits: usize) -> Result<(), Error> {
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(path)?;
    // let total_size: u64 = 8 * vec.len() as u64;
    let total_size: u64 = match num_bits {
        8 => (size_of::<u8>() * vec.len()) as u64,
        16 => (size_of::<u16>() * vec.len()) as u64,
        32 => (size_of::<u32>() * vec.len()) as u64,
        _ => { panic!("Invalid number of bits"); }
    };
    file.set_len(total_size as u64)?;
    // Write after converting usize to matched integer type
    let mut writer = BufWriter::new(file);
    match num_bits {
        8 => {
            let vec_u8 = vec.iter().map(|&x| x as u8).collect::<Vec<u8>>();
            let vec_bytes = unsafe { 
                std::slice::from_raw_parts(vec_u8.as_ptr() as *const u8, 
                total_size as usize) 
            };
            writer.write_all(vec_bytes)?;
        },
        16 => {
            let vec_u16 = vec.iter().map(|&x| x as u16).collect::<Vec<u16>>();
            let vec_bytes = unsafe { 
                std::slice::from_raw_parts(vec_u16.as_ptr() as *const u8, 
                total_size as usize) 
            };
            writer.write_all(vec_bytes)?;
        },
        32 => {
            let vec_u32 = vec.iter().map(|&x| x as u32).collect::<Vec<u32>>();
            let vec_bytes = unsafe { 
                std::slice::from_raw_parts(vec_u32.as_ptr() as *const u8, 
                total_size as usize) 
            };
            writer.write_all(vec_bytes)?;
        },
        64 => {
            let vec_bytes = unsafe { 
                std::slice::from_raw_parts(vec.as_ptr() as *const u8,
                total_size as usize) 
            };
            writer.write_all(vec_bytes)?;
        },
        _ => { panic!("Invalid number of bits"); }
    }    
    Ok(())
}

pub fn read_u8_vector(path: &str)-> Result<(Mmap, &'static [u8]), Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let vec = unsafe {
        std::slice::from_raw_parts(mmap.as_ptr() as *const u8,
        mmap.len() / size_of::<u8>())
    };
    Ok((mmap, vec))
}
pub fn read_u16_vector(path: &str)-> Result<(Mmap, &'static [u16]), Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let vec = unsafe {
        std::slice::from_raw_parts(mmap.as_ptr() as *const u16,
        mmap.len() / size_of::<u16>())
    };
    Ok((mmap, vec))
}
pub fn read_u32_vector(path: &str)-> Result<(Mmap, &'static [u32]), Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let vec = unsafe {
        std::slice::from_raw_parts(mmap.as_ptr() as *const u32,
        mmap.len() / size_of::<u32>())
    };
    Ok((mmap, vec))
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