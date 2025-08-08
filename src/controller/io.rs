// File: io.rs
// Created: 2024-02-20 20:06:01
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use dashmap::DashMap;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write, Error};
use memmap2::Mmap;
use crate::prelude::{print_log_msg, log_msg, GeometricHash, HashType, FAIL, INFO};
use crate::structure::core::{CompactStructure, Structure};
use crate::{CIFReader, PDBReader};
use std::mem::size_of;

#[cfg(feature="foldcomp")]
use crate::structure::io::fcz::FoldcompDbReader;


pub fn save_offset_map(
    path: &str, offset_map: &DashMap<GeometricHash, (usize, usize)>
) -> Result<(), Error> {
    let file = OpenOptions::new()
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
    let file = OpenOptions::new()
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

pub fn read_offset_map_single(path: &str, hash_type: HashType) -> Result<DashMap<GeometricHash, (usize, usize)>, Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let hash_encoding_type = hash_type.encoding_type();
    let offset_map = DashMap::with_capacity(mmap.len() / (16 + (hash_encoding_type / 8)) as usize);
    let mut offset = 0;
    while offset < mmap.len() {
        let (key, w): (GeometricHash, usize) = match hash_encoding_type {
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

pub fn read_offset_map(path: &str, hash_type: HashType) -> Result<DashMap<GeometricHash, (usize, usize)>, Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let hash_encoding_type = hash_type.encoding_type();
    let chunk_size = 16 + (hash_encoding_type / 8) as usize;
    let chunks: Vec<_> = mmap.chunks(chunk_size).collect();

    // Set number of threads
    // rayon::ThreadPoolBuilder::new().num_threads(8).build_global().unwrap();
    let offset_map: DashMap<_, _> = chunks.into_par_iter().filter_map(|chunk| {
        let (key, w): (GeometricHash, usize) = match hash_encoding_type {
            32 => {
                (GeometricHash::from_u32(
                    u32::from_le_bytes(chunk[0..4].try_into().unwrap()),
                    hash_type
                ), 4_usize)
            },
            64 => {
                (GeometricHash::from_u64(
                    u64::from_le_bytes(chunk[0..8].try_into().unwrap()),
                    hash_type
                ), 8_usize)
            },
            _ => { return None; }
        };
        let value = (
            usize::from_le_bytes(chunk[w..w + 8].try_into().unwrap()),
            usize::from_le_bytes(chunk[w + 8..w + 16].try_into().unwrap())
        );
        Some((key, value))
    }).collect();

    Ok(offset_map)
}



pub fn write_usize_vector(path: &str, vec: &Vec<usize>) -> Result<(), Error> {
    let file = OpenOptions::new()
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
pub fn get_values_with_offset_u24(vec: &[u8], offset: usize, length: usize) -> &[u8] {
    let offset = offset * 3;
    let length = length * 3;
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
        24 => (3 * vec.len()) as u64, // 3 bytes
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
        24 => {
            // Extract 3 bytes from usize as Vec<u8>
            let vec_u8 = vec.iter().flat_map(|&x| {
                // endianess is not considered. first 2bytes: id, last byte: grid index
                let mut bytes = Vec::new();
                bytes.push(((x >> 16) & 0xFF) as u8);
                bytes.push(((x >> 8) & 0xFF) as u8);
                bytes.push((x & 0xFF) as u8);
                bytes
            }).collect::<Vec<u8>>();
            let vec_bytes = unsafe { 
                std::slice::from_raw_parts(vec_u8.as_ptr() as *const u8, 
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

pub fn read_compact_structure(path: &str) -> Result<(CompactStructure, bool), ()> {
    #[cfg(not(feature="foldcomp"))]
    let use_foldcomp = false;
    #[cfg(feature="foldcomp")]
    let use_foldcomp = if path.contains(':') { true } else { false };

    
    #[cfg(not(feature="foldcomp"))]
    let compact_structure = read_structure_from_path(path).expect(
        &log_msg(FAIL, &format!("Failed to read structure from file: {}", &path))
    ).to_compact();
    
    #[cfg(feature="foldcomp")]
    let compact_structure = if !use_foldcomp {
        read_structure_from_path(path).expect(
            &log_msg(FAIL, &format!("Failed to read structure from file: {}", &path))
        ).to_compact()
    } else {
        let mut split = path.split(':');
        let db_path = split.next().unwrap();
        let name = split.next().unwrap();
        let mut foldcomp_db_reader = FoldcompDbReader::new(db_path);
        foldcomp_db_reader.sort_lookup_by_name();
        let structure_io_result = foldcomp_db_reader.read_single_structure(name);
        if let Ok(structure) = structure_io_result {
            structure.to_compact()
        } else {
            return Err(());
        }
    };
    Ok((compact_structure, use_foldcomp))
}


pub fn read_structure_from_path(path: &str) -> Option<Structure> {
    if path.ends_with(".gz") {
        if path.ends_with(".pdb.gz") || path.ends_with(".ent.gz") {
            let reader = PDBReader::from_file(path).expect(
                &log_msg(FAIL, format!("Failed to read PDB file: {}", path).as_str())
            );
            let structure = reader.read_structure_from_gz().expect(
                &log_msg(FAIL, format!("Failed to read structure from PDB file: {}", path).as_str())
            );
            Some(structure)
        } else if path.ends_with(".cif.gz") {
            let reader = CIFReader::from_file(path).expect(
                &log_msg(FAIL, format!("Failed to read CIF file: {}", path).as_str())
            );
            let structure = reader.read_structure_from_gz().expect(
                &log_msg(FAIL, format!("Failed to read structure from CIF file: {}", path).as_str())
            );
            Some(structure)
        } else {
            None
        }
    } else {
        if path.ends_with(".pdb") || path.ends_with(".ent") {
            let reader = PDBReader::from_file(path).expect(
                &log_msg(FAIL, format!("Failed to read PDB file: {}", path).as_str())
            );
            let structure = reader.read_structure().expect(
                &log_msg(FAIL, format!("Failed to read structure from PDB file: {}", path).as_str())
            );
            Some(structure)
        } else if path.ends_with(".cif") {
            let reader = CIFReader::from_file(path).expect(
                &log_msg(FAIL, format!("Failed to read CIF file: {}", path).as_str())
            );
            let structure = reader.read_structure().expect(
                &log_msg(FAIL, format!("Failed to read structure from CIF file: {}", path).as_str())
            );
            Some(structure)
        } else {
            None
        }
    }
}


// Functions to load index files
pub fn get_offset_value_lookup_type(index_path: String) -> (String, String, String, String) {
    let offset_path = format!("{}.offset", index_path.clone());
    let lookup_path = format!("{}.lookup", index_path.clone());
    let hash_type_path = format!("{}.type", index_path.clone());

    // Check for new format (no extension) first, then fall back to .value extension
    let value_path = if std::path::Path::new(&index_path).is_file() {
        index_path.clone()  // New format: no extension
    } else {
        format!("{}.value", index_path.clone())  // Old format: with .value extension
    };

    assert!(std::path::Path::new(&offset_path).is_file());
    assert!(std::path::Path::new(&value_path).is_file());
    assert!(std::path::Path::new(&lookup_path).is_file());
    assert!(std::path::Path::new(&hash_type_path).is_file());
    (offset_path, value_path, lookup_path, hash_type_path)
}

// Functions to load index files
pub fn check_and_get_indices(index_path: Option<String>, verbose: bool) -> Vec<String> {
    // Get path. formatting without quotation marks
    let index_path = index_path.unwrap();
    // Check if index_path_0 is a file.
    let _index_chunk_prefix = format!("{}_0", index_path.clone());
    let index_chunk_path = format!("{}_0.offset", index_path.clone());
    let mut index_paths = Vec::new();
    if std::path::Path::new(&index_chunk_path).is_file() {
        if verbose {
            print_log_msg(INFO, &format!("Index table is chunked"));
        }
        let mut i = 0;
        loop {
            let index_chunk_prefix = format!("{}_{}", index_path.clone(), i);
            let index_chunk_path = format!("{}.offset", index_chunk_prefix);
            if std::path::Path::new(&index_chunk_path).is_file() {
                index_paths.push(index_chunk_prefix);
                i += 1;
            } else {
                break;
            }
        }
    } else {
        index_paths.push(index_path.clone());
    }
    index_paths
}


#[cfg(feature = "foldcomp")]
pub fn get_foldcomp_db_path_with_prefix(prefix: &str) -> Option<String> {
    // If prefix format is like "*_folddisco", use "*" as prefix. 
    // Candidate paths with the prefix: parsed_prefix, parsed_prefix_foldcomp
    // Check if db_path, db_path.index, db_path.lookup exists.
    
    // Parse prefix - remove _folddisco if present
    let parsed_prefix = if prefix.ends_with("_folddisco") {
        prefix.trim_end_matches("_folddisco").to_string()
    } else {
        prefix.to_string()
    };
    
    // Create candidate prefixes
    let candidate_prefixes = vec![
        parsed_prefix.clone(),
        format!("{}_foldcomp", parsed_prefix)
    ];
    // Check each candidate prefix for foldcomp database files
    for candidate in candidate_prefixes {
        if is_valid_foldcomp_db(&candidate) {
            return Some(candidate);
        }
    }
    
    None
}

#[cfg(feature = "foldcomp")]
fn is_valid_foldcomp_db(db_path: &str) -> bool {
    // Check if the required foldcomp database files exist
    let db_file = std::path::Path::new(db_path);
    let index_path = format!("{}.index", db_path);
    let lookup_path = format!("{}.lookup", db_path);
    let index_file = std::path::Path::new(&index_path);
    let lookup_file = std::path::Path::new(&lookup_path);
    db_file.is_file() && index_file.is_file() && lookup_file.is_file()
}

pub fn default_index_path(input_path: &str) -> String {
    // If input_path ends with "_foldcomp", remove it and add "_folddisco"
    // Else, return input_path + "_folddisco"
    if input_path.ends_with("_foldcomp") {
        input_path.replace("_foldcomp", "_folddisco")
    } else {
        format!("{}_folddisco", input_path)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_offset_map_io() {
        let offset_map = DashMap::new();
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
        let (_mmap, vec) = read_usize_vector("test_usize_vector_io.value").unwrap();
        assert_eq!(vec, &[1, 2, 3, 4, 5]);
    }
}