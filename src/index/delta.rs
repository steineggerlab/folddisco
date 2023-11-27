use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write, Error};
use std::time::Instant;
use std::fs::OpenOptions;
use memmap2::{Mmap, MmapMut};
use std::collections::HashMap;
use std::mem::{self, size_of};

// Write/Read hashmap using memmap2 and delta encoding

pub fn delta_encoding_bytes_from_u64(num: u64) -> Vec<u8> {
    // If num is smaller than 2^15, return two bytes
    if num < 0x8000 {
        // Convert to big endian
        let mut bytes = num.to_be_bytes();
        // Set first bit of first byte to 1
        bytes[6] |= 0x80;
        return vec![bytes[6], bytes[7]];
    } else if num < 0x40000000 {
        // If num is smaller than 2^30, return four bytes
        // Convert to big endian
        // Get bytes by 15 bits
        let mut first_bytes = (num >> 15).to_be_bytes();
        // Set first bit of first byte to 0
        first_bytes[6] &= 0x7F;
        // Get last 15 bits
        let mut last_bytes = (num & 0x7FFF).to_be_bytes();
        // Set first bit to 1
        last_bytes[6] |= 0x80;
        // Combine bytes
        return vec![first_bytes[6], first_bytes[7], last_bytes[6], last_bytes[7]];
    } else if num < 0x200000000000 {
        // If num is smaller than 2^45, return six bytes
        // Convert to big endian
        // Get bytes by 15 bits
        let mut first_bytes = (num >> 30).to_be_bytes();
        // Set first bit of first byte to 0
        first_bytes[6] &= 0x7F;
        // Get bytes by 15 bits
        let mut second_bytes = ((num >> 15) & 0x7FFF).to_be_bytes();
        // Set first bit of first byte to 0
        second_bytes[6] &= 0x7F;
        // Get last 15 bits
        let mut last_bytes = (num & 0x7FFF).to_be_bytes();
        // Set first bit to 1
        last_bytes[6] |= 0x80;
        // Combine bytes
        return vec![first_bytes[6], first_bytes[7], second_bytes[6], second_bytes[7], last_bytes[6], last_bytes[7]];
    } else {
        // If num is larger than 2^45, return eight bytes
        // Convert to big endian
        // Get bytes by 15 bits
        let mut first_bytes = (num >> 45).to_be_bytes();
        // Set first bit of first byte to 0
        first_bytes[6] &= 0x7F;
        // Get bytes by 15 bits
        let mut second_bytes = ((num >> 30) & 0x7FFF).to_be_bytes();
        // Set first bit of first byte to 0
        second_bytes[6] &= 0x7F;
        // Get bytes by 15 bits
        let mut third_bytes = ((num >> 15) & 0x7FFF).to_be_bytes();
        // Set first bit of first byte to 0
        third_bytes[6] &= 0x7F;
        // Get last 15 bits
        let mut last_bytes = (num & 0x7FFF).to_be_bytes();
        // Set first bit to 1
        last_bytes[6] |= 0x80;
        // Combine bytes
        return vec![first_bytes[6], first_bytes[7], second_bytes[6], second_bytes[7], third_bytes[6], third_bytes[7], last_bytes[6], last_bytes[7]];
    }
}

pub fn as_u64_from_delta_encoding_bytes(bytes: &[u8], cursor: &mut usize) -> u64 {
    // If first bit of first byte is 1, return two bytes
    if bytes[*cursor] & 0x80 == 0x80 {
        // Remove first bit
        let mut new_bytes = [0u8; 8];
        new_bytes[6] = bytes[*cursor];
        new_bytes[7] = bytes[*cursor + 1];
        new_bytes[6] &= 0x7F;
        *cursor += 2;
        return u64::from_be_bytes(new_bytes);
    } else if bytes[*cursor + 2] & 0x80 == 0x80 {
        // If first bit of first byte is 0 and first bit of third byte is 1, return four bytes
        // Convert to big endian
        // Get bytes by 15 bits
        let mut first_bytes = [0u8; 8];
        first_bytes[6] = bytes[*cursor];
        first_bytes[7] = bytes[*cursor + 1];
        // Set first bit of first byte to 0
        first_bytes[6] &= 0x7F;
        // Get last 15 bits
        let mut last_bytes = [0u8; 8];
        last_bytes[6] = bytes[*cursor + 2];
        last_bytes[7] = bytes[*cursor + 3];
        // Set first bit to 0
        last_bytes[6] &= 0x7F;
        // Combine bytes
        *cursor += 4;
        return u64::from_be_bytes(first_bytes) << 15 | u64::from_be_bytes(last_bytes);
    } else if bytes[*cursor + 4] & 0x80 == 0x80 {
        // If first bit of first byte and third byte is 0 and first bit of fifth byte is 1, return six bytes
        // Convert to big endian
        // Get bytes by 15 bits
        let mut first_bytes = [0u8; 8];
        first_bytes[6] = bytes[*cursor];
        first_bytes[7] = bytes[*cursor + 1];
        // Set first bit of first byte to 0
        first_bytes[6] &= 0x7F;
        // Get bytes by 15 bits
        let mut second_bytes = [0u8; 8];
        second_bytes[6] = bytes[*cursor + 2];
        second_bytes[7] = bytes[*cursor + 3];
        // Set first bit of first byte to 0
        second_bytes[6] &= 0x7F;
        // Get last 15 bits
        let mut last_bytes = [0u8; 8];
        last_bytes[6] = bytes[*cursor + 4];
        last_bytes[7] = bytes[*cursor + 5];
        // Set first bit to 0
        last_bytes[6] &= 0x7F;
        // Combine bytes
        *cursor += 6;
        return u64::from_be_bytes(first_bytes) << 30 | u64::from_be_bytes(second_bytes) << 15 | u64::from_be_bytes(last_bytes);
    } else {
        // If first bit of first byte, third byte and fifth byte is 0, return eight bytes
        // Convert to big endian
        // Get bytes by 15 bits
        let mut first_bytes = [0u8; 8];
        first_bytes[6] = bytes[*cursor];
        first_bytes[7] = bytes[*cursor + 1];
        // Set first bit of first byte to 0
        first_bytes[6] &= 0x7F;
        // Get bytes by 15 bits
        let mut second_bytes = [0u8; 8];
        second_bytes[6] = bytes[*cursor + 2];
        second_bytes[7] = bytes[*cursor + 3];
        // Set first bit of first byte to 0
        second_bytes[6] &= 0x7F;
        // Get bytes by 15 bits
        let mut third_bytes = [0u8; 8];
        third_bytes[6] = bytes[*cursor + 4];
        third_bytes[7] = bytes[*cursor + 5];
        // Set first bit of first byte to 0
        third_bytes[6] &= 0x7F;
        // Get last 15 bits
        let mut last_bytes = [0u8; 8];
        last_bytes[6] = bytes[*cursor + 6];
        last_bytes[7] = bytes[*cursor + 7];
        // Set first bit to 0
        last_bytes[6] &= 0x7F;
        // Combine bytes
        *cursor += 8;
        return u64::from_be_bytes(first_bytes) << 45 | u64::from_be_bytes(second_bytes) << 30 | u64::from_be_bytes(third_bytes) << 15 | u64::from_be_bytes(last_bytes);
    }
}



pub fn write_hashmap_to_file(path: &str, map: &HashMap<u64, Vec<u64>>) -> Result<(), Error> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(path)?;
    let total_size: u64 = map.iter().map(|(&key, values)| size_of::<u64>() + values.len() * size_of::<u64>())
        .sum::<usize>() as u64;
    file.set_len(total_size as u64)?;
    let mut mmap = unsafe { MmapMut::map_mut(&file)? };

    // let mut mmap = vec![0u8; total_size as usize];
    
    let mut cursor: usize = 0;
    for (key, value) in map {
        // Write key
        let mut key_bytes = delta_encoding_bytes_from_u64(*key);
        mmap[cursor..cursor + key_bytes.len()].copy_from_slice(&key_bytes);
        cursor += key_bytes.len();
        // Write size of value as u64
        let mut value_size_bytes = delta_encoding_bytes_from_u64(value.len() as u64);
        mmap[cursor..cursor + value_size_bytes.len()].copy_from_slice(&value_size_bytes);
        cursor += value_size_bytes.len();
        // Write value
        let mut prev_num = 0u64;
        for num in value { // Value should be sorted ascending
            let diff = *num - prev_num;
            let mut num_bytes = delta_encoding_bytes_from_u64(diff);
            mmap[cursor..cursor + num_bytes.len()].copy_from_slice(&num_bytes);
            prev_num = *num;
            cursor += num_bytes.len();
        }
    }
    // Resize to actual size and add end marker
    // file.write_all(&mmap[..cursor])?;
    file.set_len(cursor as u64)?;
    mmap.flush()?;
    Ok(())
}

pub fn read_hashmap_from_file(path: &str) -> Result<HashMap<u64, Vec<u64>>, Error> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let mut cursor: usize = 0;
    let mut map = HashMap::new();
    while cursor < mmap.len() {
        // Read key
        let key = as_u64_from_delta_encoding_bytes(&mmap, &mut cursor);
        // Read size of value as u64
        let value_size = as_u64_from_delta_encoding_bytes(&mmap, &mut cursor) as usize;
        // Read value
        let mut value = Vec::new();
        let mut prev_num = 0u64;
        for _ in 0..value_size {
            let num = as_u64_from_delta_encoding_bytes(&mmap, &mut cursor) + prev_num;
            value.push(num);
            prev_num = num;
        }
        map.insert(key, value);
    }
    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delta_encoding_bytes_from_u64() {
        // Test encoding of a number smaller than 2^15
        let num = 0x7FFF;
        let encoded = delta_encoding_bytes_from_u64(num);
        assert_eq!(encoded, vec![0xFF, 0xFF]);
        
        let num = 0x0001;
        let encoded = delta_encoding_bytes_from_u64(num);
        assert_eq!(encoded, vec![0x80, 0x01]);
        
        // Test encoding of a number larger than 2^15 and smaller than 2^30
        let num = 0x00100000;
        let encoded = delta_encoding_bytes_from_u64(num);
        assert_eq!(encoded, vec![0x00, 0x20, 0x80, 0x00]);
    }

    #[test]
    fn test_as_u64_from_delta_encoding_bytes() {
        // Test decoding of a number smaller than 2^15
        let bytes = vec![0xFF, 0xFF];
        let decoded = as_u64_from_delta_encoding_bytes(&bytes, &mut 0);
        assert_eq!(decoded, 0x7FFF);
        
        let bytes = vec![0x80, 0x01];
        let decoded = as_u64_from_delta_encoding_bytes(&bytes, &mut 0);
        assert_eq!(decoded, 0x0001);
        // Test decoding of a number larger than 2^15 and smaller than 2^30
        let bytes = vec![0x00, 0x20, 0x80, 0x00];
        let decoded = as_u64_from_delta_encoding_bytes(&bytes, &mut 0);
        assert_eq!(decoded, 0x00100000);
    }

    #[test]
    fn test_write_hashmap_to_file() {
        let mut map = HashMap::new();
        map.insert(0, vec![1, 2, 3]);
        map.insert(1, vec![4, 5, 6]);
        map.insert(2, vec![7, 8, 9]);
        write_hashmap_to_file("data/test_write_hashmap_to_file", &map).unwrap();
    }
    
    #[test]
    fn test_read_hashmap_from_file() {
        let map = read_hashmap_from_file("data/test_write_hashmap_to_file").unwrap();
        assert_eq!(map.get(&0).unwrap(), &vec![1, 2, 3]);
        assert_eq!(map.get(&1).unwrap(), &vec![4, 5, 6]);
        assert_eq!(map.get(&2).unwrap(), &vec![7, 8, 9]);
    }
    #[test]
    fn test_hashmap_io_with_delta_encoding() {
        // Make a huge hashmap to test
        let mut huge_map: HashMap<u64, Vec<u64>> = HashMap::new();
        for i in 0..100 {
            huge_map.insert(i, vec![i, i + 1, i + 100000, i + 104001, i + 10000000000]);
        }
        // Write to file
        let start = Instant::now();
        write_hashmap_to_file("data/huge_map", &huge_map).unwrap();
        println!("Time to write: {:?}", start.elapsed());
        // Drop huge_map
        mem::drop(huge_map);
        // Read from file
        let start = Instant::now();
        let map = read_hashmap_from_file("data/huge_map").unwrap();
        println!("Time to read: {:?}", start.elapsed());
        
        // Check if map is correct
        for i in 0..100 {
            assert_eq!(map.get(&i).unwrap(), &vec![i, i + 1, i + 100000, i + 104001, i + 10000000000]);
        }
    }
}