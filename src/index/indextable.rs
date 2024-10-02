

use std::cell::UnsafeCell;
use std::io::Write;
use std::mem::ManuallyDrop;
use std::sync::atomic::{AtomicUsize, Ordering};

use memmap2::{Mmap, MmapMut};
// use rayon::iter::{IndexedParallelIterator, ParallelIterator};


pub struct FolddiscoIndex {
    offsets: UnsafeCell<Vec<usize>>,
    last_id: UnsafeCell<Vec<usize>>,
    loaded_offsets: ManuallyDrop<Vec<usize>>,
    total_hashes: usize,
    entries: UnsafeCell<MmapMut>,
    index_path: String,
}

unsafe impl Sync for FolddiscoIndex {}

impl FolddiscoIndex {
    pub fn new(total_hashes: usize, path: String) -> Self {
        let offsets = vec![0usize; total_hashes + 1];
        let last_id = vec![usize::MAX; total_hashes];
        let entries = MmapMut::map_anon(1024).unwrap();

        FolddiscoIndex {
            offsets: UnsafeCell::new(offsets),
            // atomic_offsets: UnsafeCell::new(atomic_offsets),
            last_id: UnsafeCell::new(last_id),
            loaded_offsets: ManuallyDrop::new(vec![]),
            total_hashes,
            entries: UnsafeCell::new(entries),
            index_path: path,
        }
    }

    #[inline(always)]
    pub fn count_single_entry(&self, hash: u32, id: usize) {
        let last_id = unsafe { &mut *self.last_id.get() };
        // let atomic_offsets = unsafe { &mut *self.atomic_offsets.get() };
        let offsets = unsafe { &mut *self.offsets.get() };
        let count = if last_id[hash as usize] == usize::MAX {
            if id == 0 {
                1usize
            } else {
                (1 + id.ilog2() / 7) as usize
            }
        } else {
            (1 + (id - last_id[hash as usize]).ilog2() / 7) as usize
        };
        let atomic_offset = unsafe { AtomicUsize::from_ptr( &mut offsets[hash as usize] ) };
        atomic_offset.fetch_add(count, Ordering::Relaxed);
            // atomic_offsets[hash as usize].fetch_add(count, Ordering::SeqCst);
        last_id[hash as usize] = id;
    }
    
    pub fn count_entries(&self, hashes: &Vec<u32>, id: usize) {
        let last_id = unsafe { &mut *self.last_id.get() };
        // let atomic_offsets = unsafe { &mut *self.atomic_offsets.get() };
        let offsets = unsafe { &mut *self.offsets.get() };
        let id_log2 = if id == 0 { 0 } else { id.ilog2() };

        for &hash in hashes {
            let count = if last_id[hash as usize] == usize::MAX {
                if id == 0 {
                    1usize
                } else {
                    (1 + id_log2 / 7) as usize
                }
            } else {
                (1 + (id - last_id[hash as usize]).ilog2() / 7) as usize
            };

            // atomic_offsets[hash as usize].fetch_add(count, Ordering::SeqCst);
            let atomic_offset = unsafe { AtomicUsize::from_ptr( &mut offsets[hash as usize] ) };
            atomic_offset.fetch_add(count, Ordering::SeqCst);
            last_id[hash as usize] = id;
        }
    }

    pub fn add_entries(&self, hashes: &[u32], id: usize, bit_container: &mut Vec<u8>) {
        let last_id = unsafe { &mut *self.last_id.get() };
        // let atomic_offsets = unsafe { &mut *self.atomic_offsets.get() };
        let offsets = unsafe { &mut *self.offsets.get() };
        let entries = unsafe { &mut *self.entries.get() };

        for &hash in hashes {

            let count = if last_id[hash as usize] == usize::MAX {
                if id == 0 {
                    1usize
                } else {
                    (1 + id.ilog2() / 7) as usize
                }
            } else {
                (1 + (id - last_id[hash as usize]).ilog2() / 7) as usize
            };

            let atomic_offset = unsafe { AtomicUsize::from_ptr( &mut offsets[hash as usize] ) };
            let offset = atomic_offset.fetch_add(count, Ordering::SeqCst);
            
            
            // let offset = atomic_offsets[hash as usize].fetch_add(count, Ordering::SeqCst);
            let prev = last_id[hash as usize];
            last_id[hash as usize] = id;

            // TODO: write entry compressed
            // Split id bytes into 7-bit chunks
            let id_to_split: usize = match prev {
                usize::MAX => id,
                _ => id - prev,
            };
            let nbit = split_by_seven_bits(id_to_split, bit_container);
            for i in 0..nbit {
                let bit = bit_container[i];
                entries[offset + i] = bit;
            }
        }
    }

    #[inline(always)]
    pub fn add_single_entry(&self, hash: u32, id: usize, bit_container: &mut Vec<u8>) {
        let last_id = unsafe { &mut *self.last_id.get() };
        // let atomic_offsets = unsafe { &mut *self.atomic_offsets.get() };
        let offsets = unsafe { &mut *self.offsets.get() };
        let entries = unsafe { &mut *self.entries.get() };

        let count = if last_id[hash as usize] == usize::MAX {
            if id == 0 {
                1usize
            } else {
                (1 + id.ilog2() / 7) as usize
            }
        } else {
            (1 + (id - last_id[hash as usize]).ilog2() / 7) as usize
        };

        // let offset = atomic_offsets[hash as usize].fetch_add(count, Ordering::SeqCst);
        let atomic_offset = unsafe { AtomicUsize::from_ptr( &mut offsets[hash as usize] ) };
        let offset = atomic_offset.fetch_add(count, Ordering::Relaxed);
        let prev = last_id[hash as usize];
        last_id[hash as usize] = id;
        let id_to_split: usize = match prev {
            usize::MAX => id,
            _ => id - prev,
        };
            
        let nbit = split_by_seven_bits(id_to_split, bit_container);
        for i in 0..nbit {
            let bit = bit_container[i];
            entries[offset + i] = bit;
        }
    }
    
    pub fn allocate_entries(&self) {
        let offsets = unsafe { &mut *self.offsets.get() };
        let last_id = unsafe { &mut *self.last_id.get() };
        let entries = unsafe { &mut *self.entries.get() };

        let mut total_entries = 0;
        for i in 0..self.total_hashes {
            let current_entry_size = offsets[i];
            offsets[i] = total_entries;
            total_entries += current_entry_size;
        }
        let index_path = format!("{}.value", self.index_path);
        let index_file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .open(index_path)
            .unwrap();
        index_file.set_len(total_entries as u64).unwrap();
        let mmap = unsafe { MmapMut::map_mut(&index_file).unwrap() };
        *entries = mmap;
        offsets[self.total_hashes] = total_entries;

        for id in last_id.iter_mut() {
            *id = usize::MAX;
        }
    }

    pub fn finish_index(&self) {
        let offsets = unsafe { &mut *self.offsets.get() };

        for i in (1..=self.total_hashes).rev() {
            offsets[i] = offsets[i - 1];
        }
        offsets[0] = 0;
    }

    pub fn get_raw_entries(&self, hash: usize) -> &[u8] {
        let offsets = if self.loaded_offsets.is_empty() {
            unsafe { &*self.offsets.get() }
        } else {
            &self.loaded_offsets
        };
        let entries = unsafe { &*self.entries.get() };

        let start = offsets[hash];
        let end = offsets[hash + 1];

        &entries[start..end]
    }
    
    pub fn get_entries(&self, hash: u32) -> Vec<usize> {
        let raw_entries = self.get_raw_entries(hash as usize);
        merge_usize_vec_from_bytes(raw_entries)
    }
    
    pub fn save_offset_to_file(&self) {
        let offsets = unsafe { &*self.offsets.get() };
        let offset_path = format!("{}.offset", self.index_path);
        let file = std::fs::File::create(&offset_path).expect("Unable to create offset file");
        let mut writer = std::io::BufWriter::new(file);
        let offset_bytes = unsafe {
            std::slice::from_raw_parts(offsets.as_ptr() as *const u8,
            offsets.len() * std::mem::size_of::<usize>())
        };
        writer.write_all(offset_bytes).expect("Unable to write offset bytes");
    }
        
}

pub fn load_big_index(index_prefix: &str) -> (FolddiscoIndex, Mmap) {
    let offset_path = format!("{}.offset", index_prefix);
    let index_path = format!("{}.value", index_prefix);
    let offset_file = std::fs::File::open(&offset_path).expect("Unable to open offset file");
    let offset_mmap = unsafe { Mmap::map(&offset_file).expect("Unable to map offset file") };
    let offsets = unsafe {
        let offsets_ptr = offset_mmap.as_ptr() as *const usize;
        ManuallyDrop::new(Vec::from_raw_parts(
            offsets_ptr as *mut usize, offset_mmap.len() / std::mem::size_of::<usize>(),
            offset_mmap.len() / std::mem::size_of::<usize>()
        ))
    };
    let total_hashes = offsets.len() - 1;

    let entries_file = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .open(&index_path)
        .expect("Unable to open index file");
    let entries_mmap = unsafe { MmapMut::map_mut(&entries_file).expect("Unable to map index file") };

    ( FolddiscoIndex {
        offsets: UnsafeCell::new(vec![]),
        last_id: UnsafeCell::new(vec![]),
        loaded_offsets: offsets,
        total_hashes: total_hashes,
        entries: UnsafeCell::new(entries_mmap),
        index_path,
    }, offset_mmap )
}


// fn split_by_seven_bits(mut id: usize) -> Vec<u8> {
//     let mut result = vec![];
//     while id > 0 {
//         let byte = (id & 0x7F) as u8;
//         id >>= 7;
//         if id > 0 {
//             result.push(byte | 0x80); // Set continuation bit
//         } else {
//             result.push(byte); // Last byte, no continuation
//         }
//     }

//     if result.is_empty() {
//         result.push(0); // Ensure at least one byte for id = 0
//     }

//     result
// }

#[inline(always)]
fn split_by_seven_bits(mut id: usize, bit_container: &mut Vec<u8>) -> usize {
    let mut length = 0usize;
    bit_container.clear();
    while id > 0 {
        let byte = (id & 0x7F) as u8;
        id >>= 7;
        if id > 0 {
            bit_container.push(byte | 0x80); // Set continuation bit
            length += 1;
        } else {
            bit_container.push(byte); // Last byte, no continuation
            length += 1;
        }
    }

    if bit_container.is_empty() {
        bit_container.push(0); // Ensure at least one byte for id = 0
        length += 1;
    }

    length
}

#[inline(always)]
fn merge_seven_bits(bytes: &[u8]) -> usize {
    let mut result = 0;
    let mut i = bytes.len() - 1;
    let mut shift = 7 * i;
    loop {
        let byte = bytes[i];
        result |= ((byte & 0x7F) as usize) << shift;
        // If first bit is set to 1, break
        if i == 0 {
            break;
        }
        shift -= 7;
        i -= 1;
    }
    result
}

#[inline(always)]
fn merge_usize_vec_from_bytes(bytes: &[u8]) -> Vec<usize> {
    let mut result = vec![];
    let mut start: usize = 0;
    let mut end: usize = 0;
    let mut prev: usize = usize::MAX;
    while end < bytes.len() {
        if bytes[end] & 0x80 != 0 {
            end += 1;
            continue;
        } else {
            // If first bit is set to 1, merge bytes
            let merged = merge_seven_bits(&bytes[start..=end]);
            if prev != usize::MAX {
                prev += merged;
                result.push(prev);
            } else {
                prev = merged;
                result.push(prev);
            }
            start = end + 1;
            end = start;
        }
    }
    result
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_folddisco_index() {
        let total_hashes = 10;
        let index = FolddiscoIndex::new(total_hashes, "test.index".to_string());

        let hashes1: Vec<u32> = (0u32..7).collect();
        let id1 = 0usize;
        let hashes4: Vec<u32> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let id4 = 10usize;
        let hashes3: Vec<u32> = vec![2, 4, 6, 8];
        let id3 = 128usize;
        let hashes2: Vec<u32>= vec![1, 3, 5, 7, 9];
        let id2 = 655345usize;

        index.count_entries(&hashes1, id1);
        index.count_entries(&hashes4, id4);
        index.count_entries(&hashes3, id3);
        index.count_entries(&hashes2, id2);

        index.allocate_entries();
        let mut bit_container = Vec::with_capacity(8);
        index.add_entries(&hashes1, id1, &mut bit_container);
        index.add_entries(&hashes4, id4, &mut bit_container);
        index.add_entries(&hashes3, id3, &mut bit_container);
        index.add_entries(&hashes2, id2, &mut bit_container);
        
        index.finish_index();

        // Print offsets
        let offsets = unsafe { &*index.offsets.get() };
        for i in 0..total_hashes+1 {
            println!("offsets[{}]: {}", i, offsets[i]);
        }
        let entries = unsafe { &*index.entries.get() };
        println!("{:?}", entries);
        
        let entries1 = index.get_raw_entries(7);
        println!("{:?}, {}", entries1, entries1.len());
        
        println!("{:?}", merge_usize_vec_from_bytes(entries1));
        
        for i in 0..total_hashes {
            let entries = index.get_entries(i as u32);
            println!("Entries for hash {}: {:?}", i, entries);
        }
    }
    
}

