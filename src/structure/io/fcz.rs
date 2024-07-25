// TODO: Need conversion from atom_t to CompactStructure

use libc;
use memmap2::Mmap;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{fs::File, io::{BufRead, Read}};

// include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
include!("../../../lib/foldcomp/bindings.rs");

pub fn read_foldcomp_db_lookup(db_path: &str) -> Result<Vec<(usize, String)>, &'static str> {
    // Check if the file exists
    let lookup_path = format!("{}.lookup", db_path);
    let mut lookup_file = match File::open(&lookup_path) {
        Ok(file) => file,
        Err(_) => return Err("Lookup file not found."),
    };
    let mut output: Vec<(usize, String)> = Vec::new();
    let reader = std::io::BufReader::new(lookup_file);
    for line in reader.lines() {
        let line = line.unwrap();
        let mut split = line.split_whitespace();
        let id = split.next().unwrap().parse::<usize>().unwrap();
        let name = split.next().unwrap().to_string();
        output.push((id, name));
    }
    
    // Return
    Ok(output)
}

pub fn read_foldcomp_db_index(db_path: &str) -> Result<Vec<(usize, usize, usize)>, &'static str> {
    // Check if the file exists
    let index_path = format!("{}.index", db_path);
    let index_file = match File::open(&index_path) {
        Ok(file) => file,
        Err(_) => return Err("Index file not found."),
    };
    let mut output: Vec<(usize, usize, usize)> = Vec::new();
    let reader = std::io::BufReader::new(index_file);
    for line in reader.lines() {
        let line = line.unwrap();
        let mut split = line.split_whitespace();
        let id = split.next().unwrap().parse::<usize>().unwrap();
        let start = split.next().unwrap().parse::<usize>().unwrap();
        let length = split.next().unwrap().parse::<usize>().unwrap();
        output.push((id, start, length));
    }
    Ok(output)
}

pub fn get_path_vector_out_of_lookup(lookup: &Vec<(usize, String)>) -> Vec<String> {
    let mut output: Vec<String> = Vec::new();
    for (_, name) in lookup {
        output.push(name.clone());
    }
    output
}

pub fn read_foldcomp_db(db_path: &str) -> Result<(Mmap, &[u8]), &'static str> {
    // Check if the file exists
    let mut db_file = match File::open(&db_path) {
        Ok(file) => file,
        Err(_) => return Err("DB file not found."),
    };
    
    // Read the file
    let mmap = unsafe { Mmap::map(&db_file).unwrap() };
    let db = unsafe { std::slice::from_raw_parts(mmap.as_ptr(), mmap.len()) };
    
    Ok((mmap, db))
}

pub fn get_foldcomp_db_entry<'a>(db: &'a [u8], index: &(usize, usize, usize)) -> &'a [u8] {
    &db[index.1..index.1 + index.2]
}

pub fn get_foldcomp_db_entry_by_id<'a>(db: &'a [u8], index_vector: &Vec<(usize, usize, usize)>, id: usize) -> Option<&'a [u8]> {
    let entry_index = index_vector.binary_search_by_key(&id, |&(id, _, _)| id);
    let entry: &(usize, usize, usize) = match entry_index {
        Ok(index) => index_vector.get(index).unwrap(),
        Err(_) => return None,
    };
    Some(get_foldcomp_db_entry(db, entry))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_foldcomp() {
        unsafe {
            // Test single FCZ file
            let instance = foldcomp_create();
            // Read data/7m0y.fcz as binary and pass it.
            let mut input = File::open("data/foldcomp/7m0y.fcz").unwrap();
            let mut content : Vec<u8> = Vec::new();
            input.read_to_end(&mut content).unwrap();
            let mut atom_count = libc::size_t::default();
            let output_ptr = foldcomp_process(instance, content.as_ptr(), content.len() as libc::size_t, &mut atom_count);
            let output: &[atom_t] = std::slice::from_raw_parts(output_ptr, atom_count as usize);
            println!("Total {} atoms", atom_count);
            println!("First element: {:?}", output[0]);
            println!("Last element: {:?}", output[atom_count as usize - 1]);
            foldcomp_destroy(instance);
            foldcomp_free(output_ptr);
        }
        println!("Single FCZ file test passed.");
        unsafe {
            let db_path = "data/foldcomp/example_db";
            // let db_path = "data/s_cerevisiae";
            let (db_mmap, db) = read_foldcomp_db(db_path).unwrap();
            let lookup = read_foldcomp_db_lookup(db_path).unwrap();
            let index = read_foldcomp_db_index(db_path).unwrap();
            let path_vector = get_path_vector_out_of_lookup(&lookup);
            
            // Test single entry
            let path1 = &path_vector[0];
            println!("Path: {}", path1);
            let entry1 = get_foldcomp_db_entry_by_id(db, &index, 0).unwrap();
            let instance = foldcomp_create();
            let mut atom_count = libc::size_t::default();
            let output_ptr = foldcomp_process(instance, entry1.as_ptr(), entry1.len() as libc::size_t, &mut atom_count);
            let output: &[atom_t] = std::slice::from_raw_parts(output_ptr, atom_count as usize);
            println!("Total {} atoms", atom_count);
            println!("First element: {:?}", output[0]);
            println!("Last element: {:?}", output[atom_count as usize - 1]);
            foldcomp_destroy(instance);
            foldcomp_free(output_ptr);
            println!("Single DB entry test passed.");
            // Iterate over all entries
            let _ = &index.par_iter().for_each(|entry_index| {
                let entry = get_foldcomp_db_entry(db, entry_index);
                let instance = foldcomp_create();
                let mut atom_count = libc::size_t::default();
                let output_ptr = foldcomp_process(instance, entry.as_ptr(), entry.len() as libc::size_t, &mut atom_count);
                let output: &[atom_t] = std::slice::from_raw_parts(output_ptr, atom_count as usize);
                println!("Total {} atoms", atom_count);
                println!("First element: {:?}", output[0]);
                println!("Last element: {:?}", output[atom_count as usize - 1]);
                foldcomp_destroy(instance);
                foldcomp_free(output_ptr);
            });
            println!("Full DB entry test passed.");
        }
    }
}
// Test passed - 2024-07-25 22:12:55