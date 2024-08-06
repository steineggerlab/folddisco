// DONE: Need conversion from &[atom_t] to Structure
// TODO: 
// Us
// include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
include!("../../../lib/foldcomp/bindings.rs");

use libc;
use memmap2::{Mmap, MmapMut};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::prelude::ParallelSliceMut;
use std::{fs::File, io::{BufRead, Read}};
use std::mem::ManuallyDrop;

use crate::structure::atom::Atom;
use crate::structure::core::Structure;
use crate::structure::io::StructureFileFormat;

/// A FCZ DB reader
#[derive(Debug)]
pub struct FoldcompDbReader {
    /// The underlying reader
    pub path: String,
    pub input_type: StructureFileFormat,
    pub db_mmap: Mmap,
    pub db: ManuallyDrop<Vec<u8>>,
    pub lookup: Vec<(usize, String)>,
    pub index: Vec<(usize, usize, usize)>,
}

impl Drop for FoldcompDbReader {
    fn drop(&mut self) {
        unsafe {
            std::mem::forget(ManuallyDrop::take(&mut self.db));
        }
    }
}

impl FoldcompDbReader {
    pub fn new(path: &str) -> Self {
        let (db_mmap, db) = read_foldcomp_db(path).expect("Error reading foldcomp db file.");
        let lookup = read_foldcomp_db_lookup(path).expect("Error reading foldcomp db lookup file.");
        let index = read_foldcomp_db_index(path).expect("Error reading foldcomp db index file.");
        let path_string_to_return = path.to_string();
        
        let mut lookup = lookup;
        lookup.par_sort_unstable_by(|a, b| a.1.cmp(&b.1));
        
        FoldcompDbReader {
            path: path_string_to_return,
            input_type: StructureFileFormat::FCZDB,
            db_mmap: db_mmap,
            db: db,
            lookup: lookup,
            index: index,
        }
    }
    
    pub fn empty() -> Self {
        FoldcompDbReader {
            path: String::new(),
            input_type: StructureFileFormat::FCZDB,
            db_mmap: MmapMut::map_anon(1).unwrap().make_read_only().unwrap(),
            db: ManuallyDrop::new(Vec::new()),
            lookup: Vec::new(),
            index: Vec::new(),
        }
    }

    pub fn read_single_structure(&self, name: &str) -> Result<Structure, String> {
        let mut structure = Structure::new(); // revise
        let mut record = (b' ', 0);
        let entry = get_foldcomp_db_entry_by_name(&self.db, &self.lookup, &self.index, name);
        match entry {
            Some(entry) => unsafe {
                let instance = foldcomp_create();
                let mut atom_count = libc::size_t::default();
                let output_ptr = foldcomp_process(instance, entry.as_ptr(), entry.len() as libc::size_t, &mut atom_count);
                let output: &[atom_t] = std::slice::from_raw_parts(output_ptr, atom_count as usize);
                for atom in output {
                    let atom = Atom::from_c(atom);
                    structure.update(atom.clone(), &mut record);
                }
                foldcomp_destroy(instance);
                foldcomp_free(output_ptr);
                Ok(structure)
            }
            None => Err(format!("Entry with name {} not found.", name)),
        }
    }
    
    pub fn read_single_structure_by_id(&self, id: usize) -> Result<Structure, String> {
        let mut structure = Structure::new(); // revise
        let mut record = (b' ', 0);
        let entry = get_foldcomp_db_entry_by_id(&self.db, &self.index, id);
        match entry {
            Some(entry) => unsafe {
                let instance = foldcomp_create();
                let mut atom_count = libc::size_t::default();
                let output_ptr = foldcomp_process(instance, entry.as_ptr(), entry.len() as libc::size_t, &mut atom_count);
                let output: &[atom_t] = std::slice::from_raw_parts(output_ptr, atom_count as usize);
                for atom in output {
                    let atom = Atom::from_c(atom);
                    structure.update(atom.clone(), &mut record);
                }
                foldcomp_destroy(instance);
                foldcomp_free(output_ptr);
                Ok(structure)
            }
            None => Err(format!("Entry with ID {} not found.", id)),
        }
    }
    
    pub fn get_paths(&self) -> Vec<String> {
        get_path_vector_out_of_lookup_and_index(&self.lookup, &self.index)
    }
    
    pub fn sort_lookup_by_id(&mut self) {
        self.lookup.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
    }
    
    pub fn sort_lookup_by_name(&mut self) {
        self.lookup.par_sort_unstable_by(|a, b| a.1.cmp(&b.1));
    }
}

// Methods to convert atom_t to Atom
impl Atom {
    pub fn from_c(atom: &atom_t) -> &Self {
        unsafe { &*(atom as *const atom_t as *const Atom) }
    }
    
    pub fn from_c_mut(atom: &mut atom_t) -> &mut Self {
        unsafe { &mut *(atom as *mut atom_t as *mut Atom) }
    }

    pub fn as_c(&self) -> &atom_t {
        unsafe { &*(self as *const Atom as *const atom_t) }
    }

    pub fn as_c_mut(&mut self) -> &mut atom_t {
        unsafe { &mut *(self as *mut Atom as *mut atom_t) }
    }
}

// Convert atom_t slice to Structure
pub unsafe fn atom_t_slice_to_structure(slice: &[atom_t]) -> Structure {
    let mut structure = Structure::new(); 
    let mut record = (b' ', 0);
    for atom in slice {
        let atom = Atom::from_c(atom);
        structure.update(atom.clone(), &mut record);
    }
    structure
}

// Functions to read foldcomp db files; db, lookup, index
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

pub fn get_path_vector_out_of_lookup_and_index(lookup: &Vec<(usize, String)>, index: &Vec<(usize, usize, usize)>) -> Vec<String> {
    let mut output: Vec<String> = Vec::new();
    for (id, _, _) in index {
        let entry_index = lookup.binary_search_by_key(id, |(id, _)| *id);
        let entry: &(usize, String) = match entry_index {
            Ok(index) => lookup.get(index).unwrap(),
            Err(_) => continue,
        };
        output.push(entry.1.clone());
    }
    output
}

pub fn get_id_vector_out_of_lookup(lookup: &Vec<(usize, String)>) -> Vec<usize> {
    let mut output: Vec<usize> = Vec::new();
    for (id, _) in lookup {
        output.push(*id);
    }
    output
}

pub fn get_id_vector_subset_out_of_lookup(lookup: &Vec<(usize, String)>, subset_names: &Vec<String>) -> Vec<usize> {
    let mut output: Vec<usize> = Vec::new();
    // Order should be maintained
    for name in subset_names {
        let entry_index = lookup.binary_search_by_key(name, |(_, name)| name.to_string());
        let entry: &(usize, String) = match entry_index {
            Ok(index) => lookup.get(index).unwrap(),
            Err(_) => continue,
        };
        output.push(entry.0);
    }
    output
}

pub fn get_name_vector_subset_out_of_lookup(lookup: &Vec<(usize, String)>, subset_ids: &Vec<usize>) -> Vec<String> {
    let mut output: Vec<String> = Vec::new();
    // Order should be maintained
    for id in subset_ids {
        let entry_index = lookup.binary_search_by_key(id, |(id, _)| *id);
        let entry: &(usize, String) = match entry_index {
            Ok(index) => lookup.get(index).unwrap(),
            Err(_) => continue,
        };
        output.push(entry.1.clone());
    }
    output
}

pub fn read_foldcomp_db(db_path: &str) -> Result<(Mmap, ManuallyDrop<Vec<u8>>), &'static str> {
    // Check if the file exists
    let mut db_file = match File::open(&db_path) {
        Ok(file) => file,
        Err(_) => return Err("DB file not found."),
    };
    
    // Read the file
    let mmap = unsafe { Mmap::map(&db_file).unwrap() };
    let db = unsafe { ManuallyDrop::new(Vec::from_raw_parts(mmap.as_ptr() as *mut u8, mmap.len(), mmap.len())) };
    Ok((mmap, db))
}

pub fn get_foldcomp_db_entry<'a>(db: &'a ManuallyDrop<Vec<u8>>, index: &(usize, usize, usize)) -> &'a [u8] {
    &db[index.1..index.1 + index.2]
}

pub fn get_foldcomp_db_entry_by_id<'a>(db: &'a ManuallyDrop<Vec<u8>>, index_vector: &Vec<(usize, usize, usize)>, id: usize) -> Option<&'a [u8]> {
    let entry_index = index_vector.binary_search_by_key(&id, |&(id, _, _)| id);
    let entry: &(usize, usize, usize) = match entry_index {
        Ok(index) => index_vector.get(index).unwrap(),
        Err(_) => return None,
    };
    Some(get_foldcomp_db_entry(db, entry))
}

pub fn get_foldcomp_db_entry_by_name<'a>(
    db: &'a ManuallyDrop<Vec<u8>>, lookup: &Vec<(usize, String)>, index: &Vec<(usize, usize, usize)>, name: &str
) -> Option<&'a [u8]> {
    let entry_index = lookup.binary_search_by_key(&name, |(_, name)| name);
    let entry: &(usize, String) = match entry_index {
        Ok(index) => lookup.get(index).unwrap(),
        Err(_) => return None,
    };
    get_foldcomp_db_entry_by_id(db, index, entry.0)
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
            let path_vector = get_path_vector_out_of_lookup_and_index(&lookup, &index);
            
            // Test single entry
            let path1 = &path_vector[0];
            println!("Path: {}", path1);
            // let entry1 = get_foldcomp_db_entry_by_id(db, &index, 0).unwrap();
            let entry1 = get_foldcomp_db_entry_by_name(&db, &lookup, &index, path1).unwrap();
            let instance = foldcomp_create();
            let mut atom_count = libc::size_t::default();
            let output_ptr = foldcomp_process(instance, entry1.as_ptr(), entry1.len() as libc::size_t, &mut atom_count);
            let output: &[atom_t] = std::slice::from_raw_parts(output_ptr, atom_count as usize);
            println!("Total {} atoms", atom_count);
            println!("First element: {:?}", output[0]);
            println!("Last element: {:?}", output[atom_count as usize - 1]);
            foldcomp_destroy(instance);
            foldcomp_free(output_ptr);
            let output_structure = atom_t_slice_to_structure(output);
            println!("Structure: {:?}", output_structure);
            let compact = output_structure.to_compact();
            println!("CompactStructure: {:?}", compact);
            println!("Single DB entry test passed.");
            // Iterate over all entries
            let _ = &index.par_iter().for_each(|entry_index| {
                let entry = get_foldcomp_db_entry(&db, entry_index);
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
    
    #[test]
    fn test_foldcomp_db_reader() {
        let db_path = "data/foldcomp/example_db";
        let reader = FoldcompDbReader::new(db_path);
        let path_vector = reader.get_paths();
        let path1 = &path_vector[0];
        let structure = reader.read_single_structure(path1).unwrap();
        let compact = structure.to_compact();
        println!("Compact Structure: {:?}", compact);
        println!("Foldcomp DB reader test passed.");
        
        let compact_structure_vector = path_vector.par_iter().map(|path| {
            let structure = reader.read_single_structure(path).unwrap();
            structure.to_compact()
        }).collect::<Vec<_>>(); 
        for compact_structure in compact_structure_vector {
            println!("Compact Structure: {:?}", compact_structure);
        }
        
        let path_vector_subset = vec![path_vector[0].clone(), path_vector[1].clone()];
        let id_vector_subset = get_id_vector_subset_out_of_lookup(&reader.lookup, &path_vector_subset);
        let compact_structure_vector_subset = id_vector_subset.par_iter().map(|id| {
            let structure = reader.read_single_structure_by_id(*id).unwrap();
            structure.to_compact()
        }).collect::<Vec<_>>();
        for compact_structure in compact_structure_vector_subset {
            println!("Compact Structure: {:?}", compact_structure);
        }
    }
}
// Test passed - 2024-07-25 22:12:55