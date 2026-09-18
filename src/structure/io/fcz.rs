// Foldcomp (FCZ) DB reading through the C bindings, with mapped lookup/index caches.
// include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
#![allow(non_camel_case_types, non_snake_case, non_upper_case_globals)]
include!("../../../lib/foldcomp/bindings.rs");

use libc;
use memmap2::{Mmap, MmapMut};
use rayon::prelude::ParallelSliceMut;
use rayon::prelude::ParallelString;
use rayon::prelude::*;
use std::fs::File;
use std::mem::ManuallyDrop;

use rustc_hash::FxHashMap as HashMap;

use crate::structure::atom::Atom;
use crate::utils::pod_cache::{store_and_map, CacheRecord, PodCache};
use crate::structure::chain_id::ChainId;
use crate::structure::core::Structure;
use crate::structure::io::StructureFileFormat;

// Mapped lookup and index tables.
//
// The MMseqs2-style `.lookup` (`key \t name \t file`) and `.index`
// (`key \t offset \t length`) are parsed and key-sorted once into
// `<db>.lookup.fdcache` / `<db>.index.fdcache`, then mapped on every open.
// Not `<db>.index.cache`: Foldcomp's own `save_cache` uses that name for another format.

/// One `.lookup` entry: a DB key and its name. 24 bytes, no padding.
#[repr(C)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FoldcompLookupRecord {
    pub key: u64,
    pub name_offset: u64,
    pub name_len: u32,
    pub _pad: u32,
}

// SAFETY: `#[repr(C)]` u64, u64, u32, u32: 24 bytes, no padding, any bit pattern valid.
unsafe impl CacheRecord for FoldcompLookupRecord {
    const MAGIC: &'static [u8; 8] = b"FDFCLKUP";
    const VERSION: u32 = 1;
}

/// One `.index` entry: a DB key and the byte range of its entry. 24 bytes.
#[repr(C)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FoldcompIndexRecord {
    pub key: u64,
    pub offset: u64,
    pub length: u64,
}

// SAFETY: `#[repr(C)]` three u64s: 24 bytes, no padding, any bit pattern valid.
unsafe impl CacheRecord for FoldcompIndexRecord {
    const MAGIC: &'static [u8; 8] = b"FDFCIDXT";
    const VERSION: u32 = 1;
}

/// A Foldcomp DB `.lookup`, mapped and sorted by key.
pub struct FoldcompLookup {
    cache: PodCache<FoldcompLookupRecord>,
}

impl FoldcompLookup {
    pub fn empty() -> Self {
        FoldcompLookup { cache: PodCache::empty() }
    }

    /// Map the cache for `<db_path>.lookup`, or build it from the text.
    pub fn load(db_path: &str) -> Result<Self, &'static str> {
        let source = format!("{}.lookup", db_path);
        let cache_path = format!("{}.fdcache", source);
        if let Some(cache) = Self::map(&cache_path, &source) {
            return Ok(FoldcompLookup { cache });
        }
        let file = File::open(&source).map_err(|_| "Lookup file not found.")?;
        let mmap = unsafe { Mmap::map(&file).map_err(|_| "Unable to mmap the lookup file.")? };
        let content = unsafe { std::str::from_utf8_unchecked(&mmap) };

        // (key, name offset within `content`, name length)
        let mut parsed: Vec<(u64, u64, u32)> = content.par_lines().map(|line| {
            let mut split = line.split('\t');
            let key = split.next().unwrap().parse::<u64>().unwrap();
            let name = split.next().unwrap();
            let offset = name.as_ptr() as usize - content.as_ptr() as usize;
            (key, offset as u64, name.len() as u32)
        }).collect();
        // Sorted at build time; per-hit reads binary-search by key.
        parsed.par_sort_unstable_by_key(|entry| entry.0);

        let names_len: usize = parsed.iter().map(|entry| entry.2 as usize).sum();
        let mut names = Vec::with_capacity(names_len);
        let mut records = Vec::with_capacity(parsed.len());
        let source_bytes = content.as_bytes();
        for (key, text_offset, name_len) in parsed {
            let start = text_offset as usize;
            let name_offset = names.len() as u64;
            names.extend_from_slice(&source_bytes[start..start + name_len as usize]);
            records.push(FoldcompLookupRecord { key, name_offset, name_len, _pad: 0 });
        }
        let cache = store_and_map(&cache_path, &source, &records, &names)
            .ok_or("Unable to build the Foldcomp lookup cache.")?;
        Ok(FoldcompLookup { cache })
    }

    /// Map and validate the cache; `None` means rebuild from the text.
    fn map(cache_path: &str, source: &str) -> Option<PodCache<FoldcompLookupRecord>> {
        let cache = PodCache::<FoldcompLookupRecord>::map(cache_path, source)?;
        let names_len = cache.names_blob().len() as u64;
        let records = cache.records();
        // Non-empty names and strictly ascending keys, checked in one pass;
        // together they reject the all-zero records of a torn write.
        let usable = |record: &FoldcompLookupRecord| {
            record.name_len > 0
                && record.name_offset.saturating_add(record.name_len as u64) <= names_len
        };
        let sound = match records.split_last() {
            None => true,
            Some((last, _)) => {
                usable(last)
                    && records.par_windows(2)
                        .all(|pair| usable(&pair[0]) && pair[0].key < pair[1].key)
            }
        };
        if sound { Some(cache) } else { None }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.cache.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }

    #[inline]
    pub fn records(&self) -> &[FoldcompLookupRecord] {
        self.cache.records()
    }

    #[inline]
    pub fn name(&self, index: usize) -> &str {
        let record = &self.records()[index];
        self.cache.name_at(record.name_offset, record.name_len)
    }

    #[inline]
    pub fn key(&self, index: usize) -> usize {
        self.records()[index].key as usize
    }

    pub fn names(&self) -> impl Iterator<Item = &str> + '_ {
        (0..self.len()).map(move |i| self.name(i))
    }

    /// The name a key belongs to, by binary search over the key-sorted records.
    pub fn name_of_key(&self, key: usize) -> Option<&str> {
        let records = self.records();
        let position = records.binary_search_by_key(&(key as u64), |record| record.key).ok()?;
        Some(self.name(position))
    }

    /// The key a name belongs to. A linear scan: names are resolved only per query
    /// (`db:name`), so a name-sorted index would not pay for its size.
    pub fn key_of_name(&self, name: &str) -> Option<usize> {
        self.records().par_iter()
            .find_any(|record| self.cache.name_at(record.name_offset, record.name_len) == name)
            .map(|record| record.key as usize)
    }

    /// Keys of several names in the given order, skipping unknown names. One pass.
    pub fn keys_of_names(&self, names: &[String]) -> Vec<usize> {
        let wanted: HashMap<&str, usize> = names.iter().enumerate()
            .map(|(i, name)| (name.as_str(), i))
            .collect();
        let mut found: Vec<Option<u64>> = vec![None; names.len()];
        for record in self.records() {
            let name = self.cache.name_at(record.name_offset, record.name_len);
            if let Some(&position) = wanted.get(name) {
                found[position] = Some(record.key);
            }
        }
        found.into_iter().flatten().map(|key| key as usize).collect()
    }
}

impl std::fmt::Debug for FoldcompLookup {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "FoldcompLookup({} entries)", self.len())
    }
}

/// A Foldcomp DB `.index`, mapped and sorted by key.
pub struct FoldcompIndex {
    cache: PodCache<FoldcompIndexRecord>,
}

impl FoldcompIndex {
    pub fn empty() -> Self {
        FoldcompIndex { cache: PodCache::empty() }
    }

    /// Map the cache for `<db_path>.index`, or build it from the text.
    pub fn load(db_path: &str) -> Result<Self, &'static str> {
        let source = format!("{}.index", db_path);
        let cache_path = format!("{}.fdcache", source);
        if let Some(cache) = Self::map(&cache_path, &source) {
            return Ok(FoldcompIndex { cache });
        }
        let file = File::open(&source).map_err(|_| "Index file not found.")?;
        let mmap = unsafe { Mmap::map(&file).map_err(|_| "Unable to mmap the index file.")? };
        let content = unsafe { std::str::from_utf8_unchecked(&mmap) };

        let mut records: Vec<FoldcompIndexRecord> = content.par_lines().map(|line| {
            let mut split = line.split('\t');
            let key = split.next().unwrap().parse::<u64>().unwrap();
            let offset = split.next().unwrap().parse::<u64>().unwrap();
            let length = split.next().unwrap().parse::<u64>().unwrap();
            FoldcompIndexRecord { key, offset, length }
        }).collect();
        records.par_sort_unstable_by_key(|record| record.key);

        let cache = store_and_map(&cache_path, &source, &records, &[])
            .ok_or("Unable to build the Foldcomp index cache.")?;
        Ok(FoldcompIndex { cache })
    }

    fn map(cache_path: &str, source: &str) -> Option<PodCache<FoldcompIndexRecord>> {
        let cache = PodCache::<FoldcompIndexRecord>::map(cache_path, source)?;
        let records = cache.records();
        // Non-zero lengths and strictly ascending keys; an all-zero record fails both.
        let sound = match records.split_last() {
            None => true,
            Some((last, _)) => {
                last.length > 0
                    && records.par_windows(2)
                        .all(|pair| pair[0].length > 0 && pair[0].key < pair[1].key)
            }
        };
        if sound { Some(cache) } else { None }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.cache.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }

    #[inline]
    pub fn records(&self) -> &[FoldcompIndexRecord] {
        self.cache.records()
    }

    /// The entry for a DB key, by binary search. This is the per-hit path.
    #[inline]
    pub fn find(&self, key: usize) -> Option<&FoldcompIndexRecord> {
        let records = self.records();
        let position = records.binary_search_by_key(&(key as u64), |record| record.key).ok()?;
        Some(&records[position])
    }

    /// Every key, in ascending order.
    pub fn keys(&self) -> impl Iterator<Item = usize> + '_ {
        self.records().iter().map(|record| record.key as usize)
    }
}

impl std::fmt::Debug for FoldcompIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "FoldcompIndex({} entries)", self.len())
    }
}

/// Memory-mapped Foldcomp DB with its lookup and index tables.
#[derive(Debug)]
pub struct FoldcompDbReader {
    pub path: String,
    pub input_type: StructureFileFormat,
    pub db_mmap: Mmap,
    pub db: ManuallyDrop<Vec<u8>>,
    pub lookup: FoldcompLookup,
    pub index: FoldcompIndex,
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
        // Both tables are mapped and key-sorted already: no parsing on open.
        let lookup = FoldcompLookup::load(path).expect("Error reading foldcomp db lookup file.");
        let index = FoldcompIndex::load(path).expect("Error reading foldcomp db index file.");
        let path_string_to_return = path.to_string();

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
            lookup: FoldcompLookup::empty(),
            index: FoldcompIndex::empty(),
        }
    }

    /// Decompress one entry by name.
    pub fn read_single_structure(&self, name: &str) -> Result<Structure, String> {
        let mut structure = Structure::new(); // revise
        let mut record = (ChainId::from_byte(b' '), 0);
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

    /// Decompress one entry by DB key.
    pub fn read_single_structure_by_id(&self, id: usize) -> Result<Structure, String> {
        let mut structure = Structure::new(); // revise
        let mut record = (ChainId::from_byte(b' '), 0);
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

    pub fn get_db_key_vector(&self) -> Vec<usize> {
        self.index.keys().collect()
    }
}

// Zero-copy casts between Foldcomp's `atom_t` and `Atom` (identical layout).
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

/// Build a `Structure` from decompressed Foldcomp atoms.
pub unsafe fn atom_t_slice_to_structure(slice: &[atom_t]) -> Structure {
    let mut structure = Structure::new(); 
    let mut record = (ChainId::from_byte(b' '), 0);
    for atom in slice {
        let atom = Atom::from_c(atom);
        structure.update(atom.clone(), &mut record);
    }
    structure
}

/// Entry names in index order, skipping keys the lookup does not name.
pub fn get_path_vector_out_of_lookup_and_index(
    lookup: &FoldcompLookup, index: &FoldcompIndex
) -> Vec<String> {
    index.records().par_iter().map(|record| {
        match lookup.name_of_key(record.key as usize) {
            Some(name) => name.to_string(),
            None => String::new(),
        }
    }).filter(|name| !name.is_empty()).collect()
}

pub fn get_id_vector_out_of_lookup(lookup: &FoldcompLookup) -> Vec<usize> {
    (0..lookup.len()).map(|i| lookup.key(i)).collect()
}

pub fn get_id_vector_subset_out_of_lookup(
    lookup: &FoldcompLookup, subset_names: &Vec<String>
) -> Vec<usize> {
    lookup.keys_of_names(subset_names)
}

pub fn get_name_vector_subset_out_of_lookup(
    lookup: &FoldcompLookup, subset_ids: &Vec<usize>
) -> Vec<String> {
    // Order should be maintained
    subset_ids.iter().filter_map(|id| lookup.name_of_key(*id).map(|n| n.to_string())).collect()
}

/// Map a Foldcomp DB file. The `Vec` aliases the mapping and must never be dropped.
pub fn read_foldcomp_db(db_path: &str) -> Result<(Mmap, ManuallyDrop<Vec<u8>>), &'static str> {
    let db_file = match File::open(&db_path) {
        Ok(file) => file,
        Err(_) => return Err("DB file not found."),
    };
    
    let mmap = unsafe { Mmap::map(&db_file).unwrap() };
    let db = unsafe { ManuallyDrop::new(Vec::from_raw_parts(mmap.as_ptr() as *mut u8, mmap.len(), mmap.len())) };
    Ok((mmap, db))
}

pub fn get_foldcomp_db_entry<'a>(
    db: &'a ManuallyDrop<Vec<u8>>, record: &FoldcompIndexRecord
) -> &'a [u8] {
    let start = record.offset as usize;
    &db[start..start + record.length as usize]
}

pub fn get_foldcomp_db_entry_by_id<'a>(
    db: &'a ManuallyDrop<Vec<u8>>, index: &FoldcompIndex, id: usize
) -> Option<&'a [u8]> {
    let record = index.find(id)?;
    Some(get_foldcomp_db_entry(db, record))
}

pub fn get_foldcomp_db_entry_by_name<'a>(
    db: &'a ManuallyDrop<Vec<u8>>, lookup: &FoldcompLookup, index: &FoldcompIndex, name: &str
) -> Option<&'a [u8]> {
    let key = lookup.key_of_name(name)?;
    get_foldcomp_db_entry_by_id(db, index, key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rayon::iter::IntoParallelRefIterator;
    use std::io::Read;
    use rayon::iter::ParallelIterator;
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
            let (_db_mmap, db) = read_foldcomp_db(db_path).unwrap();
            let lookup = FoldcompLookup::load(db_path).unwrap();
            let index = FoldcompIndex::load(db_path).unwrap();
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
            let _ = &index.records().par_iter().for_each(|entry_index| {
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
#[cfg(test)]
mod cache_tests {
    use super::*;

    /// Private copy of the example DB so test caches stay out of `data/` and unshared.
    struct TempDb {
        prefix: String,
    }

    impl TempDb {
        fn new(tag: &str) -> Self {
            let unique = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
                .unwrap().as_nanos();
            let dir = format!(
                "{}/folddisco_fcdb_{}_{}", std::env::temp_dir().to_string_lossy(), tag, unique
            );
            std::fs::create_dir_all(&dir).unwrap();
            let prefix = format!("{}/example_db", dir);
            for suffix in ["", ".index", ".lookup", ".dbtype", ".source"] {
                let from = format!("data/foldcomp/example_db{}", suffix);
                if std::path::Path::new(&from).exists() {
                    std::fs::copy(&from, format!("{}{}", prefix, suffix)).unwrap();
                }
            }
            TempDb { prefix }
        }
    }

    impl Drop for TempDb {
        fn drop(&mut self) {
            if let Some(dir) = std::path::Path::new(&self.prefix).parent() {
                let _ = std::fs::remove_dir_all(dir);
            }
        }
    }

    fn text_lookup(db: &str) -> Vec<(u64, String)> {
        let mut out: Vec<(u64, String)> = std::fs::read_to_string(format!("{}.lookup", db))
            .unwrap().lines().map(|line| {
                let mut split = line.split('\t');
                let key = split.next().unwrap().parse::<u64>().unwrap();
                (key, split.next().unwrap().to_string())
            }).collect();
        out.sort_by_key(|entry| entry.0);
        out
    }

    fn text_index(db: &str) -> Vec<(u64, u64, u64)> {
        let mut out: Vec<(u64, u64, u64)> = std::fs::read_to_string(format!("{}.index", db))
            .unwrap().lines().map(|line| {
                let mut split = line.split('\t');
                (split.next().unwrap().parse().unwrap(),
                 split.next().unwrap().parse().unwrap(),
                 split.next().unwrap().parse().unwrap())
            }).collect();
        out.sort_by_key(|entry| entry.0);
        out
    }

    /// Built and re-mapped tables both match the text files, in key order.
    #[test]
    fn mapped_tables_match_the_text_files() {
        let db = TempDb::new("roundtrip");
        let expected_lookup = text_lookup(&db.prefix);
        let expected_index = text_index(&db.prefix);

        for pass in ["built from text", "mapped from cache"] {
            let lookup = FoldcompLookup::load(&db.prefix).unwrap();
            let index = FoldcompIndex::load(&db.prefix).unwrap();
            assert_eq!(lookup.len(), expected_lookup.len(), "{pass}");
            assert_eq!(index.len(), expected_index.len(), "{pass}");
            for (i, (key, name)) in expected_lookup.iter().enumerate() {
                assert_eq!(lookup.key(i), *key as usize, "{pass}");
                assert_eq!(lookup.name(i), name, "{pass}");
                assert_eq!(lookup.name_of_key(*key as usize), Some(name.as_str()), "{pass}");
                assert_eq!(lookup.key_of_name(name), Some(*key as usize), "{pass}");
            }
            for (i, (key, offset, length)) in expected_index.iter().enumerate() {
                assert_eq!(index.records()[i], FoldcompIndexRecord {
                    key: *key, offset: *offset, length: *length,
                }, "{pass}");
                assert_eq!(index.find(*key as usize).unwrap().offset, *offset, "{pass}");
            }
            assert_eq!(index.keys().collect::<Vec<_>>(),
                       expected_index.iter().map(|e| e.0 as usize).collect::<Vec<_>>(), "{pass}");
            // The first pass has written the caches the second one maps.
            assert!(std::path::Path::new(&format!("{}.lookup.fdcache", db.prefix)).is_file());
            assert!(std::path::Path::new(&format!("{}.index.fdcache", db.prefix)).is_file());
        }
    }

    /// The cache must not use Foldcomp's own `<index>.cache` name.
    #[test]
    fn cache_paths_do_not_collide_with_foldcomps_own() {
        let db = TempDb::new("collide");
        FoldcompIndex::load(&db.prefix).unwrap();
        assert!(std::path::Path::new(&format!("{}.index.fdcache", db.prefix)).is_file());
        assert!(!std::path::Path::new(&format!("{}.index.cache", db.prefix)).exists());
    }

    #[test]
    fn a_name_not_in_the_db_is_not_found() {
        let db = TempDb::new("missing_name");
        let lookup = FoldcompLookup::load(&db.prefix).unwrap();
        assert_eq!(lookup.key_of_name("AF-NOT-A-REAL-ENTRY-F1-model_v4"), None);
        assert_eq!(lookup.name_of_key(usize::MAX), None);
        assert_eq!(FoldcompIndex::load(&db.prefix).unwrap().find(usize::MAX), None);
    }

    #[test]
    fn keys_of_names_resolves_many_in_one_pass_and_keeps_order() {
        let db = TempDb::new("many_names");
        let lookup = FoldcompLookup::load(&db.prefix).unwrap();
        let expected = text_lookup(&db.prefix);
        assert!(expected.len() >= 2, "the example DB needs at least two entries");

        // Last, first, then a name that is not there at all.
        let asked = vec![
            expected[expected.len() - 1].1.clone(),
            expected[0].1.clone(),
            "AF-NOT-A-REAL-ENTRY-F1-model_v4".to_string(),
        ];
        assert_eq!(
            lookup.keys_of_names(&asked),
            vec![expected[expected.len() - 1].0 as usize, expected[0].0 as usize]
        );
    }

    /// An all-zero (torn) record must be rejected at load.
    #[test]
    fn a_zeroed_record_is_rejected_and_the_text_is_reparsed() {
        use std::io::{Seek, SeekFrom, Write};
        let db = TempDb::new("zeroed");
        let expected = text_index(&db.prefix);
        assert!(expected.len() >= 2);
        FoldcompIndex::load(&db.prefix).unwrap();

        let cache_path = format!("{}.index.fdcache", db.prefix);
        let record_1 = crate::utils::pod_cache::HEADER_SIZE
            + std::mem::size_of::<FoldcompIndexRecord>();
        let mut cache = std::fs::OpenOptions::new().write(true).open(&cache_path).unwrap();
        cache.seek(SeekFrom::Start(record_1 as u64)).unwrap();
        cache.write_all(&[0u8; std::mem::size_of::<FoldcompIndexRecord>()]).unwrap();
        drop(cache);

        // Rejected, so the load falls back to the text and rewrites the cache.
        let index = FoldcompIndex::load(&db.prefix).unwrap();
        assert_eq!(index.len(), expected.len());
        assert_eq!(index.records()[1].key, expected[1].0);
    }

    /// A rewritten lookup must invalidate its cache even if the mtime is unchanged.
    #[test]
    fn a_changed_source_is_rejected() {
        let db = TempDb::new("stale");
        let lookup_path = format!("{}.lookup", db.prefix);
        let cache_path = format!("{}.lookup.fdcache", db.prefix);
        let original = FoldcompLookup::load(&db.prefix).unwrap().len();
        let cache_modified = std::fs::metadata(&cache_path).unwrap().modified().unwrap();

        let text = std::fs::read_to_string(&lookup_path).unwrap();
        let trimmed: String = text.lines().skip(1)
            .map(|l| format!("{}\n", l)).collect();
        std::fs::write(&lookup_path, &trimmed).unwrap();
        let file = std::fs::OpenOptions::new().write(true).open(&lookup_path).unwrap();
        file.set_times(std::fs::FileTimes::new().set_modified(cache_modified)).unwrap();
        drop(file);

        let reloaded = FoldcompLookup::load(&db.prefix).unwrap();
        assert_eq!(reloaded.len(), original - 1, "the cache should have been rebuilt");
    }
}
