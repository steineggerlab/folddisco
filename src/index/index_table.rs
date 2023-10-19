// File: index_table.rs
// Created: 2023-10-19 14:53:14
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

// IndexTable is a HashMap that maps a key to a vector of values.
// The key is a hash value for a given residue-residue pair. So, the HashMap
// should be able to handle custom hash functions.
// Value is a vector of ids that have the same hash value.
// Alternatively, value can be a vector of (id, res_id1, res_id2) tuples.
// IndexTable should be saved as / loaded from a text file or a binary file.
// Text file would be tab-separated values (TSV) file or JSON file.
// Binary file would be a serialized HashMap.

// use std::collections::HashMap;
// Using FxHashMap instead of HashMap for faster performance
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{Read, Write};
use serde::{Serialize, Deserialize};
use bincode;
use serde_json;

// Define the types for Key and Id
// 
pub type Key = u64; // Hash value for a given residue-residue pair
pub type Id = usize; // Id of a protein structure
pub type ResId = u16; // Id of a residue

// Value can be a single Id or a tuple of (id, res_id1, res_id2)
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Copy, Clone, Hash)]
pub enum Value {
    Single(Id),
    Triple(Id, ResId, ResId),
}

impl Value {
    pub fn get_id(&self) -> usize {
        match self {
            Value::Single(id) => *id,
            Value::Triple(id, _, _) => *id,
        }
    }
}


// IndexTable is a HashMap that maps a key to a vector of values
#[derive(Serialize, Deserialize, Debug)]
pub struct IndexTable(FxHashMap<Key, Vec<Value>>);

impl IndexTable {
    // Create a new IndexTable
    pub fn new() -> Self {
        IndexTable(FxHashMap::default())
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }
    
    pub fn with_capacity(capacity: usize) -> Self {
        IndexTable(FxHashMap::with_capacity_and_hasher(capacity, Default::default()))
    }

    // Save the IndexTable to a JSON file
    pub fn save_to_json(&self, filename: &str) -> std::io::Result<()> {
        let file = File::create(filename)?;
        serde_json::to_writer(file, self)?;
        Ok(())
    }

    pub fn get(&self, key: &Key) -> Option<&Vec<Value>> {
        self.0.get(key)
    }

    // Load the IndexTable from a JSON file
    pub fn load_from_json(filename: &str) -> std::io::Result<Self> {
        let mut file = File::open(filename)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;
        let table = serde_json::from_str(&contents)?;
        Ok(table)
    }

    // Save the IndexTable to a binary file
    pub fn save_to_bin(&self, filename: &str) -> std::io::Result<()> {
        let file = File::create(filename)?;
        let bin_result = bincode::serialize_into(file, self);
        if let Err(e) = bin_result {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to serialize: {}", e),
            ));
        }
        Ok(())
    }

    // Load the IndexTable from a binary file
    pub fn load_from_bin(filename: &str) -> std::io::Result<Self> {
        let file = File::open(filename)?;
        let table_load_result = bincode::deserialize_from(file);
        if let Err(e) = table_load_result {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to deserialize: {}", e),
            ));
        }
        let table = table_load_result.unwrap();
        Ok(table)
    }

    // Add a value to the IndexTable
    pub fn add(&mut self, key: Key, value: Value) {
        self.0.entry(key).or_insert_with(Vec::new).push(value);
    }
    
    pub fn query_single(&self, query: &Key) -> Option<Vec<Value>> {
        match self.get(query) {
            Some(x) => Some(x.clone()),
            None => None,
        }
    }
    
    pub fn query_multiple(&self, queries: &Vec<Key>) -> Option<Vec<Value>> {
        let mut result: Option<Vec<Value>> = None;
        for query in queries {
            let query_result = self.get(query);
            match query_result {
                Some(x) => {
                    match &mut result {
                        Some(y) => {
                            y.retain(|item| x.contains(item));
                        }
                        None => {
                            result = Some(x.clone());
                        }
                    }
                }
                None => {
                    return None;
                }
            }
        }
        result
    }

    pub fn query_multiple_with_neighbors(&self, queries_with_neighbors: &Vec<Vec<Key>>) -> Option<Vec<Value>> {
        let mut result: Option<Vec<Value>> = None;
        for query_list in queries_with_neighbors {
            let mut neighbor_result: Option<Vec<Value>> = None;
            for query in query_list {
                let query_result = self.get(&query);
                match query_result {
                    Some(x) => {
                        match &mut neighbor_result {
                            Some(y) => {
                                y.extend(x.iter().cloned());
                            }
                            None => {
                                neighbor_result = Some(x.clone());
                            }
                        }
                    }
                    None => {
                        // Do nothing
                    }
                }
            }
            match &mut result {
                Some(y) => {
                    if let Some(ref neighbor_result) = neighbor_result {
                        y.retain(|item| neighbor_result.contains(item));
                    }
                }
                None => {
                    result = neighbor_result;
                }
            }
        }
        result
    }

    pub fn concat_structure_with_res_pair(&mut self, id: Id, hashes: Vec<Key>, res_pairs: Vec<(ResId, ResId)>) {
        // This function consumes hashes and res_pairs
        assert_eq!(hashes.len(), res_pairs.len(), "Length of hashes and res_pairs must be equal.");
        for (hash, res_pair) in hashes.into_iter().zip(res_pairs.into_iter()) {
            self.add(hash, Value::Triple(id, res_pair.0, res_pair.1));
        }
    }
    
    pub fn concat_structure(&mut self, id: Id, hashes: Vec<Key>) {
        // This function consumes hashes
        for hash in hashes.into_iter() {
            self.add(hash, Value::Single(id));
        }
    }

}


pub fn build_from_ids_and_hashes(ids: Vec<Id>, hash_collections: Vec<Vec<Key>>) -> IndexTable {
    assert_eq!(ids.len(), hash_collections.len(), "Length of ids and hash_collections must be equal.");

    let total_hashes = hash_collections.iter().map(|hashes| hashes.len()).sum();
    let mut table = IndexTable(FxHashMap::with_capacity_and_hasher(total_hashes, Default::default()));
    for (id, hashes) in ids.into_iter().zip(hash_collections.into_iter()) {
        for hash in hashes {
            table.add(hash, Value::Single(id));
        }
    }
    table
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_table() {
        let mut table = IndexTable::new();
        table.add(123, Value::Single(1));
        table.add(123, Value::Triple(1, 2, 3));
        table.save_to_json("data/table.json").unwrap();
        let loaded_table = IndexTable::load_from_json("data/table.json").unwrap();
        // Assert that the key & value are the same
        assert_eq!(table.0, loaded_table.0);
    }
    #[test]
    fn test_index_table_binary() {
        let mut table = IndexTable::new();
        table.add(123, Value::Single(1));
        table.add(456, Value::Single(2));
        table.add(123, Value::Single(3));
        // table.add(123, Value::Triple(1, 2, 3));
        table.save_to_bin("data/table.bin").unwrap();
        let loaded_table = IndexTable::load_from_bin("data/table.bin").unwrap();
        assert_eq!(table.0, loaded_table.0);
    }
}