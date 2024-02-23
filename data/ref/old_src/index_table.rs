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

use std::io::{BufReader, BufWriter};
use std::io::Error;
use byteorder::{WriteBytesExt, ReadBytesExt, LittleEndian as LE};

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
pub struct IndexTable(pub FxHashMap<Key, Vec<Value>>);

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
        let mut file = File::create(filename)?;
        let encoded = bincode::serialize(self).unwrap();
        println!("Encoded length: {}", encoded.len());
        let bin_result = file.write_all(&encoded);
        // let bin_result = bincode::serialize_into(file, self);
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
        let raw: Vec<u8> = file.bytes().map(|x| x.unwrap()).collect();
        let table_load_result = bincode::deserialize::<IndexTable>(&raw);
        // let table_load_result = bincode::deserialize_from(file);

        if let Err(e) = table_load_result {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                format!("Failed to deserialize: {}", e),
            ));
        }
        let table = table_load_result.unwrap();
        Ok(table)
    }

    pub fn save_to_bin_custom(&self, path: &str) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        let mut count = 0;
        for (key, values) in &self.0 {
            count += 1;
            if count < 5 {
                println!("Key: {}, Values: {:?}", key, values.len());
            }
            writer.write_u64::<LE>(*key)?;
            writer.write_u64::<LE>(values.len() as u64)?;
            for value in values {
                match value {
                    // If single, write 0u8
                    Value::Single(value) => {
                        writer.write_u8(0u8)?;
                        writer.write_u64::<LE>(*value as u64)?;
                    },
                    Value::Triple(id, res_id1, res_id2) => {
                        writer.write_u8(1u8)?;
                        writer.write_u32::<LE>(*id as u32)?;
                        writer.write_u16::<LE>(*res_id1 as u16)?;
                        writer.write_u16::<LE>(*res_id2 as u16)?;
                    }
                }
            }
        }

        Ok(())
    }
    // TODO: IMPORTANT: Remove serde dependency
    pub fn load_from_bin_custom(path: &str) -> std::io::Result<IndexTable> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        let mut index_table = IndexTable::new();
        let mut count = 0;
        while let Ok(key) = reader.read_u64::<LE>() {
            let len = reader.read_u64::<LE>()? as u64;
            
            count += 1;
            if count < 5 {
                println!("Key: {}, Len: {}", key, len);
            }
            for _ in 0..len {
                let value_type = reader.read_u8()?;
                match value_type {
                    0u8 => {
                        let value = reader.read_u64::<LE>()? as usize;
                        index_table.add(key, Value::Single(value));
                    },
                    1u8 => {
                        let id = reader.read_u32::<LE>()? as usize;
                        let res_id1 = reader.read_u16::<LE>()? as u16;
                        let res_id2 = reader.read_u16::<LE>()? as u16;
                        index_table.add(key, Value::Triple(id, res_id1, res_id2));
                    },
                    _ => {
                        return Err(Error::new(std::io::ErrorKind::Other, "Invalid value type"));
                    }
                }
            }
        }
        Ok(index_table)
    }

    
    // Add a value to the IndexTable
    pub fn add(&mut self, key: Key, value: Value) {
        self.0.entry(key).or_insert_with(Vec::new).push(value);
    }
    pub fn remove(&mut self, key: &Key) {
        self.0.remove(key);
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
            println!("Query: {:?}, Query result: {:?}", query, query_result);
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
    pub fn query_multiple_with_count(&self, queries: &Vec<Key>, min_count: usize) -> Option<Vec<Value>> {
        let mut counts: FxHashMap<Value, usize> = FxHashMap::with_capacity_and_hasher(queries.len(), Default::default());
        for query in queries {
            let query_result = self.get(query);
            match query_result {
                Some(x) => {
                    for value in x {
                        *counts.entry(*value).or_insert(0) += 1;
                    }
                }
                None => {
                    return None;
                }
            }
        }
        let result: Vec<Value> = counts.into_iter()
            .filter(|&(_, count)| count >= min_count)
            .map(|(value, _)| value)
            .collect();
        Some(result)
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

    pub fn query_multiple_with_connectivity(&self, queries: &Vec<Key>, min_count: usize) -> Option<Vec<Value>> {
        let mut count_pair_map: FxHashMap<Value, (usize, Vec<(u16, u16)>)> = FxHashMap::with_capacity_and_hasher(queries.len(), Default::default());
        for query in queries {
            let query_result = self.get(query);
            match query_result {
                Some(x) => {
                    for value in x {
                        if let Value::Triple(id, res_id1, res_id2) = value {
                            let entry = count_pair_map.entry(*value).or_insert((0, Vec::new()));
                            entry.0 += 1;
                            entry.1.push((*res_id1, *res_id2));
                        }
                    }
                }
                None => { } // do nothing
            }
        }
        // Filter1: count over min_count
        // Filter2: res_index pair is connected with other res_index pair
        let result: Vec<Value> = count_pair_map.into_iter()
            .filter(|&(_, (count, _))| count >= min_count)
            .filter(|(_, (_, res_pairs))| {
                let mut connected = false;
                for i in 0..res_pairs.len() {
                    for j in i+1..res_pairs.len() {
                        let res_pair1 = res_pairs[i];
                        let res_pair2 = res_pairs[j];
                        if res_pair1.0 == res_pair2.0 || res_pair1.0 == res_pair2.1 || res_pair1.1 == res_pair2.0 || res_pair1.1 == res_pair2.1 {
                            connected = true;
                            break;
                        }
                    }
                    if connected {
                        break;
                    }
                }
                connected
            })
            .map(|(value, _)| value)
            .collect();
        Some(result)
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
    pub fn fill(&mut self, id_vec: &Vec<Id>, hash_collection_vec: &Vec<Vec<Key>>) {
        assert_eq!(id_vec.len(), hash_collection_vec.len(), "Length of id_vec and hash_collection_vec must be equal.");
        for i in 0..hash_collection_vec.len() {
            let hash_collection = &hash_collection_vec[i];
            let id = id_vec[i];
            for hash in hash_collection {
                self.add(*hash, Value::Single(id));
            }
        }
    }
    pub fn fill_triple(&mut self, id_vec: &Vec<Id>, hash_vec: &Vec<Vec<Value>>) {
        // assert value is triple
        assert_eq!(id_vec.len(), hash_vec.len(), "Length of id_vec and hash_vec must be equal.");
        for i in 0..hash_vec.len() {
            let hash_collection = &hash_vec[i];
            let id = id_vec[i];
            for hash in hash_collection {
                self.add(id as Key, *hash);
            }
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