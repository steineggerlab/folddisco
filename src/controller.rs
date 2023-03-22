use crate::geometry::hash::{HashCollection, HashValue};
use crate::index::IndexTable;

pub struct Controller {
    pub path: String,
}

impl Controller {
    pub fn new(path: String) -> Controller {
        Controller { path: path }
    }

    pub fn collect_hash() {
        todo!()
    }
}

pub struct GeometryHashCollector {
    pub hash_collection: HashCollection,
    // TODO: FILL IN OTHER FIELDS
}

impl GeometryHashCollector {
    pub fn new() -> GeometryHashCollector {
        GeometryHashCollector {
            hash_collection: HashCollection::new(),
        }
    }

    pub fn collect_hash(&mut self, hash_value: HashValue) {
        self.hash_collection.push(hash_value);
    }

    pub fn remove_redundancy(&mut self) {
        self.hash_collection.sort();
        self.hash_collection.dedup();
    }
}