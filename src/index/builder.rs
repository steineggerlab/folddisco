use crate::index::IndexTable;
use std::cmp::{Eq, Ord};
use std::hash::Hash;
pub struct IndexBuilder {}

// pub trait Hashable: Hash + Ord + Eq + Clone {}

impl IndexBuilder {
    pub fn new() -> IndexBuilder {
        IndexBuilder {}
    }
    pub fn concat<T: Hash + Ord + Eq + Clone, U: Clone + Copy>(
        &self,
        id_vec: &Vec<U>,
        hash_collection_vec: &Vec<Vec<T>>,
    ) -> IndexTable<T, U> {
        let mut index_table = IndexTable::new();
        for i in 0..hash_collection_vec.len() {
            let hash_collection = &hash_collection_vec[i];
            let id = id_vec[i].clone();
            for hash in hash_collection {
                let id_vec = index_table.entry(hash.clone()).or_insert(Vec::new());
                id_vec.push(id);
            }
        }
        index_table
    }
    
    
}

#[cfg(test)]
mod index_builder_tests {
    use super::*;
    use crate::index::IndexTablePrinter;

    #[test]
    fn test_concat() {
        let index_builder = IndexBuilder {};
        let id_vec: Vec<u16> = vec![1, 2, 3, 4, 5];
        let hash_collection_vec: Vec<Vec<u16>> = vec![
            vec![1000, 1001, 1002, 1003, 1004],
            vec![1002, 1004, 1006, 1008, 1010],
            vec![1001, 1003, 1005, 1006],
            vec![1000, 1002, 1004, 1006, 1008, 1010],
            vec![1001, 1003, 1005, 1007, 1009, 1011],
        ];
        let index_table = index_builder.concat(&id_vec, &hash_collection_vec);
        assert_eq!(index_table[&1000], vec![1, 4]);
        assert_eq!(index_table[&1001], vec![1, 3, 5]);
        assert_eq!(index_table[&1002], vec![1, 2, 4]);
        assert_eq!(index_table[&1003], vec![1, 3, 5]);
        assert_eq!(index_table[&1004], vec![1, 2, 4]);
        assert_eq!(index_table[&1005], vec![3, 5]);
        assert_eq!(index_table[&1006], vec![2, 3, 4]);
        assert_eq!(index_table[&1007], vec![5]);
        assert_eq!(index_table[&1008], vec![2, 4]);
        assert_eq!(index_table[&1009], vec![5]);
        assert_eq!(index_table[&1010], vec![2, 4]);
        assert_eq!(index_table[&1011], vec![5]);

        let table_printer: IndexTablePrinter = IndexTablePrinter::Text;
        table_printer.print(&index_table, "data/index_table_test.txt");
    }
}
