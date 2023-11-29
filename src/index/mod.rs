// Module: index
// Declare inner modules
pub mod alloc;
pub mod builder;
pub mod index_table;
pub mod lookup;
pub mod delta;
pub mod io;

use std::collections::HashMap;
use std::fmt::Display;
use std::hash::Hash;
use std::io::Write;

// TODO: Generalize this to handle multiple hash types

/// ## Index table saves the hash and the list of ids that have the same hash.
/// - Key: hash
/// - Value: list of ids
pub type IndexTable<T, U> = HashMap<T, Vec<U>>;

pub enum IndexTableFileType {
    Json,
    Tsv,
    Binary,
}

pub enum IndexTablePrinter {
    Text,
    Binary,
    Debug,
}

impl IndexTablePrinter {
    pub fn print<T: Hash + Eq + Display, U: Display + Hash + Ord + Eq>(
        &self,
        index_table: &IndexTable<T, U>,
        path: &str,
    ) {
        match self {
            IndexTablePrinter::Text => write_index_table_text(index_table, path),
            IndexTablePrinter::Binary => unimplemented!(),
            IndexTablePrinter::Debug => write_index_table_debug(index_table, path),
        }
    }
}

fn write_index_table_text<T: Hash + Eq + Display, U: Display + Hash + Ord + Eq>(
    index_table: &IndexTable<T, U>,
    path: &str,
) {
    let mut file = std::fs::File::create(path).expect("Unable to create file");
    for (key, value) in index_table {
        let value_comma_separated = value
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(",");
        // No quotation marks
        let line = format!("{}\t{}\n", key.to_string(), value_comma_separated);
        // Write text file
        file.write_all(line.as_bytes())
            .expect("Unable to write data");
    }
}

fn write_index_table_debug<T: Hash + Eq + Display, U: Display + Hash + Ord + Eq>(
    index_table: &IndexTable<T, U>,
    path: &str,
) {
    let mut file = std::fs::File::create(path).expect("Unable to create file");
    file.write_all(b"hash\tdist\tomega\tphi1\tphi2\tpsi1\tpsi2\tids\tidcount\tid_unique\n")
        .expect("Unable to write data");
    for (key, value) in index_table {
        // Calculate id count & add
        let value_comma_separated = value
            .iter()
            .fold(HashMap::new(), |mut acc, x| {
                *acc.entry(x).or_insert(0) += 1;
                acc
            })
            .into_iter()
            .map(|(x, y)| format!("{}:{}", x, y))
            .collect::<Vec<String>>()
            .join(",");
        let count = value.len();
        let uniq_count = value.iter().collect::<std::collections::HashSet<_>>().len();
        let line = format!(
            "{}\t{}\t{}\t{}\n",
            key.to_string(),
            value_comma_separated,
            count,
            uniq_count
        ); // Added id count for debugging
           // Write text file
        file.write_all(line.as_bytes())
            .expect("Unable to write data");
    }
}

pub fn query_single<T: Hash + Eq, U: Hash + Eq + Clone>(
    index_table: &IndexTable<T, U>,
    query: &T,
) -> Option<Vec<U>> {
    match index_table.get(query) {
        Some(x) => Some(x.clone()),
        None => None,
    }
}

pub fn query_multiple<T: Hash + Eq, U: Hash + Eq + Clone>(
    index_table: &IndexTable<T, U>,
    queries: &[T],
) -> Option<Vec<U>> {
    // query_multiple returns the intersection of the results of the queries
    // TODO: Generated. NEEDS TESTING
    let mut result: Option<Vec<U>> = None;
    for query in queries {
        let query_result = query_single(index_table, query);
        match query_result {
            Some(x) => {
                match result {
                    Some(y) => {
                        // Filter result if result is not in query_result
                        result = Some(y.into_iter().filter(|z| x.contains(z)).collect());
                    }
                    None => {
                        result = Some(x);
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


pub fn query_multiple_with_neighbors<T: Hash + Eq, U: Hash + Eq + Clone>(
    index_table: &IndexTable<T, U>,
    queries_with_neighbors: Vec<Vec<T>>,
) -> Option<Vec<U>> {
    // query_multiple returns the intersection of the results of the queries
    // TODO: Generated. NEEDS TESTING
    let mut result: Option<Vec<U>> = None;
    for query_list in queries_with_neighbors {
        let mut neighbor_result: Option<Vec<U>> = None;
        for query in query_list {
            let query_result = query_single(index_table, &query);
            match query_result {
                Some(x) => {
                    match neighbor_result {
                        Some(y) => {
                            // Union result if result is in query_result
                            // Element should be unique in neighbor_result
                            neighbor_result = Some(y.into_iter().chain(x).collect::<std::collections::HashSet<_>>().into_iter().collect());
                        }
                        None => {
                            neighbor_result = Some(x);
                        }
                    }
                }
                None => {
                    neighbor_result = neighbor_result;
                }
            }
        }
        // Intersect result if result is in neighbor_result
        match result {
            Some(y) => {
                result = Some(y.into_iter().filter(|z| neighbor_result.as_ref().unwrap().contains(z)).collect());
            }
            None => {
                result = neighbor_result;
            }
        }
    }
    result
}

pub struct OffsetTable<T: Hash>(pub HashMap<T, usize>);
