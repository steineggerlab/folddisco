// Module: index
// Declare inner modules
pub mod alloc;
pub mod builder;

use std::collections::HashMap;
use std::fmt::Display;
use std::hash::Hash;
use std::io::Write;
use std::sync::Arc;

// TODO: Generalize this to handle multiple hash types

/// ## Index table saves the hash and the list of ids that have the same hash.
/// - Key: hash
/// - Value: list of ids
pub type IndexTable<T, U> = HashMap<T, Vec<U>>;

pub enum IndexTablePrinter {
    Text,
    Binary,
    Debug,
}

impl IndexTablePrinter {
    pub fn print<T: Hash + Display, U: Display + Hash + Ord + Eq>(
        &self,
        index_table: &IndexTable<T, U>,
        path: &str,
    ) {
        match self {
            IndexTablePrinter::Text => write_index_table_text(index_table, path),
            IndexTablePrinter::Binary => write_index_table_binary(index_table, path),
            IndexTablePrinter::Debug => write_index_table_debug(index_table, path),
        }
    }
}

fn write_index_table_text<T: Hash + Display, U: Display + Hash + Ord + Eq>(
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

// WARNING: not working now
fn write_index_table_binary<T: Hash + Display, U: Display>(
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
        let line = format!("{}\t{}\n", key.to_string(), value_comma_separated);
        // Write text file
        file.write_all(line.as_bytes())
            .expect("Unable to write data");
    }
}

fn write_index_table_debug<T: Hash + Display, U: Display + Hash + Ord + Eq>(
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

pub struct OffsetTable<T: Hash>(pub HashMap<T, usize>);
