pub mod alloc;
pub mod builder;
use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Display;
use std::io::Write;

pub type IndexTable<T: Hash, U> = HashMap<T, Vec<U>>;

pub fn write_index_table<T: Hash + Display, U: Display>(index_table: &IndexTable<T, U>, path: &str) {
    let mut file = std::fs::File::create(path).expect("Unable to create file");
    // file.write_all(b"hash\tdist\tangle\tids\n").expect("Unable to write data");
    for (key, value) in index_table {
        let value_comma_separated = value.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",");
        // No quotation marks
        let line = format!("{}\t{}\n", key.to_string(), value_comma_separated);
        // Write text file
        file.write_all(line.as_bytes()).expect("Unable to write data");
    }
}
