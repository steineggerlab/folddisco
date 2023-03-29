pub mod alloc;
pub mod builder;
use std::collections::HashMap;
use std::hash::Hash;
use std::fmt::Display;
use std::io::Write;


pub type IndexTable<T: Hash, U> = HashMap<T, Vec<U>>;

pub fn write_index_table<T: Hash + Display, U: Display + Hash + Ord + Eq>(index_table: &IndexTable<T, U>, path: &str) {
    let mut file = std::fs::File::create(path).expect("Unable to create file");
    // WARNING: Temporary header
    file.write_all(b"hash\tdist\tangle\tids\tidcount\tid_unique\n").expect("Unable to write data");
    for (key, value) in index_table {
        let value_comma_separated = value.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",");
        // No quotation marks
        // let line = format!("{}\t{}\n", key.to_string(), value_comma_separated);

        let count = value.len();
        let uniq_count = value.iter().collect::<std::collections::HashSet<_>>().len();
        let line = format!(
            "{}\t{}\t{}\t{}\n",
            key.to_string(), value_comma_separated, count, uniq_count
        ); // Added id count for debugging
        // Write text file
        file.write_all(line.as_bytes()).expect("Unable to write data");
    }
}

// WARNING: not finished
pub fn write_index_table_binary<T: Hash + Display, U: Display>(index_table: &IndexTable<T, U>, path: &str) {
    let mut file = std::fs::File::create(path).expect("Unable to create file");
    for (key, value) in index_table {
        let value_comma_separated = value.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",");
        let line = format!("{}\t{}\n", key.to_string(), value_comma_separated);
        // Write text file
        file.write_all(line.as_bytes()).expect("Unable to write data");
    }
}


pub struct OffsetTable<T: Hash> {
    pub offset_table: HashMap<T, usize>,
    pub big_allocation: Box<[usize]>,
}