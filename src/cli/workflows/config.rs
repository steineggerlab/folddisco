
use std::io::{BufRead, Write};
use crate::prelude::{HashType, log_msg, FAIL};

pub fn write_config_to_file(path: &str, hash_type: HashType, num_bin_dist: usize, num_bin_angle: usize) {
    // Write as UTF-8 text file
    let mut file = std::fs::File::create(path).expect("Unable to create file");
    // First line: hash type
    file.write_all(format!("{:?}\n", hash_type).as_bytes()).unwrap();
    // second line: num_bin_dist
    file.write_all(format!("{}\n", num_bin_dist).as_bytes()).unwrap();
    // third line: num_bin_angle
    file.write_all(format!("{}\n", num_bin_angle).as_bytes()).unwrap();
}

pub fn read_config_from_file(path: &str) -> (HashType, usize, usize) {
    let file = std::fs::File::open(path).expect(
        &log_msg(FAIL, &format!("Config file not found: {}", path))
    );
    let reader = std::io::BufReader::new(file);
    let mut lines = reader.lines();
    let hash_type = lines.next().unwrap().unwrap();
    let num_bin_dist = lines.next().unwrap().unwrap();
    let num_bin_angle = lines.next().unwrap().unwrap();
    let hash_type = HashType::get_with_str(&hash_type);
    let num_bin_dist = num_bin_dist.parse::<usize>().unwrap();
    let num_bin_angle = num_bin_angle.parse::<usize>().unwrap();
    (hash_type, num_bin_dist, num_bin_angle)
}