// Save & Load the vector of file names
// Working with controller.path_vec: Vec<String>
// Lookup file format
// id\tpath\toptional

// use BufWriter;
use std::io::{Write, BufRead};
use std::fs::File;
use std::io::BufWriter;

pub fn save_lookup_to_file(path: &str, path_vec: &Vec<String>, numeric_id_vec: &Vec<usize>, optional_vec: Option<&Vec<usize>>) {
    assert_eq!(path_vec.len(), numeric_id_vec.len());
    if optional_vec.is_some() {
        assert_eq!(path_vec.len(), optional_vec.unwrap().len());
    }
    // Save the vector of file names to a file
    let mut file = BufWriter::new(File::create(path).expect("Unable to create file"));
    for i in 0..path_vec.len() {
        let line = if optional_vec.is_some() {
            format!("{}\t{}\t{}\n", numeric_id_vec[i], path_vec[i], optional_vec.unwrap()[i])
        } else {
            format!("{}\t{}\t{}\n", numeric_id_vec[i], path_vec[i], 0)
        };
        file.write_all(line.as_bytes()).expect("Unable to write data");
    }
}

pub fn load_lookup_from_file(path: &str) -> (Vec<String>, Vec<usize>, Vec<usize>) {
    let mut path_vec: Vec<String> = Vec::new();
    let mut numeric_id_vec: Vec<usize> = Vec::new();
    let mut optional_vec: Vec<usize> = Vec::new();
    let file = std::fs::File::open(path).expect("Unable to open file");
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        let line_vec: Vec<&str> = line.split("\t").collect();
        numeric_id_vec.push(line_vec[0].parse::<usize>().unwrap());
        path_vec.push(line_vec[1].to_string());
        optional_vec.push(line_vec[2].parse::<usize>().unwrap());
    }
    (path_vec, numeric_id_vec, optional_vec)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_save_and_load_lookup() {
        let path = "data/lookup_test.lookup";
        let path_vec = vec!["path1.pdb".to_string(), "path2.pdb".to_string(), "path3.pdb".to_string()];
        let numeric_id_vec = vec![0, 1, 2];
        let optional_vec = Some(vec![0, 0, 0]);

        // Save the data to a file
        save_lookup_to_file(path, &path_vec, &numeric_id_vec, optional_vec.as_ref());

        // Load the data from the file
        let (loaded_path_vec, loaded_numeric_id_vec, loaded_optional_vec) = load_lookup_from_file(path);

        // Check that the loaded data is the same as the original data
        assert_eq!(path_vec, loaded_path_vec);
        assert_eq!(numeric_id_vec, loaded_numeric_id_vec);
        assert_eq!(optional_vec.unwrap(), loaded_optional_vec);

        // Clean up the test file
        std::fs::remove_file(path).unwrap();
    }
}