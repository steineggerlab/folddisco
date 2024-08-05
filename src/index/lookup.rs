// Save & Load the vector of file names
// Working with controller.path_vec: Vec<String>
// Lookup file format
// id\tpath\tinteger\tfloat
// id\tpath\tn_res\tplddt

use std::io::{Write, BufRead};
use std::fs::File;
use std::io::BufWriter;

use crate::utils::log::{log_msg, FAIL};

pub fn save_lookup_to_file(
    path: &str, path_vec: &Vec<String>, numeric_id_vec: &Vec<usize>, 
    optional_int_vec: Option<&Vec<usize>>, optional_float_vec: Option<&Vec<f32>>
) {
    assert_eq!(path_vec.len(), numeric_id_vec.len());
    if optional_int_vec.is_some() {
        assert_eq!(path_vec.len(), optional_int_vec.unwrap().len());
    }
    if optional_float_vec.is_some() {
        assert_eq!(path_vec.len(), optional_float_vec.unwrap().len());
    }
    
    // Save the vector of file names to a file
    let mut file = BufWriter::new(File::create(path).expect(&log_msg(FAIL, "Unable to create the lookup file")));
    for i in 0..path_vec.len() {
        let line = match (optional_int_vec, optional_float_vec) {
            (Some(int_vec), Some(float_vec)) => {
                format!("{}\t{}\t{}\t{}\n", numeric_id_vec[i], path_vec[i], int_vec[i], float_vec[i])
            },
            (Some(int_vec), None) => {
                format!("{}\t{}\t{}\t{}\n", numeric_id_vec[i], path_vec[i], int_vec[i], 0.0)
            },
            (None, Some(float_vec)) => {
                format!("{}\t{}\t{}\t{}\n", numeric_id_vec[i], path_vec[i], 0, float_vec[i])
            },
            (None, None) => {
                format!("{}\t{}\t{}\t{}\n", numeric_id_vec[i], path_vec[i], 0, 0.0)
            }
        };
        file.write_all(line.as_bytes()).expect(&log_msg(FAIL, "Unable to write the lookup file"));
    }
}

pub fn load_lookup_from_file(path: &str) -> (Vec<String>, Vec<usize>, Vec<usize>, Vec<f32>) {
    let mut path_vec: Vec<String> = Vec::new();
    let mut numeric_id_vec: Vec<usize> = Vec::new();
    let mut integer_vec: Vec<usize> = Vec::new();
    let mut float_vec: Vec<f32> = Vec::new();
    let file = std::fs::File::open(path).expect(&log_msg(FAIL, "Unable to open the lookup file"));
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line = line.expect(&log_msg(FAIL, "Unable to read the lookup file"));
        let line_vec: Vec<&str> = line.split("\t").collect();
        numeric_id_vec.push(line_vec[0].parse::<usize>().unwrap());
        path_vec.push(line_vec[1].to_string());
        integer_vec.push(line_vec[2].parse::<usize>().unwrap());
        float_vec.push(line_vec[3].parse::<f32>().unwrap());
    }
    (path_vec, numeric_id_vec, integer_vec, float_vec)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_save_and_load_lookup() {
        let path = "data/lookup_test.lookup";
        let path_vec = vec!["path1.pdb".to_string(), "path2.pdb".to_string(), "path3.pdb".to_string()];
        let numeric_id_vec = vec![0, 1, 2];
        let nres_vec = Some(vec![100, 200, 5000]);
        let plddt_vec = Some(vec![50.0, 60.0, 70.0]);

        // Save the data to a file
        save_lookup_to_file(path, &path_vec, &numeric_id_vec, nres_vec.as_ref(), plddt_vec.as_ref());

        // Load the data from the file
        let (loaded_path_vec, loaded_numeric_id_vec,
             loaded_optional_vec, loaded_float_vec ) = load_lookup_from_file(path);

        // Check that the loaded data is the same as the original data
        assert_eq!(path_vec, loaded_path_vec);
        assert_eq!(numeric_id_vec, loaded_numeric_id_vec);
        assert_eq!(nres_vec.unwrap(), loaded_optional_vec);
        assert_eq!(plddt_vec.unwrap(), loaded_float_vec);

        // Clean up the test file
        // std::fs::remove_file(path).unwrap();
    }
}