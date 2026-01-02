
use std::io::{BufRead, Write};
use crate::prelude::{HashType, log_msg, FAIL};
use toml::map::Map;
use crate::structure::io::StructureFileFormat;

#[derive(Debug, Clone, PartialEq)]
pub struct IndexConfig {
    pub hash_type: HashType,
    pub num_bin_dist: usize,
    pub num_bin_angle: usize,
    pub grid_width: f32,
    pub chunk_size: usize,
    pub max_residue: usize,
    pub input_format: StructureFileFormat,
    pub foldcomp_db: Option<String>,
    pub multiple_bin: Option<Vec<(usize, usize)>>,
}

impl IndexConfig {
    pub fn new(
        hash_type: HashType, num_bin_dist: usize, num_bin_angle: usize,
        grid_width: f32, chunk_size: usize, max_residue: usize,
        input_format: StructureFileFormat, foldcomp_db: Option<String>,
        multiple_bin: Option<Vec<(usize, usize)>>,
    ) -> Self {
        Self {
            hash_type,
            num_bin_dist,
            num_bin_angle,
            grid_width,
            chunk_size,
            max_residue,
            input_format,
            foldcomp_db,
            multiple_bin,
        }
    }
    pub fn from_toml(toml: &toml::Value) -> Self {
        let hash_type = toml["hash_type"].as_str().unwrap();
        let num_bin_dist = toml["num_bin_dist"].as_integer().unwrap() as usize;
        let num_bin_angle = toml["num_bin_angle"].as_integer().unwrap() as usize;
        let grid_width = toml["grid_width"].as_float().unwrap() as f32;
        let chunk_size = toml["chunk_size"].as_integer().unwrap() as usize;
        let max_residue = toml["max_residue"].as_integer().unwrap() as usize;
        let input_format = StructureFileFormat::get_with_string(toml["input_format"].as_str().unwrap());
        let foldcomp_db = toml.get("foldcomp_db").map(|x| x.as_str().unwrap().to_string());
        let multiple_bin = toml.get("multiple_bin").map(|x| {
            x.as_array().unwrap().iter().map(|y| {
                let bin = y.as_array().unwrap();
                (bin[0].as_integer().unwrap() as usize, bin[1].as_integer().unwrap() as usize)
            }).collect()
        });
        Self {
            hash_type: HashType::get_with_str(hash_type),
            num_bin_dist,
            num_bin_angle,
            grid_width,
            chunk_size,
            max_residue,
            input_format,
            foldcomp_db,
            multiple_bin,
        }
    }
    pub fn to_toml(&self) -> toml::Value {
        let mut map = Map::new();
        map.insert("hash_type".to_string(), toml::Value::String(self.hash_type.to_string()));
        map.insert("num_bin_dist".to_string(), toml::Value::Integer(self.num_bin_dist as i64));
        map.insert("num_bin_angle".to_string(), toml::Value::Integer(self.num_bin_angle as i64));
        map.insert("grid_width".to_string(), toml::Value::Float(self.grid_width as f64));
        map.insert("chunk_size".to_string(), toml::Value::Integer(self.chunk_size as i64));
        map.insert("max_residue".to_string(), toml::Value::Integer(self.max_residue as i64));
        map.insert("input_format".to_string(), toml::Value::String(self.input_format.to_string()));
        if let Some(foldcomp_db) = &self.foldcomp_db {
            map.insert("foldcomp_db".to_string(), toml::Value::String(foldcomp_db.clone()));
        }
        if let Some(multiple_bin) = &self.multiple_bin {
            map.insert("multiple_bin".to_string(), toml::Value::Array(
                multiple_bin.iter().map(|x| {
                    toml::Value::Array(vec![toml::Value::Integer(x.0 as i64), toml::Value::Integer(x.1 as i64)])
                }).collect()
            ));
        }
        toml::Value::Table(map)
    }
}



pub fn write_index_config_to_file(path: &str, index_config: IndexConfig) {
    let mut file = std::fs::File::create(path).expect(
        &log_msg(FAIL, &format!("Unable to create config file: {}", path))
    );
    let toml = index_config.to_toml();
    file.write_all(toml::to_string(&toml).unwrap().as_bytes()).unwrap();
}

pub fn read_index_config_from_file(path: &str) -> IndexConfig {
    // Read file as read only, read from start to end including new line characters
    let file = std::fs::File::open(path).expect(
        &log_msg(FAIL, &format!("Config file not found: {}", path))
    );
    let reader = std::io::BufReader::new(file);
    let toml = toml::from_str(&reader.lines().map(|x| format!("{}\n",x.unwrap())).collect::<String>()).unwrap();
    IndexConfig::from_toml(&toml)
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_write_index_config_to_file() {
        let path = "data/index_config.toml";
        let index_config = IndexConfig::new(
            HashType::PDBTrRosetta, 10, 10,
            30.0, 65535, 4000,
            StructureFileFormat::FCZDB, Some("data/foldcomp_db".to_string()),
            Some(vec![(16, 4), (8, 3)])
        );
        write_index_config_to_file(path, index_config.clone());
        let index_config_read = read_index_config_from_file(path);
        assert_eq!(index_config, index_config_read);
    }
}