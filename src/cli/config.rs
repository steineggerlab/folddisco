
use std::{fs, io::{BufRead, Write}};
use crate::prelude::{HashType, log_msg, FAIL};
use toml::map::Map;
use crate::controller::mode::IndexMode;

#[derive(Debug, Clone, PartialEq)]
pub struct IndexConfig {
    pub hash_type: HashType,
    pub num_bin_dist: usize,
    pub num_bin_angle: usize,
    pub mode: IndexMode,
    pub grid_width: f32,
    pub chunk_size: usize,
    pub max_residue: usize,
}

impl IndexConfig {
    pub fn new(
        hash_type: HashType, num_bin_dist: usize, num_bin_angle: usize,
        mode: IndexMode, grid_width: f32, chunk_size: usize, max_residue: usize
    ) -> Self {
        Self {
            hash_type,
            num_bin_dist,
            num_bin_angle,
            mode,
            grid_width,
            chunk_size,
            max_residue,
        }
    }
    pub fn from_toml(toml: &toml::Value) -> Self {
        let hash_type = toml["hash_type"].as_str().unwrap();
        let num_bin_dist = toml["num_bin_dist"].as_integer().unwrap() as usize;
        let num_bin_angle = toml["num_bin_angle"].as_integer().unwrap() as usize;
        let mode = toml["mode"].as_str().unwrap();
        let mode = IndexMode::get_with_str(mode);
        let grid_width = toml["grid_width"].as_float().unwrap() as f32;
        let chunk_size = toml["chunk_size"].as_integer().unwrap() as usize;
        let max_residue = toml["max_residue"].as_integer().unwrap() as usize;
        Self {
            hash_type: HashType::get_with_str(hash_type),
            num_bin_dist,
            num_bin_angle,
            mode,
            grid_width,
            chunk_size,
            max_residue,
        }
    }
    pub fn to_toml(&self) -> toml::Value {
        let mut map = Map::new();
        map.insert("hash_type".to_string(), toml::Value::String(self.hash_type.to_string()));
        map.insert("num_bin_dist".to_string(), toml::Value::Integer(self.num_bin_dist as i64));
        map.insert("num_bin_angle".to_string(), toml::Value::Integer(self.num_bin_angle as i64));
        map.insert("mode".to_string(), toml::Value::String(self.mode.to_string()));
        map.insert("grid_width".to_string(), toml::Value::Float(self.grid_width as f64));
        map.insert("chunk_size".to_string(), toml::Value::Integer(self.chunk_size as i64));
        map.insert("max_residue".to_string(), toml::Value::Integer(self.max_residue as i64));
        toml::Value::Table(map)
    }
}


#[derive(Debug, Clone, PartialEq)]
pub struct QueryConfig {
    pub retrieve: bool,
    pub amino_acid: u8,
    pub dist_threshold: Vec<f32>,
    pub angle_threshold: Vec<f32>,
    pub match_cutoff: Vec<f32>,
    pub score_cutoff: f32,
    pub num_res_cutoff: usize,
    pub plddt_cutoff: f32,
}

impl QueryConfig {
    pub fn new(
        retrieve: bool, amino_acid: u8, dist_threshold: Vec<f32>,
        angle_threshold: Vec<f32>, match_cutoff: Vec<f32>,
        score_cutoff: f32, num_res_cutoff: usize, plddt_cutoff: f32
    ) -> Self {
        Self {
            retrieve,
            amino_acid,
            dist_threshold,
            angle_threshold,
            match_cutoff,
            score_cutoff,
            num_res_cutoff,
            plddt_cutoff,
        }
    }
    pub fn from_toml(toml: &toml::Value) -> Self {
        let retrieve = toml["retrieve"].as_bool().unwrap();
        let amino_acid = toml["amino_acid"].as_integer().unwrap() as u8;
        let dist_threshold = toml["dist_threshold"].as_array().unwrap().iter().map(
            |x| x.as_float().unwrap() as f32
        ).collect();
        let angle_threshold = toml["angle_threshold"].as_array().unwrap().iter().map(
            |x| x.as_float().unwrap() as f32
        ).collect();
        let match_cutoff = toml["match_cutoff"].as_array().unwrap().iter().map(
            |x| x.as_float().unwrap() as f32
        ).collect();
        let score_cutoff = toml["score_cutoff"].as_float().unwrap() as f32;
        let num_res_cutoff = toml["num_res_cutoff"].as_integer().unwrap() as usize;
        let plddt_cutoff = toml["plddt_cutoff"].as_float().unwrap() as f32;
        Self {
            retrieve,
            amino_acid,
            dist_threshold,
            angle_threshold,
            match_cutoff,
            score_cutoff,
            num_res_cutoff,
            plddt_cutoff,
        }
    }
    pub fn to_toml(&self) -> toml::Value {
        let mut map = Map::new();
        map.insert("retrieve".to_string(), toml::Value::Boolean(self.retrieve));
        map.insert("amino_acid".to_string(), toml::Value::Integer(self.amino_acid as i64));
        map.insert("dist_threshold".to_string(), toml::Value::Array(
            self.dist_threshold.iter().map(|x| toml::Value::Float(*x as f64)).collect())
        );
        map.insert("angle_threshold".to_string(), toml::Value::Array(
            self.angle_threshold.iter().map(|x| toml::Value::Float(*x as f64)).collect())
        );
        map.insert("match_cutoff".to_string(), toml::Value::Array(
            self.match_cutoff.iter().map(|x| toml::Value::Float(*x as f64)).collect())
        );
        map.insert("score_cutoff".to_string(), toml::Value::Float(self.score_cutoff as f64));
        map.insert("num_res_cutoff".to_string(), toml::Value::Integer(self.num_res_cutoff as i64));
        map.insert("plddt_cutoff".to_string(), toml::Value::Float(self.plddt_cutoff as f64));
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

pub fn write_query_config_to_file(path: &str, query_config: QueryConfig) {
    let mut file = std::fs::File::create(path).expect(
        &log_msg(FAIL, &format!("Unable to create config file: {}", path))
    );
    let toml = query_config.to_toml();
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

pub fn read_query_config_from_file(path: &str) -> QueryConfig {
    let file = std::fs::File::open(path).expect(
        &log_msg(FAIL, &format!("Config file not found: {}", path))
    );
    let reader = std::io::BufReader::new(file);
    let toml = toml::from_str(&reader.lines().map(|x| format!("{}\n",x.unwrap())).collect::<String>()).unwrap();
    QueryConfig::from_toml(&toml)
}

pub fn write_configs_to_file(path: &str, index_config: Option<IndexConfig>, query_config: Option<QueryConfig>) {
    let mut toml = Map::new();
    if let Some(index_config) = index_config {
        toml.insert("index".to_string(), index_config.to_toml());
    }
    if let Some(query_config) = query_config {
        toml.insert("query".to_string(), query_config.to_toml());
    }
    fs::write(path, toml::to_string(&toml).unwrap()).unwrap();
}

pub fn read_configs_from_file(path: &str) -> (Option<IndexConfig>, Option<QueryConfig>) {
    let file = std::fs::File::open(path).expect(
        &log_msg(FAIL, &format!("Config file not found: {}", path))
    );
    let reader = std::io::BufReader::new(file);
    let toml: Map<String, toml::Value> = toml::from_str(&reader.lines().map(|x| format!("{}\n",x.unwrap())).collect::<String>()).unwrap();
    let index_config = if let Some(index) = toml.get("index") {
        Some(IndexConfig::from_toml(index))
    } else {
        None
    };
    let query_config = if let Some(query) = toml.get("query") {
        Some(QueryConfig::from_toml(query))
    } else {
        None
    };
    (index_config, query_config)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_write_configs_to_file() {
        let path = "data/configs.toml";
        let index_config = Some(
            IndexConfig::new(HashType::PDBTrRosetta, 10, 10, IndexMode::Id, 30.0, 65535, 4000)
        );
        let query_config = Some(
            QueryConfig::new(true, 0, vec![0.0], vec![0.0], vec![0.0], 0.0, 50000, 0.0)
        );
        write_configs_to_file(path, index_config.clone(), query_config.clone());
        let (index_config_read, query_config_read) = read_configs_from_file(path);
        assert_eq!(index_config, index_config_read);
        assert_eq!(query_config, query_config_read);
    }
    
    #[test]
    fn test_write_index_config_to_file() {
        let path = "data/index_config.toml";
        let index_config = IndexConfig::new(
            HashType::PDBTrRosetta, 10, 10, IndexMode::Grid, 30.0, 65535, 4000
        );
        write_index_config_to_file(path, index_config.clone());
        let index_config_read = read_index_config_from_file(path);
        assert_eq!(index_config, index_config_read);
    }
    
    #[test]
    fn test_write_query_config_to_file() {
        let path = "data/query_config.toml";
        let query_config = QueryConfig::new(
            true, 0, vec![0.0], vec![0.0], vec![0.0], 0.0, 50000, 0.0
        );
        write_query_config_to_file(path, query_config.clone());
        let query_config_read = read_query_config_from_file(path);
        assert_eq!(query_config, query_config_read);
    }
}