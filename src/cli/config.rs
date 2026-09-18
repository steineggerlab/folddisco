//! Index configuration stored next to an index as `<index>.type` (TOML).

use std::io::{BufRead, Write};
use crate::prelude::{HashType, log_msg, FAIL};
use toml::map::Map;
use crate::structure::io::StructureFileFormat;
use crate::controller::expand::IndexExpansion;
use crate::controller::substitution::SubstitutionScheme;

/// Hashing and input settings an index was built with; queries must reuse them.
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
    /// Build-time expansion; absent in plain and older indices.
    pub expansion: Option<IndexExpansion>,
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
            expansion: None,
        }
    }
    /// Parse a config table; panics on missing required keys.
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
        let expansion = toml.get("expand_radius").and_then(|radius| IndexExpansion::new(
            radius.as_integer().unwrap_or(0) as usize,
            toml.get("expand_distance").and_then(|x| x.as_float()).unwrap_or(0.0) as f32,
            toml.get("expand_angle").and_then(|x| x.as_float()).unwrap_or(0.0) as f32,
            toml.get("aa_subst").and_then(|x| x.as_str()).and_then(SubstitutionScheme::from_str),
        ));
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
            expansion,
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
        if let Some(expansion) = &self.expansion {
            let widest = |t: &[f32]| t.iter().fold(0.0f32, |w, x| w.max(x.abs())) as f64;
            map.insert("expand_radius".to_string(), toml::Value::Integer(expansion.tolerance.radius as i64));
            map.insert("expand_distance".to_string(), toml::Value::Float(widest(&expansion.tolerance.dist_thresholds)));
            map.insert("expand_angle".to_string(), toml::Value::Float(widest(&expansion.tolerance.angle_thresholds)));
            if let Some(scheme) = expansion.scheme {
                map.insert("aa_subst".to_string(), toml::Value::String(scheme.to_string()));
            }
        }
        toml::Value::Table(map)
    }
}



/// Write `index_config` as TOML to `path`.
pub fn write_index_config_to_file(path: &str, index_config: IndexConfig) {
    let mut file = std::fs::File::create(path).expect(
        &log_msg(FAIL, &format!("Unable to create config file: {}", path))
    );
    let toml = index_config.to_toml();
    file.write_all(toml::to_string(&toml).unwrap().as_bytes()).unwrap();
}

/// Read a TOML index config; panics with a message if the file is missing.
pub fn read_index_config_from_file(path: &str) -> IndexConfig {
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
        assert_eq!(index_config_read.expansion, None);
    }

    #[test]
    fn expansion_round_trips_through_the_config_file() {
        let path = std::env::temp_dir().join(format!("folddisco_expansion_{}.type", std::process::id()));
        let path = path.to_str().unwrap();
        for expansion in [
            IndexExpansion::new(1, 0.5, 5.0, Some(SubstitutionScheme::Blosum62)),
            IndexExpansion::new(2, 1.0, 10.0, None),
            IndexExpansion::new(0, 0.5, 5.0, Some(SubstitutionScheme::Size)),
        ] {
            let mut index_config = IndexConfig::new(
                HashType::PDBTrRosetta, 16, 4, 20.0, 10, 50000, StructureFileFormat::PDB, None, None,
            );
            assert!(expansion.is_some());
            index_config.expansion = expansion;
            write_index_config_to_file(path, index_config.clone());
            assert_eq!(read_index_config_from_file(path), index_config);
        }
        std::fs::remove_file(path).ok();
    }
}