
use std::collections::HashSet;
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdType {
    Pdb,
    Afdb,
    UniProt,
    BasenameWithoutExt,
    BasenameWithExt,
    AbsPath,
    RelPath,
    Other,
}

impl IdType {
    pub fn get_with_str(id_type: &str) -> Self {
        match id_type {
            "Pdb" | "PDB" | "pdb" => Self::Pdb,
            "Afdb" | "AFDB" | "afdb" => Self::Afdb,
            "Uniprot" | "UniProt" | "uniprot" => Self::UniProt,
            "BasenameWithoutExt" | "basename_without_ext" | "basename_no_ext" | "filename" => Self::BasenameWithoutExt,
            "BasenameWithExt" | "basename_with_ext" | "basename" | "file" => Self::BasenameWithExt,
            "AbsPath" | "Abspath" | "abspath" | "absolute_path" | "path" => Self::AbsPath,
            "RelPath" | "Relpath" | "relpath" | "relative_path" | "default" => Self::RelPath,
            _ => Self::Other,
        }
    }
    pub fn to_string(&self) -> String {
        match self {
            Self::Pdb => "pdb".to_string(),
            Self::Afdb => "afdb".to_string(),
            Self::UniProt => "uniprot".to_string(),
            Self::BasenameWithoutExt => "basename_without_ext".to_string(),
            Self::BasenameWithExt => "basename_with_ext".to_string(),
            Self::AbsPath => "absolute_path".to_string(),
            Self::RelPath => "relative_path".to_string(),
            Self::Other => "other".to_string(),
        }
    }
    pub fn get_with_u8(id_type: u8) -> Self {
        match id_type {
            0 => Self::Pdb,
            1 => Self::Afdb,
            2 => Self::UniProt,
            3 => Self::BasenameWithoutExt,
            4 => Self::BasenameWithExt,
            5 => Self::AbsPath,
            6 => Self::RelPath,
            _ => Self::Other,
        }
    }
    pub fn to_u8(&self) -> u8 {
        match self {
            Self::Pdb => 0,
            Self::Afdb => 1,
            Self::UniProt => 2,
            Self::BasenameWithoutExt => 3,
            Self::BasenameWithExt => 4,
            Self::AbsPath => 5,
            Self::RelPath => 6,
            Self::Other => 7,
        }
    }
}

#[inline]
pub fn parse_path_by_id_type(path: &str, id_type: &IdType) -> String {
    // TODO: 2024-04-04 15:07:54 Fill in this function to ease benchmarking
    let afdb_regex = regex::Regex::new(r"AF-.+-model_v\d").unwrap();
    match id_type {
        IdType::Pdb => {
            // Get the basename of the path
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap();
            // Remove extension
            let file_name = file_name.to_str().unwrap();
            // Remove extension, If startswith "pdb" remove "pdb" from the start
            if file_name.starts_with("pdb") {
                file_name[3..].to_string()
            } else {
                file_name.to_string()
            }
        }
        IdType::Afdb => {
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap().to_str().unwrap();
            // Find the matching pattern
            let afdb_id = afdb_regex.find(file_name);
            if afdb_id.is_none() {
                return file_name.to_string();
            } 
            file_name[afdb_id.unwrap().start()..afdb_id.unwrap().end()].to_string()
        }
        IdType::UniProt => {
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap().to_str().unwrap();
            // Find the matching pattern
            let afdb_id = afdb_regex.find(file_name);
            if afdb_id.is_none() {
                return file_name.to_string();
            } 
            let afdb_id = file_name[afdb_id.unwrap().start()..afdb_id.unwrap().end()].to_string();
            let afdb_id = afdb_id.split("-").collect::<Vec<_>>();
            afdb_id[1].to_string()
        }
        IdType::BasenameWithoutExt => {
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap().to_str().unwrap();
            file_name.to_string()
        }
        IdType::BasenameWithExt => {
            let path = Path::new(path);
            let file_name = path.file_name().unwrap().to_str().unwrap();
            file_name.to_string()
        }
        IdType::AbsPath => {
            let path = fs::canonicalize(path).unwrap();
            path.to_str().unwrap().to_string()
        }
        IdType::RelPath => path.to_string(),
        IdType::Other => path.to_string(),
    }
}

#[inline]
pub fn parse_path_by_id_type_with_string(path: &str, id_type: &IdType, string: &mut String) {
    // TODO: 2024-04-04 15:07:54 Fill in this function to ease benchmarking
    string.clear();
    let afdb_regex = regex::Regex::new(r"AF-.+-model_v\d").unwrap();
    match id_type {
        IdType::Pdb => {
            // Get the basename of the path
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap();
            // Remove extension
            let file_name = file_name.to_str().unwrap();
            // Remove extension, If startswith "pdb" remove "pdb" from the start
            if file_name.starts_with("pdb") {
                // &file_name[3..]
                string.push_str(&file_name[3..]);
            } else {
                // file_name
                string.push_str(file_name);
            }
        }
        IdType::Afdb => {
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap().to_str().unwrap();
            // Find the matching pattern
            let afdb_id = afdb_regex.find(file_name);
            if afdb_id.is_none() {
                // return file_name;
                string.push_str(file_name);
            } else {
                // &file_name[afdb_id.unwrap().start()..afdb_id.unwrap().end()]
                string.push_str(&file_name[afdb_id.unwrap().start()..afdb_id.unwrap().end()]);
            }
        }
        IdType::UniProt => {
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap().to_str().unwrap();
            // Find the matching pattern
            let afdb_id = afdb_regex.find(file_name);
            if afdb_id.is_none() {
                // return file_name;
                string.push_str(file_name);
            } 
            let afdb_id = file_name[afdb_id.unwrap().start()..afdb_id.unwrap().end()].to_string();
            let afdb_id = afdb_id.split("-").collect::<Vec<_>>();
            // afdb_id[1]
            string.push_str(afdb_id[1]);
        }
        IdType::BasenameWithoutExt => {
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap().to_str().unwrap();
            // file_name
            string.push_str(file_name);
        }
        IdType::BasenameWithExt => {
            let path = Path::new(path);
            let file_name = path.file_name().unwrap().to_str().unwrap();
            // file_name
            string.push_str(file_name);
        }
        IdType::AbsPath => {
            let path = fs::canonicalize(path).unwrap();
            // path.to_str().unwrap()
            string.push_str(path.to_str().unwrap());
        }
        IdType::RelPath => {
            // path
            string.push_str(path);
        }
        IdType::Other => {
            // path
            string.push_str(path);
        }
    }
}


pub fn parse_path_vec_by_id_type(path_vec: &Vec<String>, id_type: &IdType) -> Vec<String> {
    let mut parsed_path_vec = Vec::with_capacity(path_vec.len());
    for path in path_vec {
        parsed_path_vec.push(parse_path_by_id_type(&path, id_type));
    }
    parsed_path_vec
}

pub fn parse_path_set_by_id_type(path_set: &HashSet<String>, id_type: &IdType) -> HashSet<String> {
    let mut parsed_path_set = HashSet::with_capacity(path_set.len());
    for path in path_set {
        parsed_path_set.insert(parse_path_by_id_type(&path, id_type));
    }
    parsed_path_set
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexMode {
    Id,
    Big,
}

impl IndexMode {
    pub fn get_with_str(mode: &str) -> Self {
        match mode {
            "Id" | "id" | "ID" => Self::Id,
            "Big" | "big" | "BIG" => Self::Big,
            _ => Self::Id,
        }
    }
    pub fn to_string(&self) -> String {
        match self {
            Self::Id => "id".to_string(),
            Self::Big => "big".to_string(),
        }
    }
    pub fn get_with_u8(mode: u8) -> Self {
        match mode {
            0 => Self::Id,
            1 => Self::Big,
            _ => Self::Id,
        }
    }
    pub fn to_u8(&self) -> u8 {
        match self {
            Self::Id => 0,
            Self::Big => 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_id_type() {
        let id_type = IdType::Pdb;
        assert_eq!(id_type.to_string(), "pdb");
        assert_eq!(id_type.to_u8(), 0);
        assert_eq!(IdType::get_with_u8(0), IdType::Pdb);
        assert_eq!(IdType::get_with_str("pdb"), IdType::Pdb);
    }

    #[test]
    fn test_index_mode() {
        let index_mode = IndexMode::Id;
        assert_eq!(index_mode.to_string(), "id");
        assert_eq!(index_mode.to_u8(), 0);
        assert_eq!(IndexMode::get_with_u8(0), IndexMode::Id);
        assert_eq!(IndexMode::get_with_str("id"), IndexMode::Id);
    }

    #[test]
    fn test_parse_path_by_id_type() {
        let pdb_path = "data/serine_peptidases_filtered/1azw.pdb";
        let afdb_path = "data/AF-P17538-F1-model_v4.pdb";

        let pdb_id = parse_path_by_id_type(pdb_path, &IdType::Pdb);
        let afdb_id = parse_path_by_id_type(afdb_path, &IdType::Afdb);
        let uniprot_id = parse_path_by_id_type(afdb_path, &IdType::UniProt);
        let basename_ext_id = parse_path_by_id_type(afdb_path, &IdType::BasenameWithExt);
        let basename_no_ext_id = parse_path_by_id_type(afdb_path, &IdType::BasenameWithoutExt);
        let abs_path = parse_path_by_id_type(afdb_path,&IdType::AbsPath);
        let rel_path = parse_path_by_id_type(pdb_path, &IdType::RelPath);
        
        assert_eq!(pdb_id, "1azw");
        assert_eq!(afdb_id, "AF-P17538-F1-model_v4");
        assert_eq!(uniprot_id, "P17538");
        assert_eq!(basename_ext_id, "AF-P17538-F1-model_v4.pdb");
        assert_eq!(basename_no_ext_id, "AF-P17538-F1-model_v4");
        println!("abs_path: {}", abs_path);
        assert_eq!(rel_path, "data/serine_peptidases_filtered/1azw.pdb");
    }
}