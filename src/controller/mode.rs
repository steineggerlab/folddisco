// Structure ID formats for index lookups, and query output modes.
use std::collections::HashSet;
use std::fs;
use std::path::Path;

/// How a structure path is turned into the ID stored in the lookup (`--id`).
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

/// Structure ID for `path`: e.g. `pdb1abc.ent` -> `1abc` (Pdb), `AF-P17538-F1-model_v4` (Afdb),
/// `P17538` (UniProt). Unmatched AFDB/UniProt names fall back to the file stem.
#[inline]
pub fn parse_path_by_id_type(path: &str, id_type: &IdType) -> String {
    let afdb_regex = regex::Regex::new(r"AF-.+-model_v\d").unwrap();
    match id_type {
        IdType::Pdb => {
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap();
            let file_name = file_name.to_str().unwrap();
            // Strip the `pdb` prefix of wwPDB file names
            if file_name.starts_with("pdb") {
                file_name[3..].to_string()
            } else {
                file_name.to_string()
            }
        }
        IdType::Afdb => {
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap().to_str().unwrap();
            let afdb_id = afdb_regex.find(file_name);
            if afdb_id.is_none() {
                return file_name.to_string();
            } 
            file_name[afdb_id.unwrap().start()..afdb_id.unwrap().end()].to_string()
        }
        IdType::UniProt => {
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap().to_str().unwrap();
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

/// `parse_path_by_id_type` writing into a reused buffer.
#[inline]
pub fn parse_path_by_id_type_with_string(path: &str, id_type: &IdType, string: &mut String) {
    string.clear();
    let afdb_regex = regex::Regex::new(r"AF-.+-model_v\d").unwrap();
    match id_type {
        IdType::Pdb => {
            let path = Path::new(path);
            let file_name = path.file_stem().unwrap();
            let file_name = file_name.to_str().unwrap();
            // Strip the `pdb` prefix of wwPDB file names
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


/// `parse_path_by_id_type` over a list.
pub fn parse_path_vec_by_id_type(path_vec: &Vec<String>, id_type: &IdType) -> Vec<String> {
    let mut parsed_path_vec = Vec::with_capacity(path_vec.len());
    for path in path_vec {
        parsed_path_vec.push(parse_path_by_id_type(&path, id_type));
    }
    parsed_path_vec
}

/// `parse_path_by_id_type` over a set.
pub fn parse_path_set_by_id_type(path_set: &HashSet<String>, id_type: &IdType) -> HashSet<String> {
    let mut parsed_path_set = HashSet::with_capacity(path_set.len());
    for path in path_set {
        parsed_path_set.insert(parse_path_by_id_type(&path, id_type));
    }
    parsed_path_set
}

/// Query output mode, derived from the output flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueryMode {
    PerMatch,      // Default: one line per residue match
    PerStructure,  // One line per structure
    SkipMatch,     // Per structure, without residue matching
    Web,           // Per match, capped and with superposition
    ContradictoryPrintError, // --per-structure and --per-match together
}

impl QueryMode {
    pub fn from_flags(
        skip_match: bool, 
        is_web: bool,
        per_structure: bool, 
        per_match: bool,
    ) -> Self {
        match (skip_match, is_web, per_structure, per_match) {
            (_, _, true, true) => Self::ContradictoryPrintError,
            (_, true, _, _) => Self::Web,
            (true, _, _, _) => Self::SkipMatch,
            (false, false, true, false) => Self::PerStructure,
            (false, false, false, _) => Self::PerMatch,
        }
    }

    pub fn to_string(&self) -> String {
        match self {
            QueryMode::PerMatch => "per match".to_string(),
            QueryMode::PerStructure => "per structure".to_string(),
            QueryMode::SkipMatch => "per structure skipping match".to_string(),
            QueryMode::Web => "for web".to_string(),
            QueryMode::ContradictoryPrintError => "contradictory print error".to_string(),
        }
    }
}

impl std::fmt::Display for QueryMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            QueryMode::PerMatch => write!(f, "per match"),
            QueryMode::PerStructure => write!(f, "per structure"),
            QueryMode::SkipMatch => write!(f, "per structure skipping match"),
            QueryMode::Web => write!(f, "for web"),
            QueryMode::ContradictoryPrintError => write!(f, "contradictory print error"),
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
    fn test_parse_path_by_id_type() {
        let pdb_path = "data/serine_peptidases/1azw.pdb";
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
        assert_eq!(rel_path, "data/serine_peptidases/1azw.pdb");
    }
}