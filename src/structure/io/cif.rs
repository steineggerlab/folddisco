// MMCIF reader
// Modified from pdbtbx
use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::structure::atom::Atom;

use super::super::core::*;
use super::*;

use pdbtbx_cif::lex_item::{DataBlock, Item, DataItem, Loop, Value};
use pdbtbx_cif::error::{PDBError, ErrorLevel, Context};


/// A PDB reader
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    /// The underlying reader
    pub reader: R,
    ///
    pub input_type: StructureFileFormat,
}

// ??? trait Read -> impl Read for __ ???
impl Reader<File> {
    pub fn new(file: File) -> Self {
        Reader {
            reader: file,
            input_type: StructureFileFormat::CIF,
        }
    }

    /// Read from a file path
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<Self, &'static str> {
        File::open(&path)
            .map(Reader::new)
            .map_err(|_e| "Error opening file")
    }

    pub fn read_structure(&self) -> Result<Structure, &str> {
        let mut reader = BufReader::new(&self.reader);
        let mut structure = Structure::new(); // revise
        let mut contents = String::new();
        if reader.read_to_string(&mut contents).is_ok() {
            match pdbtbx_cif::lex_cif(contents.as_str()) {
                Ok(data_block) => parse_mmcif_block_into_structure(&data_block, &mut structure),
                Err(e) => {
                    eprintln!("Error parsing CIF file: {:?}", e);
                    return Err("Error parsing CIF file");
                }
            }
        } else {
            return Err("Error reading file");
        }
        
        Ok(structure)
    }

    pub fn read_structure_from_gz(&self) -> Result<Structure, &str> {
        // Load whole file and close the file
        let mut decoder = GzDecoder::new(&self.reader);
        let mut binary = Vec::new();
        decoder.read_to_end(&mut binary).unwrap();
        // Close the decoder after flushing
        decoder.flush().unwrap();
        // Drop the decoder
        drop(decoder);

        // Create a new Structure
        let mut structure = Structure::new();
        
        // Read binary as a string. Conver
        let mut reader = BufReader::new(&binary[..]);
        let mut contents = String::new();
        if reader.read_to_string(&mut contents).is_ok() {
            match pdbtbx_cif::lex_cif(contents.as_str()) {
                Ok(data_block) => parse_mmcif_block_into_structure(&data_block, &mut structure),
                Err(e) => {
                    eprintln!("Error parsing CIF file: {:?}", e);
                    return Err("Error parsing CIF file");
                }
            }
        } else {
            return Err("Error reading file");
        }
        

        drop(binary);
        Ok(structure)
    }
    
}


fn parse_mmcif_block_into_structure(input: &DataBlock, structure: &mut Structure) {
    let mut errors: Vec<PDBError> = Vec::new();
    let mut record = (b' ', 0);
    for item in &input.items {
        let result = match item {
            Item::DataItem(di) => match di {
                DataItem::Loop(multiple) => {
                    if multiple.header.contains(&"atom_site.group_PDB".to_string()) {
                        parse_atoms(multiple, structure, &mut record)
                    } else {
                        None
                    }
                }
                _ => None,
            },
            _ => None,
        };
        if let Some(e) = result {
            errors.extend(e);
        }
    }
    if !errors.is_empty() {
        for error in errors {
            eprintln!("{}", error);
        }
    }
}


/// Flatten a Result of a Result with the same error type (#70142 is still unstable)
fn flatten_result<T, E>(value: Result<Result<T, E>, E>) -> Result<T, E> {
    match value {
        Ok(Ok(t)) => Ok(t),
        Ok(Err(e)) => Err(e),
        Err(e) => Err(e),
    }
}

/// Parse a loop containing atomic data
fn parse_atoms(
    input: &Loop, structure: &mut Structure, record: &mut (u8, u64)
) -> Option<Vec<PDBError>> {
    #[derive(Eq, PartialEq)]
    /// The mode of a column
    enum Mode {
        /// A required column (has to be defined)
        Required,
        /// An optional column, if undefined it will have a default value
        Optional,
    }
    use Mode::{Optional, Required};

    /// Easily define all columns
    macro_rules! define_columns {
        ($($i:expr, $name:ident, $label:expr, $req:expr);+;) => {
            $(const $name: (usize, &str, Mode) = ($i, $label, $req);)+
            const COLUMNS: &[(Mode, &str)] = &[
                $(($req, $name.1)),+
            ];
        };
    }

    define_columns!(
        0,  ATOM_ALT_ID, "atom_site.label_alt_id", Optional;
        1,  ATOM_ANISOU_1_1, "_atom_site.aniso_U[1][1]", Optional;
        2,  ATOM_ANISOU_1_2, "_atom_site.aniso_U[1][2]", Optional;
        3,  ATOM_ANISOU_1_3, "_atom_site.aniso_U[1][3]", Optional;
        4,  ATOM_ANISOU_2_1, "_atom_site.aniso_U[2][1]", Optional;
        5,  ATOM_ANISOU_2_2, "_atom_site.aniso_U[2][2]", Optional;
        6,  ATOM_ANISOU_2_3, "_atom_site.aniso_U[2][3]", Optional;
        7,  ATOM_ANISOU_3_1, "_atom_site.aniso_U[3][1]", Optional;
        8,  ATOM_ANISOU_3_2, "_atom_site.aniso_U[3][2]", Optional;
        9,  ATOM_ANISOU_3_3, "_atom_site.aniso_U[3][3]", Optional;
        10, ATOM_ASYM_ID, "atom_site.label_asym_id", Required;
        11, ATOM_AUTH_ASYM_ID, "atom_site.auth_asym_id", Optional;
        12, ATOM_B, "atom_site.B_iso_or_equiv", Optional;
        13, ATOM_CHARGE, "atom_site.pdbx_formal_charge", Optional;
        14, ATOM_COMP_ID, "atom_site.label_comp_id", Required;
        15, ATOM_GROUP, "atom_site.group_PDB", Optional;
        16, ATOM_ID, "atom_site.id", Required;
        17, ATOM_INSERTION, "atom_site.pdbx_PDB_ins_code", Optional;
        18, ATOM_MODEL, "atom_site.pdbx_PDB_model_num", Optional;
        19, ATOM_NAME, "atom_site.label_atom_id", Required;
        20, ATOM_OCCUPANCY, "atom_site.occupancy", Optional;
        21, ATOM_SEQ_ID, "atom_site.label_seq_id", Required;
        22, ATOM_AUTH_SEQ_ID, "atom_site.auth_seq_id", Optional;
        23, ATOM_TYPE, "atom_site.type_symbol", Required;
        24, ATOM_X, "atom_site.Cartn_x", Required;
        25, ATOM_Y, "atom_site.Cartn_y", Required;
        26, ATOM_Z, "atom_site.Cartn_z", Required;
    );

    let positions_: Vec<Result<Option<usize>, PDBError>> = COLUMNS
        .iter()
        .map(|tag| (input.header.iter().position(|t| t == tag.1), tag))
        .map(|(pos, tag)| match pos {
            Some(p) => Ok(Some(p)),
            None if tag.0 == Required => Err(PDBError::new(
                ErrorLevel::InvalidatingError,
                "Missing column in coordinate atoms data loop",
                "The above column is missing",
                Context::show(tag.1),
            )),
            None => Ok(None),
        })
        .collect();

    let mut errors = positions_
        .iter()
        .filter_map(|i| i.clone().err())
        .collect::<Vec<_>>();

    if !errors.is_empty() {
        return Some(errors);
    }

    // Currently, ignoring atom deduplicte check from original code

    // The previous lines make sure that there is no error in the vector.
    let positions: Vec<Option<usize>> = positions_.iter().map(|i| *i.as_ref().unwrap()).collect();
    let mut first_model_number: usize = 0;
    for (index, row) in input.data.iter().enumerate() {
        let values: Vec<Option<&Value>> = positions.iter().map(|i| i.map(|x| &row[x])).collect();
        let context = Context::show(format!("Main atomic data loop row: {index}"));

        /// Parse a column given the function to use and the column index
        macro_rules! parse_column {
            ($type:tt, $index:tt) => {
                if let Some(value) = values[$index.0] {
                    match $type(value, &context, Some($index.1)) {
                        Ok(t) => t,
                        Err(e) => {
                            errors.push(e);
                            None
                        }
                    }
                } else {
                    None
                }
            };
        }

        // Early return cases
        // let element = parse_column!(get_text, ATOM_TYPE).expect("Atom element should be provided");
        let model_number = parse_column!(get_usize, ATOM_MODEL).unwrap_or(1);
        // Use only first model
        if index == 0 {
            first_model_number = model_number;
        } else if model_number != first_model_number {
            break;
        }

        // Parse remaining fields in the order they appear in the line

        let name = parse_column!(get_four_char_array, ATOM_NAME).expect("Atom name should be provided");
        let id: u64 = parse_column!(get_isize, ATOM_ID).expect("Atom ID should be provided") as u64;
        let residue_name: [u8; 3] = parse_column!(get_three_char_array, ATOM_COMP_ID).expect("Residue name should be provided");
        let residue_number: u64 = parse_column!(get_isize, ATOM_AUTH_SEQ_ID).unwrap_or_else(|| {
            parse_column!(get_isize, ATOM_SEQ_ID)
                .expect("Residue number should be provided")
        }) as u64;
        let chain_name = parse_column!(get_one_char, ATOM_AUTH_ASYM_ID).unwrap_or_else(|| {
            parse_column!(get_one_char, ATOM_ASYM_ID).expect("Chain name should be provided")
        });
        let pos_x = parse_column!(get_f32, ATOM_X).expect("Atom X position should be provided");
        let pos_y = parse_column!(get_f32, ATOM_Y).expect("Atom Y position should be provided");
        let pos_z = parse_column!(get_f32, ATOM_Z).expect("Atom Z position should be provided");
        let b_factor = parse_column!(get_f32, ATOM_B).unwrap_or(1.0);
        // Current version does not support Occupancy, Charge and Anisotropic temperature factors 

        // NOT handling HETATM with this version
        // let atom_type = parse_column!(get_text, ATOM_GROUP).unwrap_or_else(|| "ATOM".to_string());
        // let hetero = if atom_type == "ATOM" {
        //     false
        // } else if atom_type == "HETATM" {
        //     true
        // } else {
        //     true
        // };

        let atom = Atom::new(
            pos_x, pos_y, pos_z, name, id,
            chain_name, residue_name, residue_number, b_factor
        );
        
        structure.update(atom, record);
    }

    if !errors.is_empty() {
        Some(errors)
    } else {
        None
    }
}

/// Get the Textual content of the value, if available
fn get_four_char_array(
    value: &Value,
    _context: &Context,
    _column: Option<&str>,
) -> Result<Option<[u8; 4]>, PDBError> {
    match value {
        Value::Text(t) => {
            match t.as_bytes().len() {
                1 => Ok(Some([b' ', t.as_bytes()[0], b' ', b' '])),
                2 => Ok(Some([b' ', t.as_bytes()[0], t.as_bytes()[1], b' '])),
                3 => Ok(Some([b' ', t.as_bytes()[0], t.as_bytes()[1], t.as_bytes()[2]])),
                4 => Ok(Some([t.as_bytes()[0], t.as_bytes()[1], t.as_bytes()[2], t.as_bytes()[3]])),
                _ => Err(PDBError::new(
                    ErrorLevel::InvalidatingError,
                    "Invalid atom name",
                    "Invalid atom name",
                    _context.clone(),
                )),
            }
        },
        _ => Ok(None),
    }
}

fn get_three_char_array(
    value: &Value,
    _context: &Context,
    _column: Option<&str>,
) -> Result<Option<[u8; 3]>, PDBError> {
    match value {
        Value::Text(t) => {
            match t.as_bytes().len() {
                1 => Ok(Some([t.as_bytes()[0], b' ', b' '])),
                2 => Ok(Some([t.as_bytes()[0], t.as_bytes()[1], b' '])),
                3 => Ok(Some([t.as_bytes()[0], t.as_bytes()[1], t.as_bytes()[2]])),
                _ => Err(PDBError::new(
                    ErrorLevel::InvalidatingError,
                    "Invalid residue name",
                    "Invalid residue name",
                    _context.clone(),
                )),
            }
        },
        _ => Ok(None),
    }
}

fn get_one_char(
    value: &Value,
    _context: &Context,
    _column: Option<&str>,
) -> Result<Option<u8>, PDBError> {
    match value {
        Value::Text(t) => {
            match t.as_bytes().len() {
                1 => Ok(Some(t.as_bytes()[0])),
                _ => Err(PDBError::new(
                    ErrorLevel::InvalidatingError,
                    "Invalid chain name",
                    "Currently only one character chain names are supported",
                    _context.clone(),
                )),
            }
        },
        _ => Ok(None),
    }
}

// fn get_text(
//     value: &Value,
//     _context: &Context,
//     _column: Option<&str>,
// ) -> Result<Option<String>, PDBError> {
//     match value {
//         Value::Text(t) => Ok(Some(t.to_string())),
//         Value::Inapplicable => Ok(None),
//         Value::Unknown => Ok(None),
//         Value::Numeric(n) => Ok(Some(format!("{n}"))),
//         Value::NumericWithUncertainty(n, u) => Ok(Some(format!("{n}({u})"))),
//     }
// }

/// Get the Numeric content of the value, if available, it also fails on NumericWithUncertainty
fn get_f32(
    value: &Value,
    context: &Context,
    column: Option<&str>,
) -> Result<Option<f32>, PDBError> {
    match value {
        Value::Numeric(num) => Ok(Some(*num)),
        Value::Inapplicable => Ok(None),
        Value::Unknown => Ok(None),
        _ => Err(PDBError::new(
            ErrorLevel::InvalidatingError,
            "Not a number",
            column.map_or(String::new(), |v| {
                format!("The '{v}' column should contain a number.")
            }),
            context.clone(),
        )),
    }
}

/// Get the Numeric content of the value, if available, as a usize
fn get_usize(
    value: &Value,
    context: &Context,
    column: Option<&str>,
) -> Result<Option<usize>, PDBError> {
    flatten_result(get_f32(value, context, column).map(|result| {
        if let Some(num) = result {
            if (0.0..usize::MAX as f32).contains(&num) && num.trunc() == num {
                Ok(Some(num as usize))
            } else {
                Err(PDBError::new(
                    ErrorLevel::InvalidatingError,
                    "Not an unsigned integer",
                    column.map_or(String::new(), |v| {
                        format!("The '{v}' column should contain an unsigned integer.")
                    }),
                    context.clone(),
                ))
            }
        } else {
            Ok(None)
        }
    }))
}

/// Get the Numeric content of the value, if available, as an isize
fn get_isize(
    value: &Value,
    context: &Context,
    column: Option<&str>,
) -> Result<Option<isize>, PDBError> {
    flatten_result(get_f32(value, context, column).map(|result| {
        if let Some(num) = result {
            if (isize::MIN as f32..isize::MAX as f32).contains(&num) && num.trunc() == num {
                Ok(Some(num as isize))
            } else {
                Err(PDBError::new(
                    ErrorLevel::InvalidatingError,
                    "Not an integer",
                    column.map_or(String::new(), |v| {
                        format!("The '{v}' column should a singed integer.")
                    }),
                    context.clone(),
                ))
            }
        } else {
            Ok(None)
        }
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::path::Path;

    #[test]
    fn test_read_cif_from_pdb() {
        let path = Path::new("data/io_test/cif/2wnb.cif");
        // Measure time
        let start = std::time::Instant::now();
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure().unwrap();
        let duration = start.elapsed();
        println!("Time elapsed in reading structure: {:?}", duration);
        // println!("{structure:?}");
        let _ = structure.to_compact();
        // println!("{:?}", compact);
    }

    #[test]
    fn test_read_cif_from_pdb_gz() {
        let path = Path::new("data/io_test/cif/2wnb.cif.gz");
        let start = std::time::Instant::now();
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure_from_gz().unwrap();
        let duration = start.elapsed();
        println!("Time elapsed in reading structure: {:?}", duration);
        // println!("{structure:?}");
        let compact = structure.to_compact();
        println!("{:?}", compact);
    }
    
    #[test]
    fn test_read_cif_from_afdb() {
        let path = Path::new("data/io_test/cif/AF-A0A4S3KKF6-F1-model_v4.cif");
        let start = std::time::Instant::now();
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure().unwrap();
        let duration = start.elapsed();
        println!("Time elapsed in reading structure: {:?}", duration);
        // println!("{structure:?}");
        let _ = structure.to_compact();
        // println!("{:?}", compact);
    }
    
    #[test]
    fn test_read_cif_from_afdb_gz() {
        let path = Path::new("data/io_test/cif/AF-A0A4S3KKF6-F1-model_v4.cif.gz");
        let start = std::time::Instant::now();
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure_from_gz().unwrap();
        let duration = start.elapsed();
        println!("Time elapsed in reading structure: {:?}", duration);
        // println!("{structure:?}");
        let compact = structure.to_compact();
        println!("{:?}", compact);
    }
}