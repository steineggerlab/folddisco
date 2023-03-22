use std::fs::File;
use std::path::Path;
use std::io::{self, BufRead, BufReader};

use super::super::core::*;
use super::*;
use super::parser::*;

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
            input_type: StructureFileFormat::Unknown,
        }
    }

    /// Read from a file path
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<Self, &'static str> {
        File::open(&path)
            .map(Reader::new)
            .map_err(|_e| "Error opening file")
    }

     pub fn read_structure(&self) ->  Result<Structure, &str> {
        let reader = BufReader::new(&self.reader);
        let mut structure = Structure::new();// revise
        let mut record = (b' ',0);

        // Reading each line of PDB, parse and build atomvector.
        for (idx, line) in reader.lines().enumerate() {
            if let Ok(atomline) = line {
                match &atomline[..6] {
                    "ATOM  " => {
                        let atom = parse_line(&atomline);
                        match atom {
                            Ok(atom) => {
                                structure.update(atom, &mut record);
                            },
                            Err(e) => { // Conversion error. Jusk skip the line.
                                // If verbose, print message (NOT IMPLEMENTED)
                                // println!("Skipping line{}: {}", idx, e);
                                continue;
                            },
                        }
                    },
                    _ => continue,
                }
            } else {
                return Err("Error reading line");
            };
        }
        // println!("{structure:?}");
        Ok(structure)
     }

}

