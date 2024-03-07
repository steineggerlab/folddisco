use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use flate2::read::GzDecoder;

use super::super::core::*;
use super::parser::*;
use super::*;

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
            input_type: StructureFileFormat::PDB,
        }
    }

    /// Read from a file path
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<Self, &'static str> {
        File::open(&path)
            .map(Reader::new)
            .map_err(|_e| "Error opening file")
    }

    pub fn read_structure(&self) -> Result<Structure, &str> {
        let reader = BufReader::new(&self.reader);
        let mut structure = Structure::new(); // revise
        let mut record = (b' ', 0);

        // Reading each line of PDB, parse and build atomvector.
        for (idx, line) in reader.lines().enumerate() {
            if let Ok(atomline) = line {
                match &atomline[..6] {
                    "ATOM  " => {
                        let atom = parse_line(&atomline);
                        match atom {
                            Ok(atom) => {
                                structure.update(atom, &mut record);
                            }
                            Err(e) => {
                                // Conversion error. Jusk skip the line.
                                // If verbose, print message (NOT IMPLEMENTED)
                                // println!("Skipping line{}: {}", idx, e);
                                continue;
                            }
                        }
                    }
                    _ => continue,
                }
            } else {
                return Err("Error reading line");
            };
        }
        // println!("{structure:?}");
        Ok(structure)
    }
    
    pub fn read_structure_from_gz(&self) -> Result<Structure, &str> {
        let reader = BufReader::new(GzDecoder::new(&self.reader));
        
        // Parse the string
        let mut structure = Structure::new(); // revise
        let mut record = (b' ', 0);
        // Iterate over lines of the string
        for (idx, line) in reader.lines().enumerate() {
            if let Ok(atomline) = line {
                match &atomline[..6] {
                    "ATOM  " => {
                        let atom = parse_line(&atomline);
                        match atom {
                            Ok(atom) => {
                                structure.update(atom, &mut record);
                            }
                            Err(e) => {
                                // Conversion error. Jusk skip the line.
                                // If verbose, print message (NOT IMPLEMENTED)
                                // println!("Skipping line{}: {}", idx, e);
                                continue;
                            }
                        }
                    }
                    _ => continue,
                };
            } else {
                return Err("Error reading line");
            }
        }
        // println!("{structure:?}");
        // Flush the buffer
        Ok(structure)
    }
    
}

#[cfg(test)]
mod tests {
    use crate::prelude::load_path;

    use super::*;
    use std::fs::File;
    use std::io::Read;
    use std::path::Path;
    
    
    #[test]
    fn test_read_pdb() {
        let path = Path::new("data/homeobox/1akha-.pdb");
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure().unwrap();
        let compact = structure.to_compact();
        assert_eq!(compact.num_residues, 49);
    }
    
    #[test]
    fn test_read_pdb_gz() {
        let path = Path::new("data/homeobox/inner/1akha-.pdb.gz");
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure_from_gz().unwrap();
        let compact = structure.to_compact();
        assert_eq!(compact.num_residues, 49);
    }
    
    #[test]
    fn test_loading_works() {
        let dir = "data/io_test";
        let pdb_paths = load_path(dir, true);
        for pdb_path in pdb_paths {
            let file = File::open(&pdb_path).unwrap();
            let reader = Reader::new(file);
            // If the file is gzipped, use the gzipped reader
            let structure = if pdb_path.ends_with(".gz") {
                reader.read_structure_from_gz().unwrap()
            } else {
                reader.read_structure().unwrap()
            };
            let compact = structure.to_compact();
            println!("{}:{}", pdb_path, compact.num_residues);
        }
    }
}