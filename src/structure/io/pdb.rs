use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use flate2::read::GzDecoder;


use crate::structure::chain_id::ChainId;

use super::super::core::*;
use super::parser::*;
use super::*;

/// A PDB reader. Only ATOM records of the first model are read.
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    pub reader: R,
    pub input_type: StructureFileFormat,
}

impl Reader<File> {
    pub fn new(file: File) -> Self {
        Reader {
            reader: file,
            input_type: StructureFileFormat::PDB,
        }
    }

    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<Self, &'static str> {
        File::open(&path)
            .map(Reader::new)
            .map_err(|_e| "Error opening file")
    }

    /// Parse a plain-text PDB file; unparsable ATOM lines are skipped.
    pub fn read_structure(&self) -> Result<Structure, &str> {
        let reader = BufReader::new(&self.reader);
        let mut structure = Structure::new(); // revise
        let mut record = (ChainId::from_byte(b' '), 0);
        let mut model = 0;
        for (_idx, line) in reader.lines().enumerate() {
            if let Ok(atomline) = line {
                if model > 1 {
                    // First model only
                    break;
                }
                if atomline.len() < 6 {
                    continue;
                }
                match &atomline[..6] {
                    "MODEL " => {
                        model += 1;
                    }
                    "ATOM  " => {
                        let atom = parse_line(&atomline);
                        match atom {
                            Ok(atom) => {
                                structure.update(atom, &mut record);
                            }
                            Err(_e) => {
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

    /// Parse a gzipped PDB file. Unlike `read_structure`, all models are read.
    pub fn read_structure_from_gz(&self) -> Result<Structure, &str> {
        let mut decoder = GzDecoder::new(&self.reader);
        let mut binary = Vec::new();
        decoder.read_to_end(&mut binary).unwrap();
        decoder.flush().unwrap();
        drop(decoder);

        let mut structure = Structure::new();
        let mut record = (ChainId::from_byte(b' '), 0);
        
        let reader = BufReader::new(&binary[..]);
        for (_idx, line) in reader.lines().enumerate() {
            if let Ok(atomline) = line {
                match &atomline[..6] {
                    "ATOM  " => {
                        let atom = parse_line(&atomline);
                        match atom {
                            Ok(atom) => {
                                structure.update(atom, &mut record);
                            }
                            Err(_e) => {
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
        drop(binary);
        Ok(structure)
    }
    
}

#[cfg(test)]
mod tests {
    use crate::prelude::load_path;

    use super::*;
    use std::fs::File;
    
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
    fn test_get_min_max_coords() {
        let path = Path::new("data/homeobox/1akha-.pdb");
        // let path = Path::new("analysis/h_sapiens_pdb/AF-Q02817-F3-model_v4.pdb");
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure().unwrap();
        let compact = structure.to_compact();
        let min = compact.ca_vector.min_coord();
        let max = compact.ca_vector.max_coord();
        let min = min.unwrap();
        let max = max.unwrap();
        let len = compact.ca_vector.x.len();
        // Distance between min and max
        let dist = max.distance(&min);
        let diff = max.sub(&min);
        let mut bins = diff.clone();
        bins.x = (bins.x / 30.0).ceil().min(6.0);
        bins.y = (bins.y / 30.0).ceil().min(6.0);
        bins.z = (bins.z / 30.0).ceil().min(6.0);
        println!("Bins: {:?}", bins);
        // let dist = (dist.0.powi(2) + dist.1.powi(2) + dist.2.powi(2)).sqrt();
        let x_bin = diff.x / bins.x;
        let y_bin = diff.y / bins.y;
        let z_bin = diff.z / bins.z;
        println!("Min: {:?}, Max: {:?}, len: {}, dist: {:?}, x_bin: {}, y_bin: {}, z_bin: {}", min, max, len, dist, x_bin, y_bin, z_bin);
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
        let pdb_paths = load_path(dir, false);
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