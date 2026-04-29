use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::GzDecoder;


use super::super::core::*;
use super::parser::*;
use super::*;

fn parse_pdb_lines<I>(lines: I) -> Result<Structure, &'static str>
where
    I: Iterator<Item = io::Result<String>>,
{
    let mut structure = Structure::new();
    let mut record = (b' ', 0);
    let mut model = 0;

    // Read each line, parse ATOM records, and keep only the first model.
    for line in lines {
        let atomline = line.map_err(|_| "Error reading line")?;

        if model > 1 {
            // Current version does not support multiple models in one PDB file.
            break;
        }

        // Skip short lines (e.g., END) before record slicing.
        if atomline.len() < 6 {
            continue;
        }

        match &atomline[..6] {
            "MODEL " => {
                model += 1;
            }
            "ATOM  " => {
                if let Ok(atom) = parse_line(&atomline) {
                    structure.update(atom, &mut record);
                }
            }
            _ => continue,
        }
    }

    Ok(structure)
}

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
        parse_pdb_lines(reader.lines())
    }

    pub fn read_structure_from_gz(&self) -> Result<Structure, &str> {
        // Load whole gz file and decode into memory.
        let mut decoder = GzDecoder::new(&self.reader);
        let mut binary = Vec::new();
        decoder
            .read_to_end(&mut binary)
            .map_err(|_| "Error reading gz file")?;

        // Parse decoded lines with the same logic as plain PDB reads.
        let reader = BufReader::new(&binary[..]);
        parse_pdb_lines(reader.lines())
    }
    
}

#[cfg(test)]
mod tests {
    use crate::prelude::load_path;
    use flate2::write::GzEncoder;
    use flate2::Compression;

    use super::*;
    use std::fs::File;
    use std::io::Write;
    use std::path::Path;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn write_temp_file(prefix: &str, ext: &str, contents: &str) -> std::path::PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!("{}_{}.{}", prefix, unique, ext));
        std::fs::write(&path, contents).unwrap();
        path
    }

    fn write_temp_gz_file(prefix: &str, ext: &str, contents: &str) -> std::path::PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!("{}_{}.{}", prefix, unique, ext));
        let file = File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(contents.as_bytes()).unwrap();
        encoder.finish().unwrap();
        path
    }

    fn model_test_pdb() -> &'static str {
        "MODEL        1\n\
ATOM      1  N   ILE A  77      14.206  47.471   5.277  1.00 45.79           N  \n\
ATOM      2  CA  ILE A  77      14.689  46.123   5.703  1.00 45.28           C  \n\
ENDMDL\n\
MODEL        2\n\
ATOM      3  N   SER A  78      13.087  44.282   5.538  1.00 50.16           N  \n\
ATOM      4  CA  SER A  78      11.858  43.505   5.819  1.00 50.17           C  \n\
ENDMDL\n\
END\n"
    }

    fn end_terminated_pdb() -> &'static str {
        "ATOM      1  N   ILE A  77      14.206  47.471   5.277  1.00 45.79           N  \n\
ATOM      2  CA  ILE A  77      14.689  46.123   5.703  1.00 45.28           C  \n\
END\n"
    }
    
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

    #[test]
    fn test_read_pdb_gz_with_end_record() {
        let plain_path = write_temp_file("folddisco_pdb_end", "pdb", end_terminated_pdb());
        let gz_path = write_temp_gz_file("folddisco_pdb_end", "pdb.gz", end_terminated_pdb());

        let plain_reader = Reader::new(File::open(&plain_path).unwrap());
        let plain_compact = plain_reader.read_structure().unwrap().to_compact();

        let gz_reader = Reader::new(File::open(&gz_path).unwrap());
        let gz_compact = gz_reader.read_structure_from_gz().unwrap().to_compact();

        assert_eq!(plain_compact.num_residues, gz_compact.num_residues);

        std::fs::remove_file(plain_path).unwrap();
        std::fs::remove_file(gz_path).unwrap();
    }

    #[test]
    fn test_plain_and_gz_model_handling_match() {
        let plain_path = write_temp_file("folddisco_pdb_models", "pdb", model_test_pdb());
        let gz_path = write_temp_gz_file("folddisco_pdb_models", "pdb.gz", model_test_pdb());

        let plain_reader = Reader::new(File::open(&plain_path).unwrap());
        let plain_compact = plain_reader.read_structure().unwrap().to_compact();

        let gz_reader = Reader::new(File::open(&gz_path).unwrap());
        let gz_compact = gz_reader.read_structure_from_gz().unwrap().to_compact();

        assert_eq!(plain_compact.num_residues, gz_compact.num_residues);

        std::fs::remove_file(plain_path).unwrap();
        std::fs::remove_file(gz_path).unwrap();
    }
}