use std::fs;
use std::io;
use std::path::Path;

use super::super::core::Structure;
use super::StructureFileFormat;

/// A PDB reader
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    /// The underlying reader
    reader: R,
    /// The current line number
    line_number: usize,
    /// The current line
    line: String,
    ///
    input_type: StructureFileFormat,
}

impl Reader<fs::File> {
    pub fn new(file: fs::File) -> Self {
        Reader {
            reader: file,
            line_number: 0,
            line: String::new(),
            input_type: StructureFileFormat::Unknown,
        }
    }

    /// Read from a file path
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> io::Result<Self> {
        fs::File::open(&path)
            .map(Reader::new)
            .map_err(|e| io::Error::new(e.kind(), format!("{}: {:?}", e, path)))
    }

    /// Return structure.
    pub fn parse(&self) -> Result<Structure, String> {
        todo!()
    }

}

