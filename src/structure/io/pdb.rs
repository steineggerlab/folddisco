
use std::fs;
use std::io;
use std::path::Path;

/// A PDB reader
#[derive(Debug)]
pub struct PDBReader<R: io::Read> {
    /// The underlying reader
    reader: R,
    /// The current line number
    line_number: usize,
    /// The current line
    line: String,
}

impl Reader<fs::File> {
    /// Read from a file path
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> io::Result<Self> {
        fs::File::open(&path)
            .map(Reader::new)
            .with_context(|| format!("Failed to open file {:#?}", path))
    }
}