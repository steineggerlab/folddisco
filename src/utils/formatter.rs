///! TSV Formatter Utility
///! This module provides a generic TSV formatter.

use std::io::{self, Write};
use std::sync::Arc;

pub const DEFAULT_FLOAT_PRECISION: usize = 4;

/// Runtime value with embedded formatting information.
#[derive(Clone, Debug)]
pub enum Value {
    Int(i64),
    Uint(u64),
    Float(f32, usize), // value, precision
    ScientificFloat(f64, usize), // value, precision. Print in scientific notation
    Str(String),
    Bool(bool),
    Float3DMatrix([[f32; 3]; 3], usize, &'static str), // matrix, precision, separator
    Float3DVector([f32; 3], usize, &'static str), // vector, precision, separator
    FloatVector(Vec<f32>, usize, &'static str), // vector, precision, separator
}

impl From<i64> for Value {
    fn from(v: i64) -> Self {
        Value::Int(v)
    }
}

impl From<u64> for Value {
    fn from(v: u64) -> Self {
        Value::Uint(v)
    }
}

impl From<f32> for Value {
    fn from(v: f32) -> Self {
        Value::Float(v, DEFAULT_FLOAT_PRECISION) // default precision
    }
}

impl From<f64> for Value {
    fn from(v: f64) -> Self {
        Value::ScientificFloat(v, DEFAULT_FLOAT_PRECISION) // default precision
    }
}

impl From<String> for Value {
    fn from(v: String) -> Self {
        Value::Str(v)
    }
}

impl From<&str> for Value {
    fn from(v: &str) -> Self {
        Value::Str(v.to_string())
    }
}

impl From<bool> for Value {
    fn from(v: bool) -> Self {
        Value::Bool(v)
    }
}

impl From<[[f32; 3]; 3]> for Value {
    fn from(v: [[f32; 3]; 3]) -> Self {
        Value::Float3DMatrix(v, DEFAULT_FLOAT_PRECISION, ",") // default precision and separator
    }
}

impl From<[f32; 3]> for Value {
    fn from(v: [f32; 3]) -> Self {
        Value::Float3DVector(v, DEFAULT_FLOAT_PRECISION, ",") // default precision and separator
    }
}

impl From<Vec<f32>> for Value {
    fn from(v: Vec<f32>) -> Self {
        Value::FloatVector(v, DEFAULT_FLOAT_PRECISION, ",") // default precision and separator
    }
}

/// A column definition that combines metadata and extraction logic.
pub struct Column<R> {
    /// Parameter name / CLI key (e.g. "score", "rmsd", "query_id").
    pub key: &'static str,
    /// Description of what this column represents.
    pub description: &'static str,
    /// Function to extract the value from a record.
    pub extractor: Arc<dyn Fn(&R) -> Value + Send + Sync>,
}

impl<R> Clone for Column<R> {
    fn clone(&self) -> Self {
        Column {
            key: self.key,
            description: self.description,
            extractor: Arc::clone(&self.extractor),
        }
    }
}

impl<R> Column<R> {
    pub fn new<F>(key: &'static str, description: &'static str, extractor: F) -> Self
    where
        F: Fn(&R) -> Value + Send + Sync + 'static,
    {
        Column {
            key,
            description,
            extractor: Arc::new(extractor),
        }
    }
}

/// Generic TSV formatter that can be instantiated for any record type `R`.
pub struct TsvFormatter<R> {
    pub columns: Vec<Column<R>>,
}

impl<R> TsvFormatter<R> {
    pub fn new(columns: Vec<Column<R>>) -> Self {
        TsvFormatter { columns }
    }

    /// Write header line (column keys) to the writer.
    pub fn write_header<W: Write>(&self, mut w: W) -> io::Result<()> {
        for (i, col) in self.columns.iter().enumerate() {
            if i > 0 {
                write!(w, "\t")?;
            }
            write!(w, "{}", col.key)?;
        }
        writeln!(w)
    }

    /// Write a single record as TSV.
    pub fn write_record<W: Write>(&self, mut w: W, record: &R) -> io::Result<()> {
        for (i, col) in self.columns.iter().enumerate() {
            if i > 0 {
                write!(w, "\t")?;
            }
            let value = (col.extractor)(record);
            Self::write_value(&mut w, &value)?;
        }
        writeln!(w)
    }

    fn write_value<W: Write>(w: &mut W, value: &Value) -> io::Result<()> {
        match value {
            Value::Int(v) => write!(w, "{}", v),
            Value::Uint(v) => write!(w, "{}", v),
            Value::Float(v, precision) => write!(w, "{:.1$}", v, precision),
            Value::ScientificFloat(v, precision) => write!(w, "{:.1$e}", v, precision),
            Value::Str(s) => write!(w, "{}", escape_tsv(s)),
            Value::Bool(b) => write!(w, "{}", if *b { "1" } else { "0" }),
            Value::Float3DMatrix(m, precision, separator) => {
                for (i, row) in m.iter().enumerate() {
                    if i > 0 { write!(w, "{}", separator)?; }
                    for (j, val) in row.iter().enumerate() {
                        if j > 0 { write!(w, "{}", separator)?; }
                        write!(w, "{:.1$}", val, precision)?;
                    }
                }
                Ok(())
            }
            Value::Float3DVector(v, precision, separator) => {
                for (i, val) in v.iter().enumerate() {
                    if i > 0 { write!(w, "{}", separator)?; }
                    write!(w, "{:.1$}", val, precision)?;
                }
                Ok(())
            }
            Value::FloatVector(vec, precision, separator) => {
                for (i, val) in vec.iter().enumerate() {
                    if i > 0 { write!(w, "{}", separator)?; }
                    write!(w, "{:.1$}", val, precision)?;
                }
                Ok(())
            }
        }
    }
}

/// Minimal TSV escaping (tabs and newlines).
fn escape_tsv(s: &str) -> String {
    s.replace('\t', " ").replace('\n', " ")
}

// Tests
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_tsv_formatter() {
        struct Record {
            id: i64,
            score: f32,
            score_scientific: f64,
            name: String,
            active: bool,
        }
        let columns = vec![
            Column::new("id", "Record ID", |r: &Record| Value::from(r.id)),
            Column::new("score", "Score", |r: &Record| Value::from(r.score)),
            Column::new("score_sci", "Score Scientific", |r: &Record| Value::ScientificFloat(r.score_scientific, 4)),
            Column::new("name", "Name", |r: &Record| Value::from(r.name.clone())),
            Column::new("active", "Active Status", |r: &Record| Value::from(r.active)),
        ];
        let formatter = TsvFormatter::new(columns);
        let record = Record {
            id: 42,
            score: 3.14159,
            score_scientific: 0.0000000314159,
            name: "Test\tName".to_string(),
            active: true,
        };
        let mut output = Vec::new();
        formatter.write_header(&mut output).unwrap();
        formatter.write_record(&mut output, &record).unwrap();
        let result = String::from_utf8(output).unwrap();
        let expected = "id\tscore\tscore_sci\tname\tactive\n42\t3.1416\t3.1500e-08\tTest Name\t1\n";
        println!("Result:\n{}", result);
        assert_eq!(result, expected);
        
        let record_vec = vec![
            Record {
                id: 1,
                score: 2.71828,
                score_scientific: 0.0000000271828,
                name: "Alice".to_string(),
                active: false,
            },
            Record {
                id: 2,
                score: 1.61803,
                score_scientific: 0.0000000161803,
                name: "Bob".to_string(),
                active: true,
            },
        ];

        let columns_new = vec![
            Column::new("name", "Name", |r: &Record| Value::from(r.name.clone())),
            Column::new("active", "Active Status", |r: &Record| Value::from(r.active)),
            Column::new("score", "Score", |r: &Record| Value::from(r.score)),
            Column::new("score_sci", "Score Scientific", |r: &Record| Value::ScientificFloat(r.score_scientific, 4)),
        ];
        let formatter_new = TsvFormatter::new(columns_new);
        let mut output_new = Vec::new();
        formatter_new.write_header(&mut output_new).unwrap();
        for rec in &record_vec {
            formatter_new.write_record(&mut output_new, rec).unwrap();
        }
        let result_new = String::from_utf8(output_new).unwrap();
        println!("Result with new columns:\n{}", result_new);
        let expected_new = "name\tactive\tscore\tscore_sci\nAlice\t0\t2.7183\t2.7183e-08\nBob\t1\t1.6180\t1.6180e-08\n";
        assert_eq!(result_new, expected_new);
    }
}   