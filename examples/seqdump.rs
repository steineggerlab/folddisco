// Dump one-letter sequences from a Foldcomp DB, so a benchmark answer set can be
// built from sequence patterns instead of external annotation. See
// scripts/build_serine_answer.py for the consumer.
//
//   cargo run --release --example seqdump -- index/h_sapiens > h_sapiens_seq.tsv
//
// Output: <structure name>\t<one-letter sequence>, one line per DB entry.

use folddisco::structure::io::fcz::FoldcompDbReader;

fn three_to_one(r: &[u8; 3]) -> char {
    match r {
        b"ALA" => 'A', b"ARG" => 'R', b"ASN" => 'N', b"ASP" => 'D', b"CYS" => 'C',
        b"GLN" => 'Q', b"GLU" => 'E', b"GLY" => 'G', b"HIS" => 'H', b"ILE" => 'I',
        b"LEU" => 'L', b"LYS" => 'K', b"MET" => 'M', b"PHE" => 'F', b"PRO" => 'P',
        b"SER" => 'S', b"THR" => 'T', b"TRP" => 'W', b"TYR" => 'Y', b"VAL" => 'V',
        b"MSE" => 'M', b"SEC" => 'U', b"PYL" => 'O',
        _ => 'X',
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let reader = FoldcompDbReader::new(&args[1]);
    let paths = reader.get_paths();          // needs the db-key sort new() applies
    for name in paths.iter() {
        match reader.read_single_structure(name) {
            Ok(s) => {
                let av = &s.atom_vector;
                let mut seq = String::with_capacity(s.num_residues);
                for i in 0..av.atom_name.len() {
                    if &av.atom_name[i] == b" CA " {
                        seq.push(three_to_one(&av.res_name[i]));
                    }
                }
                println!("{}\t{}", name, seq);
            }
            Err(e) => eprintln!("ERR\t{}\t{}", name, e),
        }
    }
}
