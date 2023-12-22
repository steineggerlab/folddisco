use crate::geometry::trrosetta_subfamily::{HashCollection, HashValue};
use crate::PDBReader;

pub fn make_query(path: &String, residues: &Vec<u64>) -> Vec<u64> {
    let pdb_reader = PDBReader::from_file(path).expect("PDB file not found");
    let compact = pdb_reader.read_structure().expect("Failed to read PDB file");
    let compact = compact.to_compact();
    
    let mut hash_collection = Vec::new();
    
    // Convert residue indices to vector indices
    let mut indices = Vec::new();
    for r in residues {
        indices.push(compact.residue_serial.iter().position(|&x| x == *r).unwrap());
    }
    // Make combinations
    for i in 0..indices.len() {
        for j in i+1..indices.len() {
            let trr = compact.get_trrosetta_feature(indices[i], indices[j]).unwrap_or([0.0; 6]);
            let hash_value = HashValue::perfect_hash(
                // i as u64, j as u64, 
                trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]
            );
            hash_collection.push(hash_value.as_u64());
        }
    }
    hash_collection
}