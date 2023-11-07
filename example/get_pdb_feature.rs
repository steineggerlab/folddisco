use motifsearch::geometry::trrosetta_subfamily::{HashCollection, HashValue};
use motifsearch::PDBReader;

pub fn make_query(path: &String, residues: &Vec<u64>) -> HashCollection {
    let pdb_reader = PDBReader::from_file(path).expect("PDB file not found");
    let compact = pdb_reader.read_structure().expect("Failed to read PDB file");
    let compact = compact.to_compact();
    
    let mut hash_collection = HashCollection::new();
    
    // Convert residue indices to vector indices
    let mut indices = Vec::new();
    for r in residues {
        indices.push(compact.residue_serial.iter().position(|&x| x == *r).unwrap());
    }
    // Make combinations
    for i in 0..indices.len() {
        for j in i+1..indices.len() {
            let trr = compact.get_trrosetta_feature(indices[i], indices[j]).unwrap_or([0.0; 6]);
            let hash_value = HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);
            hash_collection.push(hash_value);
        }
    }
    hash_collection
}


pub fn main() {
    let list: Vec<(String, Vec<u64>)> = vec![
        ("analysis/raw_ecoli/AF-P76176-F1-model_v4.pdb".to_string(), vec![84, 145, 223]),
        // ("analysis/raw_ecoli/AF-Q7BSW5-F1-model_v4.pdb".to_string(), vec![127, 156, 263]),
    ];
    for (path, residues) in list {
        let hash_collection = make_query(&path, &residues);
        println!("Path: {:?}", path);
        println!("Residues: {:?}", residues);
        println!("Hash collection: {:?}", hash_collection);
    }
}