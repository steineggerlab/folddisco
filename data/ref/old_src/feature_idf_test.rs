use std::{collections::HashSet, fs::File, io::Write, path::Path};
use folddisco::{controller::{feature::get_single_feature, map::SimpleHashMap}, prelude::*, utils::combination::CombinationIterator};
mod common;


fn get_idf_matrix_and_vector(
    query_pdb: &str, motif_residue: &Vec<(u8, u64)>, index: &str,
    hash_type: HashType, cutoff: f32, nbin_distance: usize, nbin_angle: usize,
    output_prefix: &str
) {
    // Load query pdb
    let pdb_file = File::open(query_pdb).expect("File not found");
    let pdb_reader = PDBReader::new(pdb_file);
    let compact= pdb_reader.read_structure().unwrap().to_compact();

    let mut idf_matrix: Vec<Vec<f32>>  = vec![vec![0.0; compact.num_residues]; compact.num_residues];    
    let mut motif_idf_vector: Vec<f32> = Vec::with_capacity(motif_residue.len() * motif_residue.len());
    
    // Load index
    let lookup_path = format!("{}.lookup", index);
    let offset_path = format!("{}.offset", index);
    let lookup = load_lookup_from_file(&lookup_path);
    let total_structure = lookup.len();
    let (index, mmap) = SimpleHashMap::load_from_disk(Path::new(&offset_path));
    let index_map = index.unwrap();

    // Motif 
    let motif_index: HashSet<usize> = motif_residue.iter().map(|(chain, res)| {
        compact.get_index(chain, res).unwrap()
    }).collect::<HashSet<usize>>();

    // Get feature
    let comb_iter = CombinationIterator::new(compact.num_residues);
    let mut feature = vec![0.0; 9];
    comb_iter.for_each(|(i, j)| {
        let is_feature = get_single_feature(i, j, &compact, hash_type, cutoff, &mut feature);
        if !is_feature {
            // No feature set to 0
            idf_matrix[i][j] = 0.0;
        } else {
            let hash = GeometricHash::perfect_hash(&feature, hash_type, nbin_distance, nbin_angle);
            index_map.get(&hash).map(|v| {
                let idf = (total_structure as f32 / v.1 as f32).log2();
                idf_matrix[i][j] = idf;
            });
            if motif_index.contains(&i) && motif_index.contains(&j) {
                motif_idf_vector.push(idf_matrix[i][j]);
            }
        }        
    }); 

    // Save matrix as csv
    let matrix_path = format!("{}.idf_matrix.csv", output_prefix);
    println!("Saving to {}", matrix_path);
    let mut file = File::create(matrix_path).unwrap();
    for i in 0..compact.num_residues {
        file.write_all(idf_matrix[i].iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",").as_bytes()).unwrap();
        file.write_all("\n".as_bytes()).unwrap();
    }
    
    // Save motif idf vector as csv
    let vector_path = format!("{}.motif_idf_vector.csv", output_prefix);
    let mut file = File::create(vector_path).unwrap();
    file.write_all(motif_idf_vector.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",").as_bytes()).unwrap();
    
}

#[test]
// #[ignore]
fn get_idf_matrix_vector_for_given_query() {
    // Read PDB
    let query_pdb_vec = vec![
        "/fast/hyunbin/folddisco/query/1G2F.pdb", // zinc finger
        "/fast/hyunbin/folddisco/query/4CHA.pdb", // serine protease
        "/fast/hyunbin/folddisco/query/1LAP.pdb",   
    ];
    
    let motif_residue_vec: Vec<Vec<(u8, u64)>> = vec![
        vec![(b'F', 207), (b'F', 212), (b'F', 225), (b'F', 229)],
        vec![(b'B', 57), (b'B', 102), (b'C', 195)],
        vec![(b'A', 250), (b'A', 255), (b'A', 273), (b'A', 332), (b'A', 334)],
    ];
    
    let index_vec = vec![
        "/mnt/scratch/hyunbin/human/index/human_index/h_sapiens_pdb",
        "/mnt/scratch/hyunbin/human/index/human_index/h_sapiens_folddisco",
        "/mnt/scratch/hyunbin/human/index/human_index/h_sapiens_trrosetta",
        "/mnt/scratch/hyunbin/human/index/human_index/h_sapiens_ppf",
        "/mnt/scratch/hyunbin/human/index/human_index/h_sapiens_d8a4",
        "/mnt/scratch/hyunbin/human/index/human_index/h_sapiens_d12a4",
        "/mnt/scratch/hyunbin/human/index/human_index/h_sapiens_d16a3",
        "/mnt/scratch/hyunbin/human/index/human_index/h_sapiens_d16a2",
    ];

    let hash_type_nbin_vec = vec![
        (HashType::PDBMotifSinCos, 20.0, 16, 4),
        (HashType::PDBTrRosetta, 20.0, 16, 4),
        (HashType::TrRosetta, 20.0, 16, 4),
        (HashType::PointPairFeature, 20.0, 16, 4),
        (HashType::PDBTrRosetta, 20.0, 8, 4),
        (HashType::PDBTrRosetta, 20.0, 12, 4),
        (HashType::PDBTrRosetta, 20.0, 16, 3),
        (HashType::PDBTrRosetta, 20.0, 16, 2),
    ];
    
    let motif_name_vec = vec![
        "zinc_finger",
        "serine_protease",
        "aminopeptidase",
    ];
    let index_name_vec = vec![
        "pdb",
        "folddisco",
        "trrosetta",
        "ppf",
        "d8a4",
        "d12a4",
        "d16a3",
        "d16a2",
    ];
    
    for i in 0..query_pdb_vec.len() {
        for j in 0..index_vec.len() {
            get_idf_matrix_and_vector(
                query_pdb_vec[i], &motif_residue_vec[i], index_vec[j],
                hash_type_nbin_vec[j].0, hash_type_nbin_vec[j].1, hash_type_nbin_vec[j].2, hash_type_nbin_vec[j].3,
                &format!("/mnt/scratch/hyunbin/idf_testing/{}_{}", motif_name_vec[i], index_name_vec[j])
            );
        }
    }


}