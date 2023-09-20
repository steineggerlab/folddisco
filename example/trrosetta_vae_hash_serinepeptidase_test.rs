// File: trrosetta_vae_hash_serinepeptidase_test.rs
// Created: 2023-09-05 16:36:23
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

use std::fmt;
use tch::{
    nn::Module,
    nn::OptimizerConfig,
    nn::{self, VarStore},
    IndexOp, Kind, Reduction, Tensor,
};

use motifsearch::{geometry::two_float::{HashValue, HashCollection}, index::IndexTablePrinter};
use motifsearch::controller::{self, Controller, GeometryHashCollector};
use motifsearch::index::builder::IndexBuilder;
use motifsearch::PDBReader;

// Use hashmap
use std::collections::HashMap;

struct Encoder {
    fc1: nn::Linear,
    fc2: nn::Linear,
    fc3: nn::Linear,
}

impl Encoder {
    fn forward(&self, x: &Tensor) -> Tensor {
        x.apply(&self.fc1).relu().apply(&self.fc2).relu().apply(&self.fc3)
    }
}

fn get_tensor_from_vec(vec: &Vec<(String, Tensor)>, name: &str) -> Tensor {
    for (n, t) in vec {
        if n == name {
            let mut out = Tensor::empty(&t.size(), (Kind::Float, t.device()));
            out.copy_(&t);
            return out;
        }
    }
    panic!("No tensor named {}", name);
}

fn read_and_build_model() -> Encoder {
    let path = "data/encoder.safetensors";
    let loaded = tch::Tensor::read_safetensors(path).expect("Failed to read safetensors");

    let mut vs = VarStore::new(tch::Device::Cpu);
    let mut fc1 = nn::linear(&vs.root() / "fc1", 7, 64, Default::default());
    let mut fc2 = nn::linear(&vs.root() / "fc2", 64, 64, Default::default());
    let mut fc3 = nn::linear(&vs.root() / "fc3", 64, 2, Default::default());
    
    let fc1_weight = get_tensor_from_vec(&loaded, "0.weight");
    let fc1_bias = get_tensor_from_vec(&loaded, "0.bias");
    let fc2_weight = get_tensor_from_vec(&loaded, "2.weight");
    let fc2_bias = get_tensor_from_vec(&loaded, "2.bias");
    let fc3_weight = get_tensor_from_vec(&loaded, "4.weight");
    let fc3_bias = get_tensor_from_vec(&loaded, "4.bias");

    fc1.ws = fc1_weight;
    fc1.bs = Some(fc1_bias);
    fc2.ws = fc2_weight;
    fc2.bs = Some(fc2_bias);
    fc3.ws = fc3_weight;
    fc3.bs = Some(fc3_bias);
    
    Encoder {
        fc1: fc1,
        fc2: fc2,
        fc3: fc3,
    }
}

fn load_homeobox_toy() -> Vec<String> {
    vec![
        "data/homeobox/1akha-.pdb".to_string(),
        "data/homeobox/1b72a-.pdb".to_string(),
        "data/homeobox/1b72b-.pdb".to_string(),
        "data/homeobox/1ba5--.pdb".to_string(),
    ]
}
fn load_path(dir: &str) -> Vec<String> {
    // Load all pdbs in given path
    let mut pdb_paths = Vec::new();
    let paths = std::fs::read_dir(dir).expect("Unable to read pdb directory");
    for path in paths {
        let path = path.expect("Unable to read path");
        let path = path.path();
        let path = path.to_str().expect("Unable to convert path to string");
        // If the path is a pdb file, add it to the list
        if path.ends_with(".pdb") {
            pdb_paths.push(path.to_string());
        }
    }
    pdb_paths
}

fn main() {
    // IMPORTANT: Model should be saved as safetensors
    // Load model
    let path = "data/encoder.safetensors";
    let loaded = tch::Tensor::read_safetensors(path);
    let enc = read_and_build_model();
    // Test if encoder works
    // let temp = enc.forward(&Tensor::ones(&[1, 7], (Kind::Float, tch::Device::Cpu)));
    // println!("temp: {:?}", temp);

    // Load dataset
    // let dataset = load_homeobox_toy();
    let dataset = load_path("data/serine_peptidases_filtered");
    let mut controller = Controller::new(dataset.clone());
    controller.fill_numeric_id_vec();
    controller.path_vec = dataset.clone();
    for pdb_path in dataset {
        // Start measure time
        let start = std::time::Instant::now();
        let pdb_reader = PDBReader::from_file(&pdb_path).expect("Failed to read PDB file");
        let structure = pdb_reader.read_structure().expect("Failed to read structure");
        let compact = structure.to_compact();

        let mut hash_collector = GeometryHashCollector::new();

        let N: usize = compact.num_residues.try_into().unwrap();
        let mut trr_input_vec = Vec::with_capacity(N * N - N);
        let mut counter: usize = 0;
        let mut inner_res_pair_vec =  Vec::with_capacity(N * N - N);
        for i in 0..compact.num_residues {
            for j in 0..compact.num_residues {
                counter += 1;
                if i == j {
                    continue;
                }
                let trr = compact
                    .get_trrosetta_feature2(i, j);
                if trr.is_none() {
                    let empty: [f32; 7] = [0.0; 7];
                    let trr = empty;
                    trr_input_vec.push(trr);
                    // Fill inner_res_pair_vec with (i, j)
                    let inner_res_pair = (
                        compact.residue_serial[i], compact.residue_serial[j],
                        compact.residue_name[i], compact.residue_name[j],
                    );
                    inner_res_pair_vec.push(inner_res_pair);
                    counter += 1;
                    continue;
                }
                // Fill trr_input_vec with trr
                trr_input_vec.push(trr.unwrap());
                
                // Fill inner_res_pair_vec with (i, j)
                let inner_res_pair = (
                    compact.residue_serial[i], compact.residue_serial[j],
                    compact.residue_name[i], compact.residue_name[j],
                );
                inner_res_pair_vec.push(inner_res_pair);
                // dbg!("{},{}, {:?}, {}", i, j, trr, trr_input_vec.len());
                // print!("{}/{} ", trr, encoded);
                    // let hash_value =
                //     HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);
                // hash_collector.collect_hash(hash_value);
                counter += 1;
            }
        }
        // Convert trr_input_vec to tensor
        // Print trr_input_vec is filled
        // println!("trr_input_vec: {:?}", trr_input_vec);
        // Flatten
        let trr_input_vec = trr_input_vec.into_iter().flatten().collect::<Vec<_>>();
        let trr_input = Tensor::of_slice(&trr_input_vec).reshape(&[(N*N - N) as i64, 7]);
        let encoded = enc.forward(&trr_input);
        // encoded.print();
        // Convert torch tensor to vector
        let encoded_vec: Vec<Vec<f32>> = TryFrom::try_from(encoded).unwrap();
        // 2023-07-11 17:36:44
        let mut hash_vec = Vec::with_capacity(encoded_vec.len());
        for i in 0..encoded_vec.len() {
            let hash_value = HashValue::perfect_hash(encoded_vec[i][0], encoded_vec[i][1]);
            hash_vec.push(hash_value);
        }
        controller.hash_new_collection_vec.push(hash_vec);
        controller.res_pair_vec.push(inner_res_pair_vec);
    }
    println!("Numeric id length {:?}", controller.numeric_id_vec.len());
    println!("Hashvec length {:?}", controller.hash_new_collection_vec.len());
    let index_builder = IndexBuilder::new();
    let index_table = index_builder.concat(
        &controller.numeric_id_vec,
        &controller.hash_new_collection_vec,
    );
    let table_printer = IndexTablePrinter::Debug;
    table_printer.print(&index_table, "data/serine_vae_hash_table_rounded.tsv");
    println!("Checking");
    let mut serine_filter: HashMap<String, Vec<u64>> = HashMap::new();
    serine_filter.insert("1aq2.pdb".to_string(), vec![250, 232, 269]);
    serine_filter.insert("1wab.pdb".to_string(), vec![47, 195, 192]);
    serine_filter.insert("1sc9.pdb".to_string(), vec![80, 235, 207]);
    serine_filter.insert("2o7r.pdb".to_string(), vec![169, 306, 276]);
    serine_filter.insert("1bs9.pdb".to_string(), vec![90, 187, 175]);
    serine_filter.insert("1ju3.pdb".to_string(), vec![117, 287, 259]);
    serine_filter.insert("1uk7.pdb".to_string(), vec![34, 252, 224]);
    serine_filter.insert("1okg.pdb".to_string(), vec![255, 75, 61]);
    serine_filter.insert("1qfm.pdb".to_string(), vec![554, 680, 641]);
    controller.save_hash_per_pair("data/serine_vae_per_pair_hash.tsv");
    controller.save_filtered_hash_pair("data/serine_vae_per_pair_hash_filtered.tsv", &serine_filter);
    println!("DONE");
}
