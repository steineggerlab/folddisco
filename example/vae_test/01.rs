/*
 * File: 01.rs
 * Project: vae_test
 * Created: 2023-06-29 14:37:33
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code is written as part of project "vae_test".
 * ---
 * Last Modified: 2023-07-11 17:33:18
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2023 Hyunbin Kim, All rights reserved
 */

use std::fmt;
use tch::{
    nn::Module,
    nn::OptimizerConfig,
    nn::{self, VarStore},
    IndexOp, Kind, Reduction, Tensor,
};

use motifsearch::controller::{self, Controller, GeometryHashCollector};

use motifsearch::index::builder::IndexBuilder;
use motifsearch::PDBReader;

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
    let dataset = load_homeobox_toy();
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
        for i in 0..compact.num_residues {
            for j in 0..compact.num_residues {
                counter += 1;
                if i == j {
                    continue;
                }
                let trr = compact
                    .get_trrosetta_feature2(i, j)
                    .expect("Failed to get trrosetta feature");
                // Fill trr_input_vec with trr
                trr_input_vec.push(trr);
                dbg!("{},{}, {:?}, {}", i, j, trr, trr_input_vec.len());
                // print!("{}/{} ", trr, encoded);
                    // let hash_value =
                //     HashValue::perfect_hash(trr[0], trr[1], trr[2], trr[3], trr[4], trr[5]);
                // hash_collector.collect_hash(hash_value);
                counter += 1;
            }
        }
        // Convert trr_input_vec to tensor
        // Print trr_input_vec is filled
        println!("trr_input_vec: {:?}", trr_input_vec);
        // Flatten
        let trr_input_vec = trr_input_vec.into_iter().flatten().collect::<Vec<_>>();
        let trr_input = Tensor::of_slice(&trr_input_vec).reshape(&[(N*N - N) as i64, 7]);
        let encoded = enc.forward(&trr_input);
        encoded.print();
    }
    
}
