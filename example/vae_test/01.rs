/*
 * File: 01.rs
 * Project: vae_test
 * Created: 2023-06-29 14:37:33
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     This code is written as part of project "vae_test".
 * ---
 * Last Modified: 2023-06-29 17:45:03
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



fn main() {
    // IMPORTANT: Model should be saved as safetensors
    let path = "data/encoder.safetensors";
    let loaded = tch::Tensor::read_safetensors(path);
    println!("loaded: {:?}", loaded);
    let enc = read_and_build_model();
    let temp = enc.forward(&Tensor::ones(&[1, 7], (Kind::Float, tch::Device::Cpu)));
    println!("temp: {:?}", temp);
}