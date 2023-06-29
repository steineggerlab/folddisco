use std::fmt;
use tch::{
    nn::Module,
    nn::OptimizerConfig,
    nn::{self, VarStore},
    IndexOp, Kind, Reduction, Tensor,
};
#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u64);

impl HashValue {
    pub fn from_u64(hashvalue: u64) -> Self {
        HashValue(hashvalue)
    }

    pub fn perfect_hash(val1: f32, val2: f32) -> Self {
        // Save two f32 as is in a u64
        let hashvalue = (val1.to_bits() as u64) << 32 | (val2.to_bits() as u64);
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> [f32; 2] {
        let val1 = f32::from_bits((self.0 >> 32) as u32);
        let val2 = f32::from_bits((self.0 & 0x00000000FFFFFFFF) as u32);
        [val1, val2]
    }
}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let values = self.reverse_hash();
        write!(f, "HashValue({}), values={:?}", self.0, values)
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let val = self.reverse_hash();
        write!(f, "{}\t{}\t{}", self.0, val[0], val[1])
    }
}

pub type HashCollection = Vec<HashValue>;

pub fn reduce_with_vae(trr: [f32; 7], vae: &VAE) -> Vec<f32> {
    let trr_tensor = Tensor::of_slice(&trr).reshape(&[1, 7]).to_kind(Kind::Float);
    let encoded = vae.encode(&trr_tensor);
    vec![val1, val2]
}

pub fn load_vae() -> VAE {
    let mut vs = VarStore::new(tch::Device::Cpu);
    vs.load("data/tr_vae.ot");
    let vae = VAE::new(&vs.root());
    vae
}

pub struct VAE {
    fc1: nn::Linear,
    fc2: nn::Linear,
    fc3: nn::Linear,
    vq: nn::Linear,
    fc4: nn::Linear,
    fc5: nn::Linear,
    mu: nn::Linear,
    var: nn::Linear,
}

impl VAE {
    pub fn new(vs: &nn::Path) -> Self {
        VAE {
            fc1: nn::linear(vs / "fc1", 6, 64, Default::default()),
            fc2: nn::linear(vs / "fc2", 64, 128, Default::default()),
            fc3: nn::linear(vs / "fc3", 128, 2, Default::default()),
            vq: nn::linear(vs / "vq", 2, 128, Default::default()),
            fc4: nn::linear(vs / "fc4", 2, 128, Default::default()),
            fc5: nn::linear(vs / "fc5", 128, 64, Default::default()),
            mu: nn::linear(vs / "mu", 64, 6, Default::default()),
            var: nn::linear(vs / "var", 64, 6, Default::default()),
        }
    }

    pub fn encode(&self, xs: &Tensor) -> Tensor {
        let h1 = xs.apply(&self.fc1).relu();
        let h2 = h1.apply(&self.fc2).relu();
        self.fc3.forward(&h2)
    }

    pub fn decode(&self, zs: &Tensor) -> (Tensor, Tensor) {
        let h1 = zs.apply(&self.fc4).relu().apply(&self.fc5).relu();
        (self.mu.forward(&h1), self.var.forward(&h1).exp())
    }

    pub fn forward(&self, xs: &Tensor) -> (Tensor, Tensor, Tensor, Tensor, Tensor) {
        let z = self.encode(&xs.view([-1, 6]));
        let (loss, quantized, perplexity, encodings) = self.vector_quantize(&z);
        let (mu, var) = self.decode(&quantized);
        (loss, mu, var, perplexity, encodings)
    }
}

pub fn normalize_angle_degree(val: f32, min: f32, max: f32) -> f32 {
    (val - min) / (max - min)
}
