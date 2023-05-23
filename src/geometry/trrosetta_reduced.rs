
use std::fmt;
use tch::{nn::{self, VarStore}, nn::Module, nn::OptimizerConfig, Kind, Reduction, Tensor, IndexOp};
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

const TRR_VEC_MEAN: [f32; 6] = [16.52766413, -1.698674698, 18.56931785, 18.56931787, 78.95890075, 78.95890075];
const TRR_VEC_STD: [f32; 6] = [8.348999265, 101.2725858, 97.42699042, 97.42699038, 38.27151932, 38.27151932];


pub fn reduce_with_vae(trr: [f32; 6], vae: &VAE) -> Vec<f32> {
    let trr = trr.iter().zip(TRR_VEC_MEAN.iter()).zip(TRR_VEC_STD.iter()).map(|((x, m), s)| (x - m) / s).collect::<Vec<f32>>();
    let trr_tensor = Tensor::of_slice(&trr).reshape(&[1, 6]).to_kind(Kind::Float);
    let encoded = vae.encode(&trr_tensor);
    let (_, quantized, _, _) = vae.vector_quantize(&encoded);

    let val1 = quantized.reshape(2).double_value(&[0]) as f32;
    let val2 = quantized.reshape(2).double_value(&[1]) as f32;
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

    pub fn vector_quantize(&self, input: &Tensor) -> (Tensor, Tensor, Tensor, Tensor) {
        let weight = self.vq.ws.t_copy();
        let distances: Tensor = input.square().sum(Kind::Float) +
            self.vq.ws.square().sum(Kind::Float) -
            2.0 * input.matmul(&weight);
        let encoding_indices = distances.argmin(1, false).unsqueeze(1);

        let mut encodings = Tensor::zeros( &[encoding_indices.size()[0], 128], (Kind::Float, input.device()) );
        // dbg!(encodings.kind());
        // dbg!(Tensor::ones_like(&encoding_indices).to_dtype(Kind::Float, true, false).kind());
        // dbg!(Tensor::ones_like(&encoding_indices).to_dtype(Kind::Float, true, false).size());
        // dbg!(encodings.size());
        // dbg!(encoding_indices.size());
        encodings.scatter_(1, &encoding_indices, &Tensor::ones_like(&encoding_indices).to_dtype(Kind::Float, true, false));

        let quantized = encodings.matmul(&self.vq.ws);
        let e_latent_loss = quantized.detach().mse_loss(input, Reduction::Mean);
        let q_latent_loss = quantized.mse_loss(&input.detach(), Reduction::Sum);
        let loss: Tensor = q_latent_loss + 0.25 * e_latent_loss;

        let quantized = input + (quantized - input).detach();
        let avg_probs = encodings.mean(Kind::Float);
        let perplexity = ((&avg_probs + 1e-10).log_() * &avg_probs).sum(Kind::Float);
        let perplexity = (perplexity * -1.0).exp();
        (loss, quantized, perplexity, encodings)
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