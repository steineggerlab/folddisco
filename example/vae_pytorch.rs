/* Variational Auto-Encoder on MNIST.
   The implementation is based on:
     https://github.com/pytorch/examples/blob/master/vae/main.py

   The 4 following dataset files can be downloaded from http://yann.lecun.com/exdb/mnist/
   These files should be extracted in the 'data' directory.
     train-images-idx3-ubyte.gz
     train-labels-idx1-ubyte.gz
     t10k-images-idx3-ubyte.gz
     t10k-labels-idx1-ubyte.gz
*/

use anyhow::Result;
use tch::{
    nn::Module,
    nn::OptimizerConfig,
    nn::{self, VarStore},
    Kind, Reduction, Tensor,
};

struct VAE {
    fc1: nn::Linear,
    fc2: nn::Linear,
    fc3: nn::Linear,
    fc4: nn::Linear,
    fc5: nn::Linear,
    mu: nn::Linear,
    var: nn::Linear,
}

impl VAE {
    fn new(vs: &nn::Path) -> Self {
        VAE {
            fc1: nn::linear(vs / "fc1", 6, 64, Default::default()),
            fc2: nn::linear(vs / "fc2", 64, 64, Default::default()),
            fc3: nn::linear(vs / "fc3", 64, 2, Default::default()),
            fc4: nn::linear(vs / "fc4", 2, 64, Default::default()),
            fc5: nn::linear(vs / "fc5", 64, 64, Default::default()),
            mu: nn::linear(vs / "mu", 64, 6, Default::default()),
            var: nn::linear(vs / "var", 64, 6, Default::default()),
        }
    }

    fn encode(&self, xs: &Tensor) -> Tensor {
        let h1 = xs.apply(&self.fc1).relu();
        let h2 = h1.apply(&self.fc2).relu();
        self.fc3.forward(&h2)
    }

    fn decode(&self, zs: &Tensor) -> (Tensor, Tensor) {
        let h1 = zs.apply(&self.fc4).relu().apply(&self.fc5).relu();
        (self.mu.forward(&h1), self.var.forward(&h1).exp())
    }

    fn forward(&self, xs: &Tensor) -> (Tensor, Tensor) {
        let z = self.encode(&xs.view([-1, 6]));
        let (mu, var) = self.decode(&z);
        // (loss, mu, var, perplexity, encodings)
        (mu, var)
    }
}

// Reconstruction + KL divergence losses summed over all elements and batch dimension.
fn loss(x: &Tensor, mu: &Tensor, var: &Tensor) -> (Tensor, Tensor) {
    // sample from the distribution having latent parameters mu, var
    let std = var.sqrt();
    let eps = Tensor::randn_like(&std);
    let z = eps * std + mu;
    let mse = z.mse_loss(x, Reduction::Mean);
    // See Appendix B from VAE paper:
    //     Kingma and Welling. Auto-Encoding Variational Bayes. ICLR, 2014
    // https://arxiv.org/abs/1312.6114
    // 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    let kld = -0.5 * (1i64 + var.log() - mu.pow_tensor_scalar(2) - var).sum(Kind::Float);
    // mse + kld
    (mse, kld)
}

fn gaussian_nll_loss(x: &Tensor, mu: &Tensor, var: &Tensor) -> Tensor {
    let loss: Tensor = 0.5 * (var.log() + ((x - mu).pow_tensor_scalar(2)) / var);
    let pi: Tensor = Tensor::of_slice(&[(2.0_f32 * std::f32::consts::PI).ln() * 0.5]);
    loss.mean(Kind::Float) + pi
}

pub fn main() -> Result<()> {
    // let device = tch::Device::cuda_if_available();
    // Get constants from command line arguments
    // let args: Vec<String> = env::args().collect();
    // let batch_size: i64 = args[1].parse().expect("Failed to parse batch size");
    // let epochs: i64 = args[2].parse().expect("Failed to parse epochs");
    // let learning_rate: f64 = args[3].parse().expect("Failed to parse learning rate");
    // let latent_dim: i64 = args[4].parse().expect("Failed to parse latent dim");
    let device = tch::Device::Cpu;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path("data/temp.vec.data.tsv")
        .expect("Failed to read data");
    // Read data
    let mut data: Vec<f32> = Vec::new();
    // Flatten data
    for result in reader.records() {
        let record = result.expect("Failed to read record");
        for field in record.iter() {
            let field: f32 = field.parse().expect("Failed to parse field");
            data.push(field);
        }
    }
    println!("{:?}", data.len());
    // Convert data to tensor
    let m = Tensor::of_slice(&data).view([-1, 6]);

    // let m = tch::vision::mnist::load_dir("data")?;
    let vs = nn::VarStore::new(device);
    let vae = VAE::new(&vs.root());
    let mut opt = nn::Adam::default().build(&vs, 1e-2)?;
    // let mut opt = nn::Sgd::default().build(&vs, 1e-2)?;

    for epoch in 1..201 {
        let mut train_loss = 0f64;
        let mut samples = 0f64;
        let mut mse_loss = 0f64;
        let mut kld_loss = 0f64;
        let mut vq_loss = 0f64;

        for x_batch in m.tensor_split(60, 0).iter() {
            // let y_batch = todo!();
            let (mu, var) = vae.forward(x_batch);
            let loss = gaussian_nll_loss(&x_batch, &mu, &var);

            opt.backward_step(&loss);
            train_loss += f64::try_from(&loss).expect("what??");
            samples += x_batch.size()[0] as f64;
        }
        println!("Epoch: {}, loss: {}", epoch, train_loss);
    }

    // Save weight
    vs.save("data/vae_weights.ot")?;

    Ok(())
}
