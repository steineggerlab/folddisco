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
use tch::{nn, nn::Module, nn::OptimizerConfig, Kind, Reduction, Tensor};

struct VAE {
    fc1: nn::Linear,
    fc21: nn::Linear,
    fc22: nn::Linear,
    fc3: nn::Linear,
    fc4: nn::Linear,
}

impl VAE {
    fn new(vs: &nn::Path) -> Self {
        VAE {
            fc1: nn::linear(vs / "fc1", 6, 64, Default::default()),
            fc21: nn::linear(vs / "fc21", 64, 2, Default::default()),
            fc22: nn::linear(vs / "fc22", 64, 2, Default::default()),
            fc3: nn::linear(vs / "fc3", 2, 64, Default::default()),
            fc4: nn::linear(vs / "fc4", 64, 6, Default::default()),
        }
    }

    fn encode(&self, xs: &Tensor) -> (Tensor, Tensor) {
        let h1 = xs.apply(&self.fc1).relu();
        (self.fc21.forward(&h1), self.fc22.forward(&h1))
    }

    fn decode(&self, zs: &Tensor) -> Tensor {
        zs.apply(&self.fc3).relu().apply(&self.fc4).sigmoid()
    }

    fn forward(&self, xs: &Tensor) -> (Tensor, Tensor, Tensor) {
        let (mu, logvar) = self.encode(&xs.view([-1, 6]));
        let std = (&logvar * 0.5).exp();
        let eps = std.randn_like();
        (self.decode(&(&mu + eps * std)), mu, logvar)
    }
}

// Reconstruction + KL divergence losses summed over all elements and batch dimension.
fn loss(recon_x: &Tensor, x: &Tensor, mu: &Tensor, logvar: &Tensor) -> Tensor {
    let mse = recon_x.mse_loss(x, Reduction::Mean);
    // See Appendix B from VAE paper:
    //     Kingma and Welling. Auto-Encoding Variational Bayes. ICLR, 2014
    // https://arxiv.org/abs/1312.6114
    // 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    let kld = -0.5 * (1i64 + logvar - mu.pow_tensor_scalar(2) - logvar.exp()).sum(Kind::Float);
    // mse + kld
    mse + kld
}

pub fn main() -> Result<()> {
    // let device = tch::Device::cuda_if_available();
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
    let mut opt = nn::Adam::default().build(&vs, 1e-1)?;
    // let mut opt = nn::Sgd::default().build(&vs, 1e-2)?;
    for epoch in 1..51 {
        let mut train_loss = 0f64;
        let mut samples = 0f64;

        for bimages in m.tensor_split(100, 0).iter(){
            // for bimage in bimages.iter(){
                let (recon_batch, mu, logvar) = vae.forward(bimages);
                let loss = loss(&recon_batch, &bimages, &mu, &logvar);
                // let loss = recon_batch.l1_loss(bimages, Reduction::Sum);
                opt.backward_step(&loss);
                train_loss += f64::try_from(&loss).expect("what??");
            // }
        }
        println!("Epoch: {}, loss: {}", epoch, train_loss);
        let s = Tensor::randn(&[1, 2], tch::kind::FLOAT_CPU).to(device);
        let s = vae.decode(&s).to(tch::Device::Cpu).view([1, 6]);
        println!("{:?}", s);
        // tch::vision::image::save(&image_matrix(&s, 8)?, format!("s_{}.png", epoch))?
    }
    Ok(())
}