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

use std::any::Any;

use anyhow::Result;
use rand::distributions::uniform::SampleBorrow;
use tch::{nn, nn::Module, nn::OptimizerConfig, Kind, Reduction, Tensor};

struct VAE {
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
    fn new(vs: &nn::Path) -> Self {
        VAE {
            fc1: nn::linear(vs / "fc1", 6, 64, Default::default()),
            fc2: nn::linear(vs / "fc2", 64, 32, Default::default()),
            fc3: nn::linear(vs / "fc3", 32, 2, Default::default()),
            vq: nn::linear(vs / "vq", 2, 32, Default::default()),
            fc4: nn::linear(vs / "fc4", 2, 32, Default::default()),
            fc5: nn::linear(vs / "fc5", 32, 64, Default::default()),
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

    fn vector_quantize(&self, input: &Tensor) -> (Tensor, Tensor, Tensor, Tensor) {
        let weight = self.vq.ws.t_copy();
        let distances: Tensor = input.square().sum(Kind::Float) +
            self.vq.ws.square().sum(Kind::Float) -
            2.0 * input.matmul(&weight);
        let encoding_indices = distances.argmin(1, false).unsqueeze(1);
        
        let mut encodings = Tensor::zeros( &[encoding_indices.size()[0], 32], (Kind::Float, input.device()) );
        // dbg!(encodings.kind());
        // dbg!(Tensor::ones_like(&encoding_indices).to_dtype(Kind::Float, true, false).kind());
        // dbg!(Tensor::ones_like(&encoding_indices).to_dtype(Kind::Float, true, false).size());
        // dbg!(encodings.size());
        // dbg!(encoding_indices.size());
        encodings.scatter_(1, &encoding_indices, &Tensor::ones_like(&encoding_indices).to_dtype(Kind::Float, true, false));
        
        let quantized = encodings.matmul(&self.vq.ws);
        
        let e_latent_loss = quantized.detach().mse_loss(input, Reduction::Mean);
        let q_latent_loss = quantized.mse_loss(&input.detach(), Reduction::Mean);
        let loss = q_latent_loss + 0.25 * e_latent_loss;
        let quantized = input + (quantized - input).detach();
        let avg_probs = encodings.mean(Kind::Float);
        let perplexity = ((&avg_probs + 1e-10).log_() * &avg_probs).sum(Kind::Float);
        let perplexity = (perplexity * -1.0).exp();
        (loss, quantized, perplexity, encodings)
    }

    fn forward(&self, xs: &Tensor) -> (Tensor, Tensor, Tensor, Tensor, Tensor) {
        let z = self.encode(&xs.view([-1, 6]));
        let (loss, quantized, perplexity, encodings) = self.vector_quantize(&z);
        let (mu, var) = self.decode(&quantized);
        (loss, mu, var, perplexity, encodings)
    }
}

// Reconstruction + KL divergence losses summed over all elements and batch dimension.
fn loss(recon_x: &Tensor, x: &Tensor, mu: &Tensor, var: &Tensor) -> Tensor {
    let mse = recon_x.mse_loss(x, Reduction::Mean);
    // See Appendix B from VAE paper:
    //     Kingma and Welling. Auto-Encoding Variational Bayes. ICLR, 2014
    // https://arxiv.org/abs/1312.6114
    // 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    let kld = -0.5 * (1i64 + var.log() - mu.pow_tensor_scalar(2) - var).sum(Kind::Float);
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
                let (recon_batch, mu, var, perp, emb) = vae.forward(bimages);
                let loss = loss(&recon_batch, &bimages, &mu, &var);
                // let loss = recon_batch.l1_loss(bimages, Reduction::Sum);
                opt.backward_step(&loss);
                train_loss += f64::try_from(&loss).expect("what??");
                samples += bimages.size()[0] as f64;
                if epoch == 50 {
                    let (recon_batch, mu, var, perp, emb) = vae.forward(&m);
                    emb.save("emb.pt")?;
                    perp.save("perp.pt")?;
                    recon_batch.save("recon_batch.pt")?;
                    mu.save("mu.pt")?;
                    var.save("var.pt")?;
                }

            // }
        }
        println!("Epoch: {}, loss: {}", epoch, train_loss / samples);


        // tch::vision::image::save(&image_matrix(&s, 8)?, format!("s_{}.png", epoch))?
    }
    Ok(())
}