// use motifsearch::*;

// type EncoderLayer = (
//     (Linear<6, 32>, ReLU),
//     (Linear<32, 64>, ReLU),
//     (Linear<64, 6>, ReLU),
// );

// type DecoderLayer = (
//     (Linear<6, 64>, ReLU),
//     (Linear<64, 32>, ReLU),
//     (Linear<32, 6>, Sigmoid),
// );

use dfdx::optim::{Momentum, Sgd, SgdConfig, WeightDecay};
use dfdx::prelude::*;

fn main() {
    // Log start time
    let start = std::time::Instant::now();
    type EncoderLayer = (
        (Linear<6, 32>, BatchNorm1D<32>, ReLU),
        (Linear<32, 32>, BatchNorm1D<32>, ReLU),
        Linear<32, 2>,
    );
    type DecoderLayer = (
        (Linear<2, 32>, BatchNorm1D<32>, ReLU),
        (Linear<32, 32>, BatchNorm1D<32>, ReLU),
        Linear<32, 6>,
    );

    type Mlp = (EncoderLayer, DecoderLayer);

    let dev = AutoDevice::default();
    let mut mlp = dev.build_module::<Mlp, f32>();
    let mut grads = mlp.alloc_grads();
    let mut sgd = Sgd::new(
        &mlp,
        SgdConfig {
            lr: 1e-2,
            momentum: Some(Momentum::Nesterov(0.9)),
            weight_decay: Some(WeightDecay::L2(1e-4)),
        },
    );

    let x: Tensor<Rank2<4, 6>, f32, _> = dev.sample_normal();
    // let y: Tensor<Rank2<4, 2>, f32, _> = dev.sample_normal();

    for i in 0..50 {
        let prediction = mlp.forward_mut(x.trace(grads));
        let loss = mse_loss(prediction, x.clone());
        println!("Epoch {:03}: {:?}", i, loss.array());
        grads = loss.backward();
        sgd.update(&mut mlp, &grads).expect("SGD failed");
        mlp.zero_grads(&mut grads);
    }
    let end = std::time::Instant::now();
    println!("{:?}", mlp.forward(x).array());
    println!("Time elapsed: {:?}", end - start);
    mlp.save("test.npz").expect("Failed to save");
}
