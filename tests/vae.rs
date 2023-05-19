// // VQ-VQE written with dfdx crate

// // Import

// use dfdx::nn::modules::{BatchNorm1D, Linear, Module, ModuleVisitor, ReLU, TensorCollection};
// use dfdx::optim::{Momentum, Sgd, SgdConfig, WeightDecay};
// use dfdx::shapes::{Dtype, Rank1, Rank2};
// use dfdx::tensor::{AutoDevice, SampleTensor, Tape, Tensor, Trace};
// use dfdx::tensor_ops::Device;

// //
// struct VAE<
//     const IN: usize,
//     const HIDDEN1: usize,
//     const HIDDEN2: usize,
//     const OUT: usize,
//     E: Dtype,
//     D: Device<E>,
// > {
//     encoder: (
//         (Linear<IN, HIDDEN1, E, D>, BatchNorm1D<HIDDEN1, E, D>, ReLU),
//         (Linear<HIDDEN1, HIDDEN2, E, D>, BatchNorm1D<HIDDEN2, E, D>, ReLU),
//     ),
//     mean: Linear<HIDDEN2, OUT, E, D>,
//     logvar: Linear<HIDDEN2, OUT, E, D>,
//     decoder: (
//         (Linear<OUT, HIDDEN2, E, D>, BatchNorm1D<HIDDEN2, E, D>, ReLU),
//         ( Linear<HIDDEN2, HIDDEN1, E, D>, BatchNorm1D<HIDDEN1, E, D>, ReLU),
//         Linear<HIDDEN1, IN, E, D>,
//     ),
// }

// impl<
//         const IN: usize,
//         const HIDDEN1: usize,
//         const HIDDEN2: usize,
//         const OUT: usize,
//         E: Dtype,
//         D: Device<E>,
//     > TensorCollection<E, D> for VAE<IN, HIDDEN1, HIDDEN2, OUT, E, D>
// {
//     type To<E2: Dtype, D2: Device<E2>> = VAE<IN, HIDDEN1, HIDDEN2, OUT, E2, D2>;

//     fn iter_tensors<V: ModuleVisitor<Self, E, D>>(
//             visitor: &mut V,
//         ) -> Result<Option<Self::To<V::E2, V::D2>>, V::Err> {
//         visitor.visit_fields(
//             // Define name of each field and how to access it,
//             // using ModuleField for Modules, and TensorField for Tensors
//             (
//                 Self::module("encoder", |s| &s.encoder, |s| &mut s.encoder),
//                 Self::module("mean", |s| &s.mean, |s| &mut s.mean),
//                 Self::module("logvar", |s| &s.logvar, |s| &mut s.logvar),
//                 Self::module("decoder", |s| &s.decoder, |s| &mut s.decoder)
//             ),
//             |(encoder, mean, logvar, decoder)| VAE {
//                 encoder,
//                 mean,
//                 logvar,
//                 decoder,
//             },
//         )
//     }




// }
