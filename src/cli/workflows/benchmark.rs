use std::collections::HashSet;

use crate::cli::*;
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::*;

//TODO: Change this to benchmark
// usage: motifsearch benchmark -r <result.tsv> -a <answer.tsv> -f tsv
// usage: motifsearch benchmark -r <result.tsv> -a <answer.tsv> -f default
pub fn benchmark(env: AppArgs) {
    match env {
        AppArgs::Benchmark {
            result,
            answer,
        } => {
            todo!("Benchmarking is not implemented yet");
            // let result = load_path(&result);
            // let answer = load_path(&answer);
            // let result = read_tsv(&result);
            // let answer = read_tsv(&answer);
            // let result = result.iter().map(|x| x[0].clone()).collect::<Vec<_>>();
            // let answer = answer.iter().map(|x| x[0].clone()).collect::<Vec<_>>();
            // let mut set = HashSet::new();
            // let mut dedup_result = Vec::new();
            // for i in result {
            //     if !set.contains(&i) {
            //         set.insert(i.clone());
            //         dedup_result.push(i);
            //     }
            // }
            // let metric = compare_target_answer(&dedup_result, &answer);
            // println!("Precision: {}", metric.precision());
            // println!("Recall: {}", metric.recall());
            // println!("Accuracy: {}", metric.accuracy());
            // println!("F1 score: {}", metric.f1_score());
            // println!("TP: {}", metric.true_pos);
            // println!("TN: {}", metric.true_neg);
            // println!("FP: {}", metric.false_pos);
            // println!("FN: {}", metric.false_neg);
        }
        _ => {
            eprintln!("Invalid subcommand");
            std::process::exit(1);
        }
    }
}


// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_benchmark() {
//         let result = "tests/data/benchmark/result.tsv";
//         let answer = "tests/data/benchmark/answer.tsv";
//         let env = AppArgs::Benchmark {
//             result: result.to_string(),
//             answer: answer.to_string(),
//         };
//         benchmark(env);
//     }
// }
