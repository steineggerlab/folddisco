use std::collections::HashSet;
use std::io::BufRead;

use crate::cli::*;
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::*;

use crate::cli::config::read_index_config_from_file;

// usage: folddisco benchmark -r <result.tsv> -a <answer.tsv> -i <index> -f tsv
// usage: folddisco benchmark -r <result.tsv> -a <answer.tsv> -i <index> -f default
pub fn benchmark(env: AppArgs) {
    match env {
        AppArgs::Benchmark {
            result,
            answer,
            index,
            format,
        } => {
            if result.is_none() || answer.is_none() || index.is_none() {
                print_log_msg(FAIL, "Result, answer, and index files must be provided");
                std::process::exit(1);
            }
            let result_path = result.unwrap();
            let answer_path = answer.unwrap();
            let index_path = index.unwrap();
            let lookup_path = format!("{}.lookup", index_path);
            let config_path = format!("{}.type", index_path);

            let format = format.as_str();
            // let result = read_one_column_of_tsv(&result_path, 0);
            let result = read_one_column_of_tsv_as_vec(&result_path, 0);
            let answer = read_one_column_of_tsv(&answer_path, 0);
            let lookup = load_lookup_from_file(&lookup_path).0;
            let lookup = lookup.into_iter().collect::<HashSet<_>>();
            let config = read_index_config_from_file(&config_path);

            let result = HashSet::from_iter(result);
            let metric = compare_target_answer_set(&result, &answer, &lookup);
            // let metric = measure_up_to_k_fp(&result, &answer, &lookup, 5.0);
            
            match format {
                "tsv" => {
                    // lookup, result, answer, lookup_len, result_len, answer_len,
                    // hash_type, num_bin_dist, num_bin_angle,
                    // true_pos, true_neg, false_pos, false_neg, precision, recall, accuracy, f1_score, 
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                        index_path,
                        result_path,
                        answer_path,
                        lookup.len(),
                        result.len(),
                        answer.len(),
                        config.hash_type.to_string(),
                        config.num_bin_dist,
                        config.num_bin_angle,
                        metric.true_pos,
                        metric.true_neg,
                        metric.false_pos,
                        metric.false_neg,
                        metric.precision(),
                        metric.recall(),
                        metric.accuracy(),
                        metric.f1_score(),
                    );
                }
                "default" => {
                    println!("Index: {}", lookup_path);
                    println!("Result: {}", result_path);
                    println!("Answer: {}", answer_path);
                    println!("Total ids: {}", lookup.len());
                    println!("Result length: {}", result.len());
                    println!("Answer length: {}", answer.len());
                    println!("Hash type: {}", config.hash_type.to_string());
                    println!("Number of distance bins: {}", config.num_bin_dist);
                    println!("Number of angle bins: {}", config.num_bin_angle);
                    println!("TP: {}", metric.true_pos);
                    println!("TN: {}", metric.true_neg);
                    println!("FP: {}", metric.false_pos);
                    println!("FN: {}", metric.false_neg);
                    // Print float with 4 decimal places
                    println!("Precision: {:.4}", metric.precision());
                    println!("Recall: {:.4}", metric.recall());
                    println!("Accuracy: {:.4}", metric.accuracy());
                    println!("F1 score: {:.4}", metric.f1_score());
                }
                _ => {
                    print_log_msg(FAIL, "Invalid format");
                    std::process::exit(1);
                }
            }
        }
        _ => {
            eprintln!("Invalid subcommand");
            std::process::exit(1);
        }
    }
}

fn read_one_column_of_tsv(file: &str, col_index: usize) -> HashSet<String> {
    let mut set = HashSet::new();
    // Open file and get specific column
    let file = std::fs::File::open(file).expect(
        &log_msg(FAIL, &format!("Failed to open tsv file: {}", file))
    );
    let reader = std::io::BufReader::new(file);
    reader.lines().for_each(|line| {
        let line = line.expect(&log_msg(FAIL, "Failed to read line"));
        let row = line.split('\t').collect::<Vec<_>>();
        set.insert(row[col_index].to_string());
    });
    set
}

fn read_one_column_of_tsv_as_vec(file: &str, col_index: usize) -> Vec<String> {
    let mut vec = Vec::new();
    // Open file and get specific column
    let file = std::fs::File::open(file).expect(
        &log_msg(FAIL, &format!("Failed to open tsv file: {}", file))
    );
    let reader = std::io::BufReader::new(file);
    reader.lines().for_each(|line| {
        let line = line.expect(&log_msg(FAIL, "Failed to read line"));
        let row = line.split('\t').collect::<Vec<_>>();
        vec.push(row[col_index].to_string());
    });
    vec
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore]
    fn test_benchmark() {
        let result = Some("data/zinc_folddisco.tsv".to_string());
        let answer = Some("data/zinc_answer.tsv".to_string());
        let index = Some("analysis/h_sapiens/d16a4/index_id".to_string());
        let format = "tsv";
        let env = AppArgs::Benchmark {
            result,
            answer,
            index,
            format: format.to_string(),
        };
        benchmark(env);
    }
}
