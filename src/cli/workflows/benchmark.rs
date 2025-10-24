use std::collections::HashSet;
use std::io::BufRead;

use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::cli::*;
use crate::index::lookup::load_lookup_from_file;
use crate::prelude::*;

use crate::cli::config::read_index_config_from_file;
use crate::utils::benchmark::{compare_target_answer_neutral_vec, compare_target_answer_vec, measure_up_to_k_fp_vec, measure_up_to_k_fp_with_neutral_vec};

// usage: folddisco benchmark -r <result.tsv> -a <answer.tsv> -i <index> -f tsv
// usage: folddisco benchmark -r <result.tsv> -a <answer.tsv> -i <index> -f default
// 2025-01-24 17:57:19
// TODO list
// 1. DONE: [ ] consider sep character based on extension (csv -> ',', tsv -> '\t')
// 2. DONE: [ ] id column index for all input files; result, answer, neutral
// 3. TODO: [ ] header option for all input files; result, answer, neutral; Default: false
// 4. TODO: [ ] add documentation.

// EXAMPLE pyScoMotif result
// ,matched_motif,similar_motif_found,RMSD,n_mutations,PDB_ID,header_description
// 0,A49L A143V A77C A103F A136R,A49L A143V A77C A103F A136R,0.002,0,d1iapa_,
// ,matched_motif,similar_motif_found,RMSD,n_mutations,PDB_ID,header_description
// 0,A360L A397N A365G A363E A400I A391R A402G A388G A413K,A360L A397N A365G A363E A400I A391R A402G A388G A413K,0.005,0,d6c3ma3,

pub const HELP_BENCHMARK: &str = "\
usage: folddisco benchmark -r <result.tsv> -a <answer.tsv> -i <index> [options]

input:
    -r, --result <FILE>          Path to result file (tsv or csv)
    -a, --answer <FILE>          Path to answer file (tsv or csv)
    -i, --index <INDEX_PREFIX>   Path to index prefix
    -n, --neutral <FILE>         Path to neutral file (tsv or csv) [optional]
    --input <FILE>               Path to input file containing result and answer file paths [optional]
    
options:
    --fp <int>                   Measure up to k false positives [optional]
    --afdb-to-uniprot            Parse AFDB IDs to UniProt IDs to match answer IDs
    --column-result <int>        Column index for result IDs (default: 0)
    --column-answer <int>        Column index for answer IDs (default: 0)
    --column-neutral <int>       Column index for neutral IDs (default: 0)
    --header-result              Whether result file has header (default: false)
    --header-answer              Whether answer file has header (default: false)
    --header-neutral             Whether neutral file has header (default: false)
";


pub fn benchmark(env: AppArgs) {
    match env {
        AppArgs::Benchmark {
            result,
            answer,
            neutral,
            index,
            input,
            format,
            fp,
            threads,
            afdb_to_uniprot,
            column_result,
            column_answer,
            column_neutral,
            header_result,
            header_answer,
            header_neutral,
        } => {
            if input.is_none() {
                if result.is_none() || answer.is_none() || index.is_none() {
                    print_log_msg(FAIL, "Result, answer, and index files must be provided");
                    eprintln!("{}", HELP_BENCHMARK);
                    std::process::exit(1);
                }
            }
            // If input is given, read from file
            let input_vector = if input.is_none() {
                if neutral.is_some() {
                    vec![(result.unwrap(), answer.unwrap(), neutral)]
                } else {
                    vec![(result.unwrap(), answer.unwrap(), None)]
                }
            } else {
                let input = input.unwrap();
                let file = std::fs::File::open(&input).expect(
                    &log_msg(FAIL, &format!("Failed to open input file: {}", input))
                );
                let reader = std::io::BufReader::new(file);
                reader.lines().map(|line| {
                    let line = line.expect(&log_msg(FAIL, "Failed to read line"));
                    let row = line.split('\t').collect::<Vec<_>>();
                    if row.len() == 2 {
                        (row[0].to_string(), row[1].to_string(), None)
                    } else if row.len() >= 3 {
                        (row[0].to_string(), row[1].to_string(), Some(row[2].to_string()))
                    } else {
                        print_log_msg(FAIL, "Invalid input format");
                        std::process::exit(1);
                    }
                    // (row[0].to_string(), row[1].to_string())
                }).collect::<Vec<_>>()
            };
            
            // TODO: Add false list (3rd column), neutral list (4th column)
            
            let _pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            
            let index_path = index.unwrap();
            let lookup_path = format!("{}.lookup", index_path);
            let config_path = format!("{}.type", index_path);
            let format = format.as_str();
            let raw_lookup = load_lookup_from_file(&lookup_path);
            let mut lookup = raw_lookup.into_iter().map(|(id, _, _, _, _)| parse_path(&id, afdb_to_uniprot).to_string()).collect::<Vec<_>>();
            lookup.sort();
            lookup.dedup();
            let config = read_index_config_from_file(&config_path);

            input_vector.par_iter().for_each(|(result_path, answer_path, neutral)| {
                // Parse path by id type
                let result = read_one_column_as_vec(&result_path, column_result, header_result, afdb_to_uniprot);

                let answer = read_one_column_as_vec(&answer_path, column_answer, header_answer, afdb_to_uniprot);

                let metric = if let Some(neutral) = neutral {
                    let neutral = read_one_column_as_vec(&neutral, column_neutral, header_neutral, afdb_to_uniprot);

                    if let Some(fp) = fp {
                        measure_up_to_k_fp_with_neutral_vec(&result, &answer, &neutral, &lookup, fp)
                    } else {
                        compare_target_answer_neutral_vec(&result, &answer, &neutral, &lookup)
                    }
                } else {
                    if let Some(fp) = fp {
                        measure_up_to_k_fp_vec(&result, &answer, &lookup, fp)
                    } else {
                        compare_target_answer_vec(&result, &answer, &lookup)
                    }
                };

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
            });

        }
        _ => {
            eprintln!("{}", HELP_BENCHMARK);
            std::process::exit(1);
        }
    }
}

fn read_one_column_as_vec(file: &str, col_index: usize, header: bool, afdb_to_uniprot: bool) -> Vec<String> {
    let mut vec = Vec::new();
    let mut seen = HashSet::new();
    let sep = get_sep_character_from_filename(file);
    // Open file and get specific column
    let file = std::fs::File::open(file).expect(
        &log_msg(FAIL, &format!("Failed to open tsv file: {}", file))
    );
    let reader = std::io::BufReader::new(file);
    
    // If header is true, skip first line
    let lines = reader.lines().skip(if header { 1 } else { 0 });
    for line in lines {
        let line = line.expect(&log_msg(FAIL, "Failed to read line"));
        let row = line.split(sep).collect::<Vec<_>>();
        let value = row[col_index];
        let parsed_value = parse_path(value, afdb_to_uniprot);
        if !seen.contains(parsed_value) {
            vec.push(parsed_value.to_string());
            seen.insert(parsed_value.to_string());
        }
    }

    vec
}

#[inline]
fn get_sep_character_from_filename(file: &str) -> char {
    let ext = file.split('.').last().unwrap();
    match ext {
        "csv" => ',',
        "tsv" | "txt" => '\t',
        _ => {
            print_log_msg(WARN, "Unexpected file extension. Using tab as separator");
            '\t'
        }
    }
}

#[inline]
fn parse_path(path: &str, afdb_to_uniprot: bool) -> &str {
    let path = path.split('/').last().unwrap();
    let path = if path.ends_with(".pdb") || path.ends_with(".cif") || path.ends_with(".fcz") || path.ends_with(".ent") {
        // Return slice of string from start to end-4
        &path[..path.len()-4]
    } else if path.ends_with(".pdb.gz") || path.ends_with(".cif.gz") || path.ends_with(".fcz.gz") || path.ends_with(".ent.gz") {
        // Return slice of string from start to end-7
        &path[..path.len()-7]
    } else {
        path
    };
    if afdb_to_uniprot {
        let split = path.split("-").collect::<Vec<_>>();
        if split.len() >= 2 {
            split[1]
        } else {
            path
        }
    } else {
        path
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore]
    fn test_benchmark() {
        let env = AppArgs::Benchmark {
            result: Some("data/zinc_folddisco.tsv".to_string()),
            answer: Some("data/zinc_answer.tsv".to_string()),
            neutral: None,
            index: Some("index/h_sapiens_folddisco".to_string()),
            input: None,
            format: "tsv".to_string(),
            fp: None,
            threads: 1,
            afdb_to_uniprot: false,
            column_result: 0,
            column_answer: 0,
            column_neutral: 0,
            header_result: false,
            header_answer: false,
            header_neutral: false,
        };
        benchmark(env);
    }
}
