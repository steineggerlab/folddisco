use std::collections::HashSet;

use crate::cli::*;
use crate::index::lookup::load_lookup_from_file;
use rayon::prelude::*;
use crate::prelude::*;

pub fn query_test_for_swissprot(env: AppArgs) {
    match env {
        AppArgs::Test {
            index_path,
            verbose,
        } => {
            let start = std::time::Instant::now();
            let lookup_path = format!("{}.lookup", index_path);
            let (path_vec, _, _) = load_lookup_from_file(&lookup_path);
            let lap1 = std::time::Instant::now();
            if verbose { println!("[INFO ] Time elapsed for loading lookup table {:?}", lap1 - start); }
            let index_table = IndexTable::load_from_bin_custom(&index_path).expect("[ERROR] Failed to load index table");
            let lap2 = std::time::Instant::now();
            if verbose { println!("[INFO ] Time elapsed for loading index table {:?}", lap2 - lap1); }
            // Measure time
            let queries = vec![
                // HashValue::from_u64(6597069965577u64),
                // HashValue::from_u64(4398046577927u64),
                // HashValue::from_u64(8800438323719u64),
                // HashValue::from_u64(5501903505417u64),
                // HashValue::from_u64(5497558337804u64),
                // HashValue::from_u64(7696631925513u64),
                HashValue::from_u64(2203335001092u64),
                HashValue::from_u64(3298551727108u64),
                HashValue::from_u64(3298551661828u64),
            ];
            let queries_u64 = queries.iter().map(|x| x.as_u64()).collect::<Vec<u64>>();
            let result = index_table.query_multiple_with_connectivity(&queries_u64, 2);

            let mut str_result = Vec::new();
            if let Some(result) = result {
                let mut set = HashSet::new();
                let mut dedup_result = Vec::new();
                println!("Total result: {}", result.len());
                for i in result {
                    if !set.contains(&i) {
                        set.insert(i);
                        dedup_result.push(i);
                    }
                }
                // println!("Queried ind (3) {:?}", dedup_result);
                str_result = dedup_result.iter().map(|x| path_vec.get(x.get_id()).unwrap().clone()).collect::<Vec<_>>();
                println!("[INFO ] Queried {:?}", str_result);
            }
            let end = std::time::Instant::now();
            println!("[INFO ] Time elapsed for quering {:?}", end - lap2);

            let answer = vec![
                "/fast/hyunbin/foldcomp/ecoli_pdb/AF-P76176-F1-model_v4.pdb", "/fast/hyunbin/foldcomp/ecoli_pdb/AF-P0C0V0-F1-model_v4.pdb", "/fast/hyunbin/foldcomp/ecoli_pdb/AF-P39099-F1-model_v4.pdb", 
                "/fast/hyunbin/foldcomp/ecoli_pdb/AF-P0AEE3-F1-model_v4.pdb", 
                "/fast/hyunbin/foldcomp/ecoli_pdb/AF-Q7BSW5-F1-model_v4.pdb", "/fast/hyunbin/foldcomp/ecoli_pdb/AF-Q47692-F1-model_v4.pdb", 
                "/fast/hyunbin/foldcomp/ecoli_pdb/AF-O68900-F1-model_v4.pdb", "/fast/hyunbin/foldcomp/ecoli_pdb/AF-Q7BS42-F1-model_v4.pdb", "/fast/hyunbin/foldcomp/ecoli_pdb/AF-Q84GK0-F1-model_v4.pdb", 
                "/fast/hyunbin/foldcomp/ecoli_pdb/AF-Q9EZE7-F1-model_v4.pdb", "/fast/hyunbin/foldcomp/ecoli_pdb/AF-Q8FDW4-F1-model_v4.pdb",
             ];
            let answer  = answer.iter().map(|x| x.to_string()).collect::<Vec<_>>();

            let metric = compare_target_answer( &str_result, &answer, &path_vec);
            // let metric2 = compare_target_answer( &str_result, &answer2, &path_vec);
            let metric_vec = vec![metric];
            for metric in metric_vec {
                println!("Precision: {}", metric.precision());
                println!("Recall: {}", metric.recall()); 
                println!("Accuracy: {}", metric.accuracy());
                println!("F1 score: {}", metric.f1_score());
                println!("TP: {}", metric.true_pos);
                println!("TN: {}", metric.true_neg);
                println!("FP: {}", metric.false_pos);
                println!("FN: {}", metric.false_neg);
            }
        }
        _ => {
            eprintln!("Invalid subcommand");
            std::process::exit(1);
        }
    }
    
}