// File: benchmark.rs
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

// Requirements
// 1. Folder of PDB files
// 2. TSV for wanted PDB files & active sites

pub struct Metrics {
    pub true_pos: f64,
    pub true_neg: f64,
    pub false_pos: f64,
    pub false_neg: f64,
}

impl Metrics {
    pub fn new(true_pos: f64, true_neg: f64, false_pos: f64, false_neg: f64) -> Self {
        Self { 
            true_pos, 
            true_neg, 
            false_pos, 
            false_neg 
        }
    }

    pub fn precision(&self) -> f64 {
        self.true_pos / (self.true_pos + self.false_pos)
    }

    pub fn recall(&self) -> f64 {
        self.true_pos / (self.true_pos + self.false_neg)
    }

    pub fn accuracy(&self) -> f64 {
        (self.true_pos + self.true_neg) / (self.true_pos + self.true_neg + self.false_pos + self.false_neg)
    }

    pub fn f1_score(&self) -> f64 {
        let precision = self.precision();
        let recall = self.recall();
        2.0 * (precision * recall) / (precision + recall)
    }
}

pub fn calculate_metrics(query: Vec<usize>, answer: Vec<usize>) -> Metrics {
    assert_eq!(query.len(), answer.len());

    let mut true_pos = 0.0;
    let mut true_neg = 0.0;
    let mut false_pos = 0.0;
    let mut false_neg = 0.0;

    for (q, a) in query.iter().zip(answer.iter()) {
        match (*q, *a) {
            (1, 1) => true_pos += 1.0,
            (0, 0) => true_neg += 1.0,
            (1, 0) => false_pos += 1.0,
            (0, 1) => false_neg += 1.0,
            _ => panic!("Invalid values in query or answer"),
        }
    }

    Metrics::new(true_pos, true_neg, false_pos, false_neg)
}

pub fn compare_target_answer<T: Eq + PartialEq>(target: &Vec<T>, answer: &Vec<T>, all: &Vec<T>) -> Metrics {
    // True Positive: Elements in both target and answer
    let true_pos = target.iter().filter(|&x| answer.contains(x)).count() as f64;
    // True Negative: Elements in neither target nor answer
    let true_neg = all.iter().filter(|&x| !target.contains(x) && !answer.contains(x)).count() as f64;
    // False Positive: Elements in target but not in answer
    let false_pos = target.iter().filter(|&x| !answer.contains(x)).count() as f64;
    // False Negative: Elements in answer but not in target
    let false_neg = answer.iter().filter(|&x| !target.contains(x)).count() as f64;

    Metrics::new(true_pos, true_neg, false_pos, false_neg)
}
