// File: benchmark.rs
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2023 Hyunbin Kim, All rights reserved

// Confusion-matrix metrics for comparing search results against answer sets.

/// Confusion-matrix counts.
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

/// Metrics of `target` (results) against `answer` (true set) within `all`.
pub fn compare_target_answer_vec<T: Eq + PartialEq>(target: &Vec<T>, answer: &Vec<T>, all: &Vec<T>) -> Metrics {
    let true_pos = target.iter().filter(|&x| answer.contains(x)).count() as f64;
    let true_neg = all.iter().filter(|&x| !target.contains(x) && !answer.contains(x)).count() as f64;
    let false_pos = target.iter().filter(|&x| !answer.contains(x)).count() as f64;
    let false_neg = answer.iter().filter(|&x| !target.contains(x)).count() as f64;

    Metrics::new(true_pos, true_neg, false_pos, false_neg)
}

/// Like `compare_target_answer_vec`, but `neutral` items are neither FP nor TN.
pub fn compare_target_answer_neutral_vec<T: Eq + PartialEq>(target: &Vec<T>, answer: &Vec<T>, neutral: &Vec<T>, all: &Vec<T>) -> Metrics {
    let true_pos = target.iter().filter(|&x| answer.contains(x)).count() as f64;
    let true_neg = all.iter().filter(|&x| !target.contains(x) && !answer.contains(x) && !neutral.contains(x)).count() as f64;
    let false_pos = target.iter().filter(|&x| !answer.contains(x) && !neutral.contains(x)).count() as f64;
    let false_neg = answer.iter().filter(|&x| !target.contains(x)).count() as f64;

    Metrics::new(true_pos, true_neg, false_pos, false_neg)
}

/// Metrics over the ranked `target` prefix that ends at the `k`-th false positive.
pub fn measure_up_to_k_fp_vec<T: Eq + PartialEq>(target: &Vec<T>, answer: &Vec<T>, all: &Vec<T>, k: f64) -> Metrics {
    let mut true_pos = 0.0;
    let mut false_pos = 0.0;
    for t in target {
        if answer.contains(t) {
            true_pos += 1.0;
        } else {
            false_pos += 1.0;
        }
        if false_pos >= k {
            break;
        }
    }
    let false_neg = answer.len() as f64 - true_pos;
    let true_neg = all.len() as f64 - (true_pos + false_pos + false_neg);
    Metrics::new(true_pos, true_neg, false_pos, false_neg)
}

/// Like `measure_up_to_k_fp_vec`, skipping `neutral` items.
pub fn measure_up_to_k_fp_with_neutral_vec<T: Eq + PartialEq>(target: &Vec<T>, answer: &Vec<T>, neutral: &Vec<T>, all: &Vec<T>, k: f64) -> Metrics {
    let mut true_pos = 0.0;
    let mut false_pos = 0.0;
    for t in target {
        if answer.contains(t) {
            true_pos += 1.0;
        } else if !neutral.contains(t) {
            false_pos += 1.0;
        }
        if false_pos >= k {
            break;
        }
    }
    let false_neg = answer.len() as f64 - true_pos;
    let true_neg = all.len() as f64 - (true_pos + false_pos + false_neg);
    Metrics::new(true_pos, true_neg, false_pos, false_neg)
}
