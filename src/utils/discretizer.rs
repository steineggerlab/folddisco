/*
 * File: discretizer.rs
 * Created: 2023-01-18 14:31:12
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Copyright © 2023 Hyunbin Kim, All rights reserved
 */

#![allow(dead_code, unused_variables)]

use std::result::Result;
// use rand::Rng;
// use test::Bencher; // unstable

/* Structs for Discretizer */

#[derive(Debug)]
/// DiscParams is a struct that contains the parameters for discretization.
pub struct DiscParams {
    /// Minimum value of the range to be discretized.
    pub min: f32,
    /// Maximum value of the range to be discretized.
    pub max: f32,
    /// Number of bins to be discretized.
    pub n_bins: u32,
    /// Discretization factor.
    pub disc_f: f32,
    /// Continuization factor.
    pub cont_f: f32,
}

#[derive(Debug)]
/// Discretizer is a struct that contains the parameters for discretization and
/// the raw f32 vector to be discretized, and the discretized vector.
/// Discretizer should be created as a mutable instance to call the discretize() method.
pub struct Discretizer {
    /// Discretization parameters.
    pub disc_params: DiscParams,
    /// Raw f32 vector to be discretized.
    pub raw_vec: Vec<f32>,
    /// Discretized vector.
    pub disc_vec: Vec<Option<u32>>,
}

/* Methods for Discretizer */

impl DiscParams {
    pub fn new(min: f32, max: f32, n_bins: u32, disc_f: f32, cont_f: f32) -> DiscParams {
        DiscParams {
            min,
            max,
            n_bins,
            disc_f,
            cont_f,
        }
    }
}

impl Discretizer {
    /// Discretizer::new() returns a new Discretizer instance.
    /// new() consumes raw_vec.
    pub fn new(raw_vec: Vec<f32>, n_bins: u32) -> Result<Discretizer, &'static str> {
        let disc_param_result = get_disc_params_from_vec(&raw_vec, n_bins);
        match disc_param_result {
            Ok(disc_params) => return Ok(Discretizer {
                disc_params,
                raw_vec,
                disc_vec: Vec::new(),
            }),
            Err(e) => return Err(e),
        }
    }
    /// Discretizer::discretize() discretizes the raw_vec and stores the result in disc_vec.
    pub fn discretize(&mut self) {
        self.disc_vec = discretize_f32_vec(&self.raw_vec, &self.disc_params);
    }
    /// Discretizer::continuize() continuizes the disc_vec and returns the result.
    pub fn continuize(&self) -> Vec<f32> {
        continuize_u32_vec(&self.disc_vec, &self.disc_params)
    }
}


/* Functions for Discretizer */

pub fn get_disc_params_from_vec(vec: &Vec<f32>, n_bins: u32) -> Result<DiscParams, &'static str> {
    // Error handling.
    // case 1: raw_vec.len() == 0
    if vec.len() == 0 {
        return Err("raw_vec cannot be empty.");
    }
    // case 2: n_bins == 0
    if n_bins == 0 || n_bins == 1 {
        return Err("n_bins cannot be zero or one.");
    }
    // iter().fold() is a functional programming style.
    // .fold() is a method that accumulates the result of applying a function to each element of an iterator.
    // https://doc.rust-lang.org/std/iter/trait.Iterator.html#method.fold
    let min = vec.iter().fold(f32::INFINITY, |acc, &x| acc.min(x));
    let max = vec.iter().fold(f32::NEG_INFINITY, |acc, &x| acc.max(x));
    let cont_f = (max - min) / ((n_bins - 1) as f32);
    let disc_f = 1.0 / cont_f;
    // Check infinities.
    if min.is_infinite() || max.is_infinite() {
        return Err("min or max is infinite.");
    }
    Ok(DiscParams::new(min, max, n_bins, disc_f, cont_f))
}

pub fn discretize_f32(val: f32, disc_params: &DiscParams) -> Result<u32, &'static str> {
    // Check if val is in the range of disc_params.
    if val < disc_params.min || val > disc_params.max {
        return Err("val is out of the range of disc_params.");
    }
    // Check if disc_params.disc_f is not zero.
    if disc_params.disc_f == 0.0 {
        return Err("disc_params.disc_f cannot be zero.");
    }
    Ok((((val - disc_params.min) * disc_params.disc_f) + 0.5) as u32)
}

pub fn discretize_f32_vec(vec: &Vec<f32>, disc_params: &DiscParams) -> Vec<Option<u32>> {
    // Avoid reallocation by using Vec::with_capacity().
    let mut disc_vec = Vec::with_capacity(vec.len());
    for i in 0..vec.len() {
        match discretize_f32(vec[i], disc_params) {
            Ok(val) => disc_vec.push(Some(val)),
            Err(e) => disc_vec.push(None),
        }
    }
    disc_vec
}

pub fn continuize_u32(val: u32, disc_params: &DiscParams) -> f32 {
    // Check if disc_params.cont_f is not zero.
    if disc_params.cont_f == 0.0 {
        panic!("disc_params.cont_f cannot be zero.");
    }
    (val as f32) * disc_params.cont_f + disc_params.min
}

pub fn continuize_u32_vec(vec: &Vec<Option<u32>>, disc_params: &DiscParams) -> Vec<f32> {
    let mut cont_vec = Vec::new();
    for i in 0..vec.len() {
        match vec[i] {
            Some(val) => cont_vec.push(continuize_u32(val, disc_params)),
            None => cont_vec.push(f32::NAN),
        }
    }
    cont_vec
}

/* TESTS */
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn discretizer_get_disc_params_from_vec() {
        let vec = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        let disc_params = get_disc_params_from_vec(&vec, 5).unwrap();
        assert_eq!(disc_params.min, 2.0);
        assert_eq!(disc_params.max, 10.0);
        assert_eq!(disc_params.n_bins, 5);
        assert_eq!(disc_params.disc_f, 0.5);
        assert_eq!(disc_params.cont_f, 2.0);
    }
    #[test]
    fn discretizer_discretize_f32() {
        let vec = vec![36.7, 34.1, 28.9, 16.8, 24.8, 29.0, 38.2];
        let mut disc = Discretizer::new(vec, 4).unwrap();
        disc.discretize();
        assert_eq!(disc.disc_vec[0], Some(3));
        assert_eq!(disc.disc_vec[1], Some(2));
        assert_eq!(disc.disc_vec[2], Some(2));
        assert_eq!(disc.disc_vec[3], Some(0));
        assert_eq!(disc.disc_vec[4], Some(1));
        assert_eq!(disc.disc_vec[5], Some(2));
        assert_eq!(disc.disc_vec[6], Some(3));
    }
    #[test]
    fn discretizer_continuize_u32() {
        let vec = vec![36.7, 34.1, 28.9, 16.8, 24.8, 29.0, 38.2];
        let mut disc = Discretizer::new(vec, 4).unwrap();
        disc.discretize();
        let cont_vec = disc.continuize();
        assert_eq!(cont_vec[0], 38.2);
        assert_eq!(cont_vec[1], 31.066666);
        assert_eq!(cont_vec[2], 31.066666);
        assert_eq!(cont_vec[3], 16.8);
        assert_eq!(cont_vec[4], 23.933332);
        assert_eq!(cont_vec[5], 31.066666);
        assert_eq!(cont_vec[6], 38.2);
    }
    #[test]
    fn discretizer_empty_vec() {
        let vec = vec![];
        let disc = Discretizer::new(vec, 4);
        assert_eq!(disc.is_err(), true);
    }
    #[test]
    fn discretizer_one_zero_bins() {
        let vec = vec![36.7, 34.1, 28.9, 16.8, 24.8, 29.0, 38.2];
        let vec2 = vec.clone();
        let disc = Discretizer::new(vec, 0);
        assert_eq!(disc.is_err(), true);
        // check error message too
        let err = disc.unwrap_err();
        assert_eq!(err, "n_bins cannot be zero or one.");
        let disc = Discretizer::new(vec2, 1);
        assert_eq!(disc.is_err(), true);
        // check error message too
        let err = disc.unwrap_err();
        assert_eq!(err, "n_bins cannot be zero or one.");
    }
    #[test]
    fn discretizer_check_infinite() {
        let vec = vec![36.7, 34.1, 28.9, 16.8, 24.8, 29.0, 38.2, f32::INFINITY];
        let disc = Discretizer::new(vec, 4);
        assert_eq!(disc.is_err(), true);
        // check error message too
        let err = disc.unwrap_err();
        assert_eq!(err, "min or max is infinite.");
    }

    // #[bench]
    // fn bench_discretizer(b: &mut Bencher) {
    //     // Randomly generated vector of 1000 elements.
    //     let vec = Rng.gen_iter::<f32>().take(1000).collect();
    //     b.iter(|| {
    //         let mut disc = Discretizer::new(vec.clone(), 4).unwrap();
    //         disc.discretize();
    //     });
    // }

}

