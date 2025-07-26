// File: kabsch.rs
// Created: 2025-07-11 10:16:18
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2025 Hyunbin Kim, All rights reserved
// Kabsch algorithm implementation for optimal superposition of two coordinate sets
// Original code from TM-align ()
// Translated to Rust by Martin Steinegger & Hyunbin Kim


use crate::structure::coordinate::Coordinate;

// Higher precision constants for f64 version
const EPSILON: f64 = 1.0e-8;
const TOLERANCE: f64 = 0.01;
const SQRT3: f64 = 1.7320508075688772; // sqrt(3)

const IP: [usize; 9] = [0, 1, 3, 1, 2, 4, 3, 4, 5]; // Index permutation for 3x3 matrix
const IP2312: [usize; 4] = [1, 2, 0, 1]; // Index permutation for 2x2 matrix


#[derive(Debug)]
pub struct KabschSuperimposer {
    pub reference_coords: Option<Vec<[f32; 3]>>,
    pub coords: Option<Vec<[f32; 3]>>,
    pub transformed_coords: Option<Vec<[f32; 3]>>,
    pub rot: Option<[[f32; 3]; 3]>,
    pub tran: Option<[f32; 3]>,
    pub rms: Option<f32>,
    pub init_rms: Option<f32>,
    pub tm_score: Option<f32>,
    pub natoms: usize,
}

impl KabschSuperimposer {
    pub fn new() -> Self {
        Self {
            reference_coords: None,
            coords: None,
            transformed_coords: None,
            rot: None,
            tran: None,
            rms: None,
            init_rms: None,
            tm_score: None,
            natoms: 0,
        }
    }

    /// Calculate TM-score between two aligned structures
    /// TM-score is normalized by the target structure length
    pub fn get_tm_score(&self) -> f32 {
        self.tm_score.expect("TM-score not calculated. Run superposition first.")
    }

    /// Calculate TM-score normalized by a specific length (useful for comparing structures of different lengths)
    pub fn get_tm_score_normalized(&self, norm_length: usize) -> f32 {
        let coords = self.coords.as_ref().expect("Coordinates not set");
        let reference_coords = self.reference_coords.as_ref().expect("Reference coordinates not set");
        let rot = self.rot.expect("Superposition not performed");
        let tran = self.tran.expect("Superposition not performed");
        
        calculate_tm_score_with_length(coords, reference_coords, rot, tran, norm_length)
    }

    pub fn set_atoms(&mut self, fixed: &[Coordinate], moving: &[Coordinate]) {
        assert!(fixed.len() == moving.len(), "Fixed and moving atom lists differ in size");

        let fixed_coords: Vec<[f32; 3]> = fixed.iter().map(|c| c.to_array()).collect();
        let moving_coords: Vec<[f32; 3]> = moving.iter().map(|c| c.to_array()).collect();

        self.set(&fixed_coords, &moving_coords);
        self.run();
        self.rms = Some(self.get_rms());
        
        // Calculate and cache TM-score
        let coords = self.coords.as_ref().unwrap();
        let reference_coords = self.reference_coords.as_ref().unwrap();
        let rot = self.rot.unwrap();
        let tran = self.tran.unwrap();
        self.tm_score = Some(calculate_tm_score(coords, reference_coords, rot, tran));
    }

    pub fn set(&mut self, reference_coords: &[[f32; 3]], coords: &[[f32; 3]]) {
        self.reference_coords = Some(reference_coords.to_vec());
        self.coords = Some(coords.to_vec());
        self.natoms = coords.len();

        assert!(
            reference_coords.len() == coords.len(),
            "Coordinates must have the same dimensions."
        );
        assert!(
            coords[0].len() == 3,
            "Coordinates must be Nx3 arrays."
        );
    }

    pub fn run(&mut self) {
        let coords = self.coords.clone().unwrap();
        let reference_coords = self.reference_coords.clone().unwrap();
        
        let (rot, tran, rmsd) = match kabsch(&coords, &reference_coords, 2) {
            Some(result) => result,
            None => ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [0.0, 0.0, 0.0], f32::MAX),
        };
        
        self.rot = Some(rot);
        self.tran = Some(tran);
        self.rms = Some(rmsd);
        
        // Calculate TM-score
        self.tm_score = Some(calculate_tm_score(&coords, &reference_coords, rot, tran));
        
        // Get transformed coordinates
        self.transformed_coords = Some(
            coords.iter()
                .map(|&coord| {
                    let rotated = matrix_vector_multiply(rot, coord);
                    add_vec(rotated, tran)
                })
                .collect()
        );
    }

    pub fn get_transformed(&mut self) -> Vec<[f32; 3]> {
        if self.transformed_coords.is_none() {
            let coords = self.coords.clone().unwrap();
            let rot = self.rot.unwrap();
            let tran = self.tran.unwrap();

            self.transformed_coords = Some(
                coords.iter()
                    .map(|&coord| {
                        let rotated = matrix_vector_multiply(rot, coord);
                        
                        add_vec(rotated, tran)
                    })
                    .collect()
            );
        }
        self.transformed_coords.clone().unwrap()
    }

    pub fn get_rotran(&self) -> ([[f32; 3]; 3], [f32; 3]) {
        (self.rot.unwrap(), self.tran.unwrap())
    }

    pub fn get_init_rms(&mut self) -> f32 {
        if self.init_rms.is_none() {
            let coords = self.coords.clone().unwrap();
            let reference_coords = self.reference_coords.clone().unwrap();
            // without superposition, calculate initial RMSD
            self.init_rms = Some(
                (0..coords.len())
                    .map(|i| {
                        let dx = coords[i][0] - reference_coords[i][0];
                        let dy = coords[i][1] - reference_coords[i][1];
                        let dz = coords[i][2] - reference_coords[i][2];
                        dx * dx + dy * dy + dz * dz
                    })
                    .sum::<f32>()
                    .sqrt() / (coords.len() as f32).sqrt()
            );
        }
        self.init_rms.unwrap()
    }

    pub fn get_rms(&self) -> f32 {
        self.rms.unwrap()
    }
}

fn matrix_vector_multiply(matrix: [[f32; 3]; 3], vector: [f32; 3]) -> [f32; 3] {
    [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    ]
}
fn add_vec(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// High-precision version of the Kabsch algorithm using f64 for better numerical stability
fn kabsch(x: &[[f32; 3]], y: &[[f32; 3]], mode: u8) -> Option<([[f32; 3]; 3], [f32; 3], f32)> {
    // Convert input f32 arrays to f64 for computation
    let x_f64: Vec<[f64; 3]> = x.iter().map(|arr| [arr[0] as f64, arr[1] as f64, arr[2] as f64]).collect();
    let y_f64: Vec<[f64; 3]> = y.iter().map(|arr| [arr[0] as f64, arr[1] as f64, arr[2] as f64]).collect();
    
    // -----------------------------------------------------------------------
    // 1. Basic sanity checks -------------------------------------------------
    // -----------------------------------------------------------------------
    let n = x_f64.len();
    if n == 0 || y_f64.len() != n {
        return Some(([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [0.0, 0.0, 0.0], f32::MAX));
    }

    // -----------------------------------------------------------------------
    // 3. Scratch scalars / arrays -------------------------------------------
    // -----------------------------------------------------------------------
    let mut rms = 0.0_f64;
    let mut e0 = 0.0_f64;
    let mut _rms1 = 0.0_f64;
    let sigma: f64;

    // vectors (f64)
    let mut s1: [f64; 3] = [0.0; 3];
    let mut s2: [f64; 3] = [0.0; 3];
    let mut sx: [f64; 3] = [0.0; 3];
    let mut sy: [f64; 3] = [0.0; 3];
    let mut sz: [f64; 3] = [0.0; 3];
    let mut xc: [f64; 3] = [0.0; 3];
    let mut yc: [f64; 3] = [0.0; 3];
    let mut t: [f64; 3] = [0.0; 3];
    let mut e: [f64; 3] = [0.0; 3];

    // 3×3 matrices (f64)
    let mut r: [[f64; 3]; 3] = [[0.0; 3]; 3];
    let mut a: [[f64; 3]; 3] = [[0.0; 3]; 3];
    let mut b: [[f64; 3]; 3] = [[0.0; 3]; 3];
    let mut u: [[f64; 3]; 3] = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ];

    // symmetric 3×3 stored as 6 upper‑triangular elements (f64)
    let mut rr: [f64; 6] = [0.0; 6];
    let mut ss: [f64; 6] = [0.0; 6];

    // -----------------------------------------------------------------------
    // 4. Accumulate first‑ and second‑order sums ----------------------------
    // -----------------------------------------------------------------------
    for i in 0..n {
        let c1 = [x_f64[i][0], x_f64[i][1], x_f64[i][2]];
        let c2 = [y_f64[i][0], y_f64[i][1], y_f64[i][2]];

        for j in 0..3 {
            s1[j] += c1[j];
            s2[j] += c2[j];
        }

        sx[0] += c1[0] * c2[0];
        sx[1] += c1[0] * c2[1];
        sx[2] += c1[0] * c2[2];

        sy[0] += c1[1] * c2[0];
        sy[1] += c1[1] * c2[1];
        sy[2] += c1[1] * c2[2];

        sz[0] += c1[2] * c2[0];
        sz[1] += c1[2] * c2[1];
        sz[2] += c1[2] * c2[2];
    }

    // centres of mass
    for j in 0..3 {
        xc[j] = s1[j] / (n as f64);
        yc[j] = s2[j] / (n as f64);
    }

    // pre‑compute e0 if the caller asked for an RMSD
    if mode == 0 || mode == 2 {
        for i in 0..n {
            e0 += (x_f64[i][0] - xc[0]).powi(2) + (y_f64[i][0] - yc[0]).powi(2);
            e0 += (x_f64[i][1] - xc[1]).powi(2) + (y_f64[i][1] - yc[1]).powi(2);
            e0 += (x_f64[i][2] - xc[2]).powi(2) + (y_f64[i][2] - yc[2]).powi(2);
        }
    }

    // -----------------------------------------------------------------------
    // 5. Build cross‑covariance matrix r ------------------------------------
    // -----------------------------------------------------------------------
    for j in 0..3 {
        r[j][0] = sx[j] - s1[0] * s2[j] / (n as f64);
        r[j][1] = sy[j] - s1[1] * s2[j] / (n as f64);
        r[j][2] = sz[j] - s1[2] * s2[j] / (n as f64);
    }

    // determinant of r
    let det_r = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
        - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
        + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    sigma = det_r;

    // -----------------------------------------------------------------------
    // 6. Build transpose(r) * r (upper triangular in rr[0..6]) --------------
    // -----------------------------------------------------------------------
    let mut m = 0;
    for j in 0..3 {
        for i in 0..=j {
            rr[m] = r[0][i] * r[0][j] + r[1][i] * r[1][j] + r[2][i] * r[2][j];
            m += 1;
        }
    }

    // -----------------------------------------------------------------------
    // 7. Characteristic equation of transpose(r) * r -------------------------
    // -----------------------------------------------------------------------
    let spur = (rr[0] + rr[2] + rr[5]) / 3.0;
    let cof = (
        ((rr[2] * rr[5] - rr[4] * rr[4])
            + rr[0] * rr[5] - rr[3] * rr[3]
            + rr[0] * rr[2] - rr[1] * rr[1])
    ) / 3.0;
    let det = det_r * det_r;

    e.fill(spur);

    if spur > 0.0 {
        // -------------------------------------------------------------------
        // 7a. Eigen‑decomposition branch ----
        // -------------------------------------------------------------------
        let d = spur * spur;
        let h = d - cof;
        let g = (spur * cof - det) / 2.0 - spur * h;

        if h > 0.0 {
            let sqrth = h.sqrt();
            let mut disc = h * h * h - g * g;
            if disc < 0.0 {
                disc = 0.0;
            }
            let sqrt_disc = disc.sqrt();
            
            // Robust angle computation to avoid NaN from extreme values
            let d_ang = if g.abs() > 1e18 {  // Higher threshold for f64
                if g > 0.0 {
                    std::f64::consts::PI / 3.0
                } else {
                    0.0
                }
            } else {
                sqrt_disc.atan2(-g) / 3.0
            };
            
            let cth = sqrth * d_ang.cos();
            let sth = sqrth * SQRT3 * d_ang.sin();

            e[0] = spur + 2.0 * cth;
            e[1] = spur - cth + sth;
            e[2] = spur - cth - sth;

            if mode != 0 {
                let mut a_failed = false;
                let mut b_failed = false;
                
                // ------- assemble eigen‑vectors in A ------------------------
                for &l in &[0usize, 2usize] {
                    let d_local = e[l];
                    ss[0] = (d_local - rr[2]) * (d_local - rr[5]) - rr[4] * rr[4];
                    ss[1] = (d_local - rr[5]) * rr[1] + rr[3] * rr[4];
                    ss[2] = (d_local - rr[0]) * (d_local - rr[5]) - rr[3] * rr[3];
                    ss[3] = (d_local - rr[2]) * rr[3] + rr[1] * rr[4];
                    ss[4] = (d_local - rr[0]) * rr[4] + rr[1] * rr[3];
                    ss[5] = (d_local - rr[0]) * (d_local - rr[2]) - rr[1] * rr[1];

                    for s in &mut ss {
                        if s.abs() <= EPSILON {
                            *s = 0.0;
                        }
                    }

                    let j = match (ss[0].abs(), ss[2].abs(), ss[5].abs()) {
                        (a, b, c) if a >= b && a >= c => 0,
                        (_a, b, c) if b >= c => 1,
                        _ => 2,
                    };

                    let mut d_norm = 0.0;
                    for i in 0..3 {
                        let k = IP[3 * j + i];
                        a[i][l] = ss[k];
                        d_norm += ss[k] * ss[k];
                    }

                    if d_norm > EPSILON {
                        d_norm = 1.0 / d_norm.sqrt();
                    } else {
                        d_norm = 0.0;
                    }
                    for i in 0..3 {
                        a[i][l] *= d_norm;
                    }
                }

                // ---------- orthogonalise the middle column -----------------
                let dot = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];

                let (m1, m) = if e[0] - e[1] > e[1] - e[2] {
                    (2usize, 0usize)
                } else {
                    (0usize, 2usize)
                };

                let mut p = 0.0;
                for i in 0..3 {
                    a[i][m1] = a[i][m1] - dot * a[i][m];
                    p += a[i][m1] * a[i][m1];
                }

                if p <= TOLERANCE {
                    // Fallback identical to the C++ block
                    let mut j = 0usize;
                    p = 1.0;
                    for i in 0..3 {
                        if p < a[i][m].abs() {
                            p = a[i][m].abs();
                            j = i;
                        }
                    }
                    let k = IP2312[j];
                    let l = IP2312[j + 1];
                    p = (a[k][m] * a[k][m] + a[l][m] * a[l][m]).sqrt();
                    if p > TOLERANCE {
                        a[j][m1] = 0.0;
                        a[k][m1] = -a[l][m] / p;
                        a[l][m1] = a[k][m] / p;
                    } else {
                        a_failed = true;
                    }
                } else {
                    p = 1.0 / p.sqrt();
                    for i in 0..3 {
                        a[i][m1] *= p;
                    }
                }

                if !a_failed {
                    a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
                    a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
                    a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
                }

                // ------------- Build B --------------------------------------
                if !a_failed {
                    for l in 0..2 {
                        let mut d_b = 0.0;
                        for i in 0..3 {
                            b[i][l] = r[i][0] * a[0][l]
                                + r[i][1] * a[1][l]
                                + r[i][2] * a[2][l];
                            d_b += b[i][l] * b[i][l];
                        }
                        if d_b > EPSILON {
                            d_b = 1.0 / d_b.sqrt();
                        } else {
                            d_b = 0.0;
                        }
                        for i in 0..3 {
                            b[i][l] *= d_b;
                        }
                    }

                    let mut dot_b = 0.0;
                    for i in 0..3 {
                        dot_b += b[i][0] * b[i][1];
                    }

                    let mut p_b = 0.0;
                    for i in 0..3 {
                        b[i][1] -= dot_b * b[i][0];
                        p_b += b[i][1] * b[i][1];
                    }

                    if p_b <= TOLERANCE {
                        p_b = 1.0;
                        let mut j = 0usize;
                        for i in 0..3 {
                            if p_b < b[i][0].abs() {
                                p_b = b[i][0].abs();
                                j = i;
                            }
                        }
                        let k = IP2312[j];
                        let l = IP2312[j + 1];
                        p_b = (b[k][0] * b[k][0] + b[l][0] * b[l][0]).sqrt();
                        if p_b > TOLERANCE {
                            b[j][1] = 0.0;
                            b[k][1] = -b[l][0] / p_b;
                            b[l][1] = b[k][0] / p_b;
                        } else {
                            b_failed = true;
                        }
                    } else {
                        p_b = 1.0 / p_b.sqrt();
                        for i in 0..3 {
                            b[i][1] *= p_b;
                        }
                    }

                    if !b_failed {
                        b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
                        b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
                        b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];

                        // U = B · Aᵀ
                        for i in 0..3 {
                            for j in 0..3 {
                                u[i][j] = b[i][0] * a[j][0]
                                    + b[i][1] * a[j][1]
                                    + b[i][2] * a[j][2];
                            }
                        }

                        // translation
                        for i in 0..3 {
                            t[i] = yc[i] - (u[i][0] * xc[0] + u[i][1] * xc[1] + u[i][2] * xc[2]);
                        }
                    }
                }
            }
        }
    } else {
        // purely translational fit
        for i in 0..3 {
            t[i] = yc[i] - (u[i][0] * xc[0] + u[i][1] * xc[1] + u[i][2] * xc[2]);
        }
    }

    // -----------------------------------------------------------------------
    // 8. Convert eigen‑values to singular values; compute RMSD --------------
    // -----------------------------------------------------------------------
    for ev in &mut e {
        if *ev < 0.0 {
            *ev = 0.0;
        }
        *ev = ev.sqrt();
    }

    let mut d_s = e[2];
    if sigma < 0.0 {
        d_s = -d_s;
    }
    d_s += e[1] + e[0];

    if mode == 0 || mode == 2 {
        _rms1 = e0 - 2.0 * d_s;
        if _rms1 < 0.0 {
            _rms1 = 0.0;
        }
    }

    // Simple RMSD calculation as fallback
    if mode == 0 || mode == 2 {
        let mut sum_sq = 0.0;
        
        for i in 0..n {
            let transformed = [
                u[0][0] * x_f64[i][0] + u[0][1] * x_f64[i][1] + u[0][2] * x_f64[i][2] + t[0],
                u[1][0] * x_f64[i][0] + u[1][1] * x_f64[i][1] + u[1][2] * x_f64[i][2] + t[1],
                u[2][0] * x_f64[i][0] + u[2][1] * x_f64[i][1] + u[2][2] * x_f64[i][2] + t[2],
            ];

            for j in 0..3 {
                let diff = transformed[j] - y_f64[i][j];
                sum_sq += diff.powi(2);
            }
        }
        rms = (sum_sq / (n as f64)).sqrt();
    }

    // rms = rms1.sqrt();

    // -----------------------------------------------------------------------
    // 9. Convert back to f32 and deliver results ----------------------------
    // -----------------------------------------------------------------------
    
    // Convert f64 results back to f32
    let u_f32: [[f32; 3]; 3] = [
        [u[0][0] as f32, u[0][1] as f32, u[0][2] as f32],
        [u[1][0] as f32, u[1][1] as f32, u[1][2] as f32],
        [u[2][0] as f32, u[2][1] as f32, u[2][2] as f32],
    ];
    let t_f32: [f32; 3] = [t[0] as f32, t[1] as f32, t[2] as f32];
    let mut rms_f32 = rms as f32;
    if rms_f32.is_nan() {
        rms_f32 = f32::MAX;
    }
    
    Some((u_f32, t_f32, rms_f32))
}

/// Calculate TM-score between two aligned coordinate sets
/// 
/// # Arguments
/// * `coords1` - Moving coordinates (will be transformed)
/// * `coords2` - Target coordinates (reference)
/// * `rotation` - Rotation matrix from Kabsch alignment
/// * `translation` - Translation vector from Kabsch alignment
/// 
/// # Returns
/// TM-score normalized by the target structure length
pub fn calculate_tm_score(
    coords1: &[[f32; 3]], 
    coords2: &[[f32; 3]], 
    rotation: [[f32; 3]; 3], 
    translation: [f32; 3]
) -> f32 {
    calculate_tm_score_with_length(coords1, coords2, rotation, translation, coords2.len())
}

/// Calculate TM-score with custom normalization length
/// 
/// # Arguments
/// * `coords1` - Moving coordinates (will be transformed)
/// * `coords2` - Target coordinates (reference)
/// * `rotation` - Rotation matrix from Kabsch alignment
/// * `translation` - Translation vector from Kabsch alignment
/// * `norm_length` - Length to use for normalization (typically target length)
/// 
/// # Returns
/// TM-score normalized by the specified length
pub fn calculate_tm_score_with_length(
    coords1: &[[f32; 3]], 
    coords2: &[[f32; 3]], 
    rotation: [[f32; 3]; 3], 
    translation: [f32; 3],
    norm_length: usize
) -> f32 {
    assert_eq!(coords1.len(), coords2.len(), "Coordinate arrays must have the same length");
    
    let n = coords1.len();
    if n == 0 || norm_length == 0 {
        return 0.0;
    }
    
    // Calculate d0 normalization factor: d0 = 1.24 * (L_norm - 15)^(1/3) - 1.8
    let d0 = if norm_length <= 15 {
        0.5 // Minimum d0 for very short sequences
    } else {
        1.24 * ((norm_length - 15) as f32).powf(1.0/3.0) - 1.8
    }.max(0.5); // Ensure d0 is at least 0.5
    
    let d0_sq = d0 * d0;
    
    // Calculate TM-score
    let mut tm_sum = 0.0;
    
    for i in 0..n {
        // Transform coords1[i] using rotation and translation
        let transformed = matrix_vector_multiply(rotation, coords1[i]);
        let transformed = add_vec(transformed, translation);
        
        // Calculate squared distance to coords2[i]
        let dx = transformed[0] - coords2[i][0];
        let dy = transformed[1] - coords2[i][1];
        let dz = transformed[2] - coords2[i][2];
        let dist_sq = dx * dx + dy * dy + dz * dz;
        
        // Add to TM-score sum: 1 / (1 + (di/d0)^2)
        tm_sum += 1.0 / (1.0 + dist_sq / d0_sq);
    }
    
    // Normalize by target length
    tm_sum / (norm_length as f32)
}

/// Calculate Global Distance Test (GDT) scores
/// GDT measures the percentage of residues under different distance thresholds
/// 
/// # Arguments
/// * `coords1` - Moving coordinates (will be transformed)
/// * `coords2` - Target coordinates (reference)
/// * `rotation` - Rotation matrix from Kabsch alignment
/// * `translation` - Translation vector from Kabsch alignment
/// * `thresholds` - Distance thresholds to test (typically [1.0, 2.0, 4.0, 8.0])
/// 
/// # Returns
/// Vector of GDT percentages for each threshold
pub fn calculate_gdt_scores(
    coords1: &[[f32; 3]], 
    coords2: &[[f32; 3]], 
    rotation: [[f32; 3]; 3], 
    translation: [f32; 3],
    thresholds: &[f32]
) -> Vec<f32> {
    assert_eq!(coords1.len(), coords2.len(), "Coordinate arrays must have the same length");
    
    let n = coords1.len();
    if n == 0 {
        return vec![0.0; thresholds.len()];
    }
    
    // Calculate distances for all residues
    let mut distances = Vec::with_capacity(n);
    
    for i in 0..n {
        // Transform coords1[i] using rotation and translation
        let transformed = matrix_vector_multiply(rotation, coords1[i]);
        let transformed = add_vec(transformed, translation);
        
        // Calculate distance to coords2[i]
        let dx = transformed[0] - coords2[i][0];
        let dy = transformed[1] - coords2[i][1];
        let dz = transformed[2] - coords2[i][2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        
        distances.push(dist);
    }
    
    // Calculate GDT scores for each threshold
    thresholds.iter().map(|&threshold| {
        let count = distances.iter().filter(|&&dist| dist <= threshold).count();
        (count as f32) / (n as f32) * 100.0
    }).collect()
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kabsch_superimposer() {
        //
        let source = vec![
            Coordinate::new(6.994, 8.354, 42.405),
            Coordinate::new(9.429, 7.479, 48.266),
            Coordinate::new(5.547, 0.158, 42.050),
        ];
        
        let target1 = vec![
            Coordinate::new(-13.958, -1.741, -4.223),
            Coordinate::new(-12.833, 3.134, -7.780),
            Coordinate::new(-5.720, -2.218, -3.368),
        ];
        
        let target2 = vec![
            Coordinate::new(-4.924, 5.813, -9.485),
            Coordinate::new(-0.499, 10.073, -8.059),
            Coordinate::new(-0.792, 0.658, -4.430),
        ];

        let mut superimposer = KabschSuperimposer::new();
        superimposer.set_atoms(&source, &target1);
        let start = std::time::Instant::now();
        println!("RMSD after superimposition: {}", superimposer.get_rms());
        println!("TM-score: {}", superimposer.get_tm_score());
        let duration = start.elapsed();
        println!("Time taken for superimposition: {:?}", duration);
        println!("Transformed coordinates: {:?}", superimposer.get_transformed());
        println!("Rotation matrix: {:?}", superimposer.get_rotran().0);
        println!("Translation vector: {:?}", superimposer.get_rotran().1);
        assert!(superimposer.get_rms() < 0.2);
        
        superimposer.set_atoms(&source, &target2);
        assert!(superimposer.get_rms() < 0.2);
        let start = std::time::Instant::now();
        println!("RMSD after superimposition: {}", superimposer.get_rms());
        println!("TM-score: {}", superimposer.get_tm_score());
        let duration = start.elapsed();
        println!("Time taken for superimposition: {:?}", duration);
        println!("Transformed coordinates: {:?}", superimposer.get_transformed());
        println!("Rotation matrix: {:?}", superimposer.get_rotran().0);
        println!("Translation vector: {:?}", superimposer.get_rotran().1);
    }

    #[test]
    fn test_identical_structures() {
        let coords = vec![
            Coordinate::new(1.0, 2.0, 3.0),
            Coordinate::new(4.0, 5.0, 6.0),
            Coordinate::new(7.0, 8.0, 9.0),
        ];

        let mut superimposer = KabschSuperimposer::new();
        superimposer.set_atoms(&coords, &coords);
        assert!(superimposer.get_rms() < 1e-6);
        
        // TM-score should be 1.0 for identical structures
        let tm_score = superimposer.get_tm_score();
        println!("TM-score for identical structures: {}", tm_score);
        assert!((tm_score - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_tm_score_calculation() {
        // Test TM-score calculation with known structures
        let coords1 = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        
        let coords2 = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        
        // Identity transformation
        let rotation = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let translation = [0.0, 0.0, 0.0];
        
        let tm_score = calculate_tm_score(&coords1, &coords2, rotation, translation);
        println!("TM-score for identical aligned structures: {}", tm_score);
        assert!((tm_score - 1.0).abs() < 1e-6);
        
        // Test with small displacement
        let coords2_displaced = vec![
            [0.1, 0.0, 0.0],
            [1.1, 0.0, 0.0],
            [2.1, 0.0, 0.0],
            [3.1, 0.0, 0.0],
        ];
        
        let tm_score_displaced = calculate_tm_score(&coords1, &coords2_displaced, rotation, translation);
        println!("TM-score with 0.1Å displacement: {}", tm_score_displaced);
        assert!(tm_score_displaced > 0.8); // Should still be high
        assert!(tm_score_displaced < 1.0); // But less than perfect
    }

    #[test]
    fn test_gdt_calculation() {
        let coords1 = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        
        let coords2 = vec![
            [0.5, 0.0, 0.0], // 0.5Å distance
            [1.0, 1.5, 0.0], // 1.5Å distance
            [2.0, 0.0, 2.5], // 2.5Å distance
            [3.0, 0.0, 5.0], // 5.0Å distance
        ];
        
        let rotation = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let translation = [0.0, 0.0, 0.0];
        let thresholds = [1.0, 2.0, 4.0, 8.0];
        
        let gdt_scores = calculate_gdt_scores(&coords1, &coords2, rotation, translation, &thresholds);
        println!("GDT scores: {:?}", gdt_scores);
        
        // Expected: 25% under 1Å, 50% under 2Å, 75% under 4Å, 100% under 8Å
        assert!((gdt_scores[0] - 25.0).abs() < 1e-6);
        assert!((gdt_scores[1] - 50.0).abs() < 1e-6);
        assert!((gdt_scores[2] - 75.0).abs() < 1e-6);
        assert!((gdt_scores[3] - 100.0).abs() < 1e-6);
    }
}