// File: qcp.rs
// Created: 2024-05-09 17:50:26
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

// Original code from BioPython. Implemented in Rust.
// reference: https://github.com/biopython/biopython/blob/master/Bio/PDB/qcprot.py
// C reference: https://theobald.brandeis.edu/qcp/

use crate::structure::coordinate::Coordinate;


#[derive(Debug)]
pub struct QCPSuperimposer {
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

impl QCPSuperimposer {
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

        let com_coords = mean(&coords);
        let com_ref = mean(&reference_coords);

        let centered_coords: Vec<[f32; 3]> = coords.iter()
            .map(|&coord| {
                [
                    coord[0] - com_coords[0],
                    coord[1] - com_coords[1],
                    coord[2] - com_coords[2]
                ]
            })
            .collect();

        let centered_ref: Vec<[f32; 3]> = reference_coords.iter()
            .map(|&coord| {
                [
                    coord[0] - com_ref[0],
                    coord[1] - com_ref[1],
                    coord[2] - com_ref[2]
                ]
            })
            .collect();

        // Get rotation from QCP (aligns centered_coords to centered_ref)
        let (_rms_from_qcp, rot, _) = qcp(&centered_coords, &centered_ref, self.natoms);

        self.rot = Some(rot);
        
        // Correct translation calculation
        let rotated_com_coords = rotate(com_coords, rot);
        let tran = [
            com_ref[0] - rotated_com_coords[0],
            com_ref[1] - rotated_com_coords[1],
            com_ref[2] - rotated_com_coords[2]
        ];
        self.tran = Some(tran);

        // Recalculate rmsd with rot, tran
        self.transformed_coords = Some(
            coords.iter()
                .map(|&coord| {
                    let rotated = rotate(coord, rot);
                    [
                        rotated[0] + tran[0],
                        rotated[1] + tran[1],
                        rotated[2] + tran[2]
                    ]
                })
                .collect()
        );
        
        // Calculate final RMSD with transformed coordinates and reference coordinates
        let transformed_coords = self.transformed_coords.clone().unwrap();
        let diff: Vec<f32> = transformed_coords.iter()
            .zip(reference_coords.iter())
            .map(|(&c1, &c2)| {
                (c1[0] - c2[0]).powi(2) + (c1[1] - c2[1]).powi(2) + (c1[2] - c2[2]).powi(2)
            })
            .collect();
        self.rms = Some((diff.iter().sum::<f32>() / self.natoms as f32).sqrt());        
        
        // Calculate TM-score
        self.tm_score = Some(calculate_tm_score(&coords, &reference_coords, rot, tran));
        
    }

    pub fn get_transformed(&mut self) -> Vec<[f32; 3]> {
        // Always recalculate transformed_coords
        let coords = self.coords.clone().unwrap();
        let rot = self.rot.unwrap();
        let tran = self.tran.unwrap();

        self.transformed_coords = Some(
            coords.iter()
                .map(|&coord| {
                    let rotated = rotate(coord, rot);
                    [
                        rotated[0] + tran[0],
                        rotated[1] + tran[1],
                        rotated[2] + tran[2]
                    ]
                })
                .collect()
        );
    self.transformed_coords.clone().unwrap()
    }

    pub fn get_rotran(&self) -> ([[f32; 3]; 3], [f32; 3]) {
        (self.rot.unwrap(), self.tran.unwrap())
    }

    pub fn get_init_rms(&mut self) -> f32 {
        if self.init_rms.is_none() {
            let coords = self.coords.clone().unwrap();
            let reference_coords = self.reference_coords.clone().unwrap();

            let diff: Vec<f32> = coords.iter()
                .zip(reference_coords.iter())
                .map(|(&c1, &c2)| {
                    (c1[0] - c2[0]).powi(2) + (c1[1] - c2[1]).powi(2) + (c1[2] - c2[2]).powi(2)
                })
                .collect();

            self.init_rms = Some((diff.iter().sum::<f32>() / self.natoms as f32).sqrt());
        }

        self.init_rms.unwrap()
    }

    pub fn get_rms(&self) -> f32 {
        self.rms.unwrap()
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
}

/// Calculate TM-score between two aligned coordinate sets
/// 
/// # Arguments
/// * `coords1` - Moving coordinates (will be transformed)
/// * `coords2` - Target coordinates (reference)
/// * `rotation` - Rotation matrix from alignment
/// * `translation` - Translation vector from alignment
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
/// * `rotation` - Rotation matrix from alignment
/// * `translation` - Translation vector from alignment
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
        let transformed = rotate(coords1[i], rotation);
        let transformed = [
            transformed[0] + translation[0],
            transformed[1] + translation[1],
            transformed[2] + translation[2]
        ];
        
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

fn mean(coords: &[[f32; 3]]) -> [f32; 3] {
    let sum: [f32; 3] = coords.iter().fold([0.0; 3], |acc, &coord| {
        [acc[0] + coord[0], acc[1] + coord[1], acc[2] + coord[2]]
    });

    [sum[0] / coords.len() as f32, sum[1] / coords.len() as f32, sum[2] / coords.len() as f32]
}

fn rotate(coord: [f32; 3], rot: [[f32; 3]; 3]) -> [f32; 3] {
    [
        coord[0] * rot[0][0] + coord[1] * rot[0][1] + coord[2] * rot[0][2],
        coord[0] * rot[1][0] + coord[1] * rot[1][1] + coord[2] * rot[1][2],
        coord[0] * rot[2][0] + coord[1] * rot[2][1] + coord[2] * rot[2][2]
    ]
}

fn qcp(coords1: &[[f32; 3]], coords2: &[[f32; 3]], natoms: usize) -> (f32, [[f32; 3]; 3], [f32; 4]) {
    let g1 = trace(dot(coords2, &transpose(coords2)));
    let g2 = trace(dot(coords1, &transpose(coords1)));
    let a = dot(&transpose(coords2), coords1);
    let e0 = (g1 + g2) * 0.5;

    let sxx = a[0][0];
    let sxy = a[0][1];
    let sxz = a[0][2];
    let syx = a[1][0];
    let syy = a[1][1];
    let syz = a[1][2];
    let szx = a[2][0];
    let szy = a[2][1];
    let szz = a[2][2];

    let sxx2 = sxx * sxx;
    let syy2 = syy * syy;
    let szz2 = szz * szz;
    let sxy2 = sxy * sxy;
    let syz2 = syz * syz;
    let sxz2 = sxz * sxz;
    let syx2 = syx * syx;
    let szy2 = szy * szy;
    let szx2 = szx * szx;

    let syz_szy_m_syy_szz2 = 2.0 * (syz * szy - syy * szz);
    let sxx2_syy2_szz2_syz2_szy2 = syy2 + szz2 - sxx2 + syz2 + szy2;

    let c2 = -2.0 * (sxx2 + syy2 + szz2 + sxy2 + syx2 + sxz2 + szx2 + syz2 + szy2);
    let c1 = 8.0 * (sxx * syz * szy + syy * szx * sxz + szz * sxy * syx - sxx * syy * szz - syz * szx * sxy - szy * syx * sxz);

    let sxz_p_szx = sxz + szx;
    let syz_p_szy = syz + szy;
    let sxy_p_syx = sxy + syx;
    let syz_m_szy = syz - szy;
    let sxz_m_szx = sxz - szx;
    let sxy_m_syx = sxy - syx;
    let sxx_p_syy = sxx + syy;
    let sxx_m_syy = sxx - syy;
    let sxy2_sxz2_syx2_szx2 = sxy2 + sxz2 - syx2 - szx2;

    let neg_sxz_p_szx = -sxz_p_szx;
    let neg_sxz_m_szx = -sxz_m_szx;
    let neg_sxy_m_syx = -sxy_m_syx;
    let sxx_p_syy_p_szz = sxx_p_syy + szz;

    let c0 = sxy2_sxz2_syx2_szx2 * sxy2_sxz2_syx2_szx2
        + (sxx2_syy2_szz2_syz2_szy2 + syz_szy_m_syy_szz2)
        * (sxx2_syy2_szz2_syz2_szy2 - syz_szy_m_syy_szz2)
        + (neg_sxz_p_szx * (syz_m_szy) + (sxy_m_syx) * (sxx_m_syy - szz))
        * (neg_sxz_m_szx * (syz_p_szy) + (sxy_m_syx) * (sxx_m_syy + szz))
        + (neg_sxz_p_szx * (syz_p_szy) - (sxy_p_syx) * (sxx_p_syy - szz))
        * (neg_sxz_m_szx * (syz_m_szy) - (sxy_p_syx) * sxx_p_syy_p_szz)
        + ((sxy_p_syx) * (syz_p_szy) + (sxz_p_szx) * (sxx_m_syy + szz))
        * (neg_sxy_m_syx * (syz_m_szy) + (sxz_p_szx) * sxx_p_syy_p_szz)
        + ((sxy_p_syx) * (syz_m_szy) + (sxz_m_szx) * (sxx_m_syy - szz))
        * (neg_sxy_m_syx * (syz_p_szy) + (sxz_m_szx) * (sxx_p_syy - szz));

    let mut mx_eigenv = e0; // starting guess (x in eqs above)
    let eval_prec = 1e-11; // convergence criterion
    let mut _converged = false;
    let iteration = 5; // 5 is known to be enough according to BioPython but following the original
    for _ in 0..iteration {
        let oldg = mx_eigenv;

        let x2 = mx_eigenv * mx_eigenv;
        let b = (x2 + c2) * mx_eigenv;
        let a = b + c1;

        let f = a * mx_eigenv + c0;
        let f_prime = 2.0 * x2 * mx_eigenv + b + a;

        let delta = f / (f_prime + eval_prec); // avoid division by zero
        mx_eigenv = (mx_eigenv - delta).abs();
        if (mx_eigenv - oldg).abs() < (eval_prec * mx_eigenv) {
            _converged = true;
            break; // convergence
        }
    }
    // if !converged {
    //     println!("Newton-Rhapson did not converge after 50 iterations");
    // }

    let rmsd = (2.0 * (e0 - mx_eigenv).abs() / natoms as f32).sqrt();

    let a11 = sxx_p_syy + szz - mx_eigenv;
    let a12 = syz_m_szy;
    let a13 = neg_sxz_m_szx;
    let a14 = sxy_m_syx;
    let a21 = syz_m_szy;
    let a22 = sxx_m_syy - szz - mx_eigenv;
    let a23 = sxy_p_syx;
    let a24 = sxz_p_szx;
    let a31 = a13;
    let a32 = a23;
    let a33 = syy - sxx - szz - mx_eigenv;
    let a34 = syz_p_szy;
    let a41 = a14;
    let a42 = a24;
    let a43 = a34;
    let a44 = szz - sxx_p_syy - mx_eigenv;

    let a3344_4334 = a33 * a44 - a43 * a34;
    let a3244_4234 = a32 * a44 - a42 * a34;
    let a3243_4233 = a32 * a43 - a42 * a33;
    let a3143_4133 = a31 * a43 - a41 * a33;
    let a3144_4134 = a31 * a44 - a41 * a34;
    let a3142_4132 = a31 * a42 - a41 * a32;

    let mut q1 = a22 * a3344_4334 - a23 * a3244_4234 + a24 * a3243_4233;
    let mut q2 = -a21 * a3344_4334 + a23 * a3144_4134 - a24 * a3143_4133;
    let mut q3 = a21 * a3244_4234 - a22 * a3144_4134 + a24 * a3142_4132;
    let mut q4 = -a21 * a3243_4233 + a22 * a3143_4133 - a23 * a3142_4132;

    let qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

    let evec_prec = 1e-6;
    if qsqr < evec_prec {
        q1 = a12 * a3344_4334 - a13 * a3244_4234 + a14 * a3243_4233;
        q2 = -a11 * a3344_4334 + a13 * a3144_4134 - a14 * a3143_4133;
        q3 = a11 * a3244_4234 - a12 * a3144_4134 + a14 * a3142_4132;
        q4 = -a11 * a3243_4233 + a12 * a3143_4133 - a13 * a3142_4132;

        let qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

        if qsqr < evec_prec {
            let a1324_1423 = a13 * a24 - a14 * a23;
            let a1224_1422 = a12 * a24 - a14 * a22;
            let a1223_1322 = a12 * a23 - a13 * a22;
            let a1124_1421 = a11 * a24 - a14 * a21;
            let a1123_1321 = a11 * a23 - a13 * a21;
            let a1122_1221 = a11 * a22 - a12 * a21;

            q1 = a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
            q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
            q3 = a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
            q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;

            let qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

            if qsqr < evec_prec {
                q1 = a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                q3 = a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;

                let qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

                if qsqr < evec_prec {
                    let rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
                    return (rmsd, rot, [q1, q2, q3, q4]);
                }
            }
        }
    }

    let normq = qsqr.sqrt();
    q1 /= normq;
    q2 /= normq;
    q3 /= normq;
    q4 /= normq;

    let a2 = q1 * q1;
    let x2 = q2 * q2;
    let y2 = q3 * q3;
    let z2 = q4 * q4;

    let xy = q2 * q3;
    let az = q1 * q4;
    let zx = q4 * q2;
    let ay = q1 * q3;
    let yz = q3 * q4;
    let ax = q1 * q2;

    let mut rot = [[0.0; 3]; 3];

    rot[0][0] = a2 + x2 - y2 - z2;
    rot[0][1] = 2.0 * (xy + az);
    rot[0][2] = 2.0 * (zx - ay);
    rot[1][0] = 2.0 * (xy - az);
    rot[1][1] = a2 - x2 + y2 - z2;
    rot[1][2] = 2.0 * (yz + ax);
    rot[2][0] = 2.0 * (zx + ay);
    rot[2][1] = 2.0 * (yz - ax);
    rot[2][2] = a2 - x2 - y2 + z2;


    (rmsd, rot, [q1, q2, q3, q4])
}



fn transpose(matrix: &[[f32; 3]]) -> [[f32; 3]; 3] {
    let mut transposed = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            transposed[j][i] = matrix[i][j];
        }
    }
    transposed
}

fn dot(matrix1: &[[f32; 3]], matrix2: &[[f32; 3]]) -> [[f32; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] = matrix1[i][0] * matrix2[0][j]
                + matrix1[i][1] * matrix2[1][j]
                + matrix1[i][2] * matrix2[2][j];
        }
    }
    result
}

fn trace(matrix: [[f32; 3]; 3]) -> f32 {
    matrix[0][0] + matrix[1][1] + matrix[2][2]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qcp_superimposer() {
        let source = vec![
            Coordinate::new(6.994, 8.354, 42.405),
            Coordinate::new(9.429, 7.479, 48.266),
            Coordinate::new(5.547, 0.158, 42.050),
        ];
        // E.coli
        let target1 = vec![
            Coordinate::new(-13.958, -1.741, -4.223),
            Coordinate::new(-12.833, 3.134, -7.780),
            Coordinate::new(-5.720, -2.218, -3.368),
        ];
        // Human
        let target2 = vec![
            Coordinate::new(-4.924, 5.813, -9.485),
            Coordinate::new(-0.499, 10.073, -8.059),
            Coordinate::new(-0.792, 0.658, -4.430),
        ];

        let mut superimposer = QCPSuperimposer::new();
        superimposer.set_atoms(&source, &target1);
        // assert!(superimposer.get_rms() < 0.2);
        let start  = std::time::Instant::now();
        println!("RMSD after superimposition: {}", superimposer.get_rms());
        println!("TM-score: {}", superimposer.get_tm_score());
        let duration = start.elapsed();
        println!("Time taken for superimposition: {:?}", duration);
        println!("Transformed coordinates: {:?}", superimposer.get_transformed());
        println!("Rotation matrix: {:?}", superimposer.get_rotran().0);
        println!("Translation vector: {:?}", superimposer.get_rotran().1);
        superimposer.set_atoms(&source, &target2);
        // assert!(superimposer.get_rms() < 0.2);
        let start = std::time::Instant::now();
        println!("RMSD after superimposition: {}", superimposer.get_rms());
        println!("TM-score: {}", superimposer.get_tm_score());
        let duration = start.elapsed();
        println!("Time taken for superimposition: {:?}", duration);
        println!("Transformed coordinates: {:?}", superimposer.get_transformed());
        println!("Rotation matrix: {:?}", superimposer.get_rotran().0);
        println!("Translation vector: {:?}", superimposer.get_rotran().1);
        assert!(superimposer.get_rms() < 0.2);
    }

    #[test]
    fn test_tm_score_calculation() {
        // Test TM-score calculation with identical structures
        let coords = vec![
            Coordinate::new(1.0, 2.0, 3.0),
            Coordinate::new(4.0, 5.0, 6.0),
            Coordinate::new(7.0, 8.0, 9.0),
        ];

        let mut superimposer = QCPSuperimposer::new();
        superimposer.set_atoms(&coords, &coords);
        
        println!("QCP RMSD for identical structures: {}", superimposer.get_rms());
        println!("QCP TM-score for identical structures: {}", superimposer.get_tm_score());
        
        // TM-score should be 1.0 for identical structures
        let tm_score = superimposer.get_tm_score();
        assert!((tm_score - 1.0).abs() < 1e-6);
        assert!(superimposer.get_rms() < 1e-6);
    }
}
