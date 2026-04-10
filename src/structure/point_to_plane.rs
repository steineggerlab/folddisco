// File: point_to_plane.rs
// Symmetric point-to-plane ICP superimposer.
// Coordinate convention: even indices (0,2,4,...) = Cα, odd indices (1,3,5,...) = Cβ.
// API is compatible with KabschSuperimposer.
//
// Algorithm:
//   1. Kabsch warm-start      – provides a good initial rotation/translation.
//   2. Iterative linearized refinement – minimises the Rusinkiewicz (2019)
//      symmetric point-to-plane objective:
//        E = Σ_i [ (p_i − q_i) · (np_i + nq_i) ]²
//      where p_i = transformed Cβ_i, q_i = reference Cβ_i,
//            np_i = transformed (Cβ_i − Cα_i), nq_i = reference (Cβ_i − Cα_i).
//
// The "rms" field stores sqrt( mean E_i ) after convergence.

use crate::structure::coordinate::Coordinate;
use crate::structure::kabsch::kabsch;

// ─── tunables ──────────────────────────────────────────────────────────────
const MAX_ITER: usize       = 50;
const CONV_TOL: f64         = 1e-10; // convergence: ‖[δα, δt]‖² < tol
const NORMAL_EPS: f64       = 1e-8;  // skip pairs with ‖np + nq‖² < eps
const GAUSS_EPS: f64        = 1e-12; // pivot threshold for Gaussian elimination

// ─── public struct ─────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct PointToPlaneSupimposer {
    pub reference_coords:  Option<Vec<[f32; 3]>>,
    pub coords:            Option<Vec<[f32; 3]>>,
    pub transformed_coords: Option<Vec<[f32; 3]>>,
    pub rot:               Option<[[f32; 3]; 3]>,
    pub tran:              Option<[f32; 3]>,
    /// sqrt( mean symmetric-point-to-plane residual² ) after superposition.
    pub rms:               Option<f32>,
    /// sqrt( mean symmetric-point-to-plane residual² ) before superposition.
    pub init_rms:          Option<f32>,
    pub natoms:            usize,
}

impl PointToPlaneSupimposer {
    pub fn new() -> Self {
        Self {
            reference_coords:  None,
            coords:            None,
            transformed_coords: None,
            rot:               None,
            tran:              None,
            rms:               None,
            init_rms:          None,
            natoms:            0,
        }
    }

    /// Convenience: set from `Coordinate` slices, run, and cache rms.
    pub fn set_atoms(&mut self, fixed: &[Coordinate], moving: &[Coordinate]) {
        assert_eq!(fixed.len(), moving.len(), "fixed and moving differ in size");
        let fixed_c:  Vec<[f32; 3]> = fixed.iter().map(|c| c.to_array()).collect();
        let moving_c: Vec<[f32; 3]> = moving.iter().map(|c| c.to_array()).collect();
        self.set(&fixed_c, &moving_c);
        self.run();
        self.rms = Some(self.get_rms());
    }

    /// Set coordinate arrays (reference = fixed, coords = moving).
    /// `coords.len()` must be even (CA/CB interleaved).
    pub fn set(&mut self, reference_coords: &[[f32; 3]], coords: &[[f32; 3]]) {
        assert_eq!(reference_coords.len(), coords.len(),
                   "reference and mobile coords differ in length");
        assert!(coords.len() >= 2 && coords.len() % 2 == 0,
                "natoms must be even (interleaved Cα/Cβ pairs)");
        self.reference_coords  = Some(reference_coords.to_vec());
        self.coords            = Some(coords.to_vec());
        self.natoms            = coords.len();
        self.transformed_coords = None;
        self.rot               = None;
        self.tran              = None;
        self.rms               = None;
        self.init_rms          = None;
    }

    /// Compute optimal superposition and store rot, tran, rms.
    pub fn run(&mut self) {
        let coords    = self.coords.clone().unwrap();
        let reference = self.reference_coords.clone().unwrap();
        let n_res     = self.natoms / 2;

        // ── 1. Kabsch warm-start ──────────────────────────────────────────
        let (rot0, tran0, _) = kabsch(&coords, &reference, 2)
            .unwrap_or(([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], [0.,0.,0.], f32::MAX));

        let mut rot64  = f32_to_f64_mat3(rot0);
        let mut tran64 = [tran0[0] as f64, tran0[1] as f64, tran0[2] as f64];

        // pre-convert references to f64 (done once)
        let ref64: Vec<[f64; 3]> = reference.iter()
            .map(|c| [c[0] as f64, c[1] as f64, c[2] as f64])
            .collect();
        let src64: Vec<[f64; 3]> = coords.iter()
            .map(|c| [c[0] as f64, c[1] as f64, c[2] as f64])
            .collect();

        // ── 2. Iterative linearized point-to-plane refinement ────────────
        for _ in 0..MAX_ITER {
            // transform all source atoms with current (rot64, tran64)
            let tf: Vec<[f64; 3]> = src64.iter()
                .map(|&c| {
                    let r = matvec64(rot64, c);
                    [r[0]+tran64[0], r[1]+tran64[1], r[2]+tran64[2]]
                })
                .collect();

            // accumulate 6×6 normal-equation system
            let mut ata = [[0.0f64; 6]; 6];
            let mut atb = [0.0f64; 6];

            for i in 0..n_res {
                let ca_t = tf[2*i];
                let cb_t = tf[2*i + 1];
                let ca_r = ref64[2*i];
                let cb_r = ref64[2*i + 1];

                // np = transformed CA→CB direction
                let np = sub64(cb_t, ca_t);
                // nq = reference CA→CB direction
                let nq = sub64(cb_r, ca_r);
                // symmetric combined normal
                let m  = add64(np, nq);

                // skip degenerate (anti-parallel) cases
                let norm_sq = dot64(m, m);
                if norm_sq < NORMAL_EPS { continue; }

                // residual  r₀ = (p − q) · m
                let d_pq = sub64(cb_t, cb_r); // p − q
                let r0   = dot64(d_pq, m);

                // linearised row: a = [np×(p−q) + p×m,  m]
                let c1  = cross64(np, d_pq);       // np × (p−q)
                let c2  = cross64(cb_t, m);         // p  × m
                let row = [
                    c1[0]+c2[0], c1[1]+c2[1], c1[2]+c2[2],
                    m[0], m[1], m[2],
                ];

                // AᵀA·x = −Aᵀr₀
                for j in 0..6 {
                    for k in 0..6 { ata[j][k] += row[j] * row[k]; }
                    atb[j] -= row[j] * r0;
                }
            }

            // solve the 6×6 system
            let x = match solve_6x6(ata, atb) {
                Some(v) => v,
                None    => break,
            };

            let dalpha = [x[0], x[1], x[2]];
            let dt     = [x[3], x[4], x[5]];

            // incremental rotation from δα via Rodrigues
            let dr = rodrigues64(dalpha);

            // compose: R_new = dR · R_old,  t_new = dR · t_old + δt
            rot64  = matmul64(dr, rot64);
            let rt = matvec64(dr, tran64);
            tran64 = [rt[0]+dt[0], rt[1]+dt[1], rt[2]+dt[2]];

            // convergence check
            let dn: f64 = dalpha.iter().chain(dt.iter()).map(|v| v*v).sum();
            if dn < CONV_TOL { break; }
        }

        // ── 3. Store results ──────────────────────────────────────────────
        let rot  = f64_to_f32_mat3(rot64);
        let tran = [tran64[0] as f32, tran64[1] as f32, tran64[2] as f32];
        self.rot  = Some(rot);
        self.tran = Some(tran);

        let transformed: Vec<[f32; 3]> = coords.iter()
            .map(|&c| {
                let r = matvec_f32(rot, c);
                [r[0]+tran[0], r[1]+tran[1], r[2]+tran[2]]
            })
            .collect();
        self.transformed_coords = Some(transformed.clone());

        // point-to-plane score on final aligned coords
        self.rms = Some(calc_ptp_score(&reference, &transformed, n_res));
    }

    pub fn get_transformed(&mut self) -> Vec<[f32; 3]> {
        if self.transformed_coords.is_none() {
            let coords = self.coords.clone().unwrap();
            let rot    = self.rot.unwrap();
            let tran   = self.tran.unwrap();
            self.transformed_coords = Some(
                coords.iter()
                    .map(|&c| {
                        let r = matvec_f32(rot, c);
                        [r[0]+tran[0], r[1]+tran[1], r[2]+tran[2]]
                    })
                    .collect()
            );
        }
        self.transformed_coords.clone().unwrap()
    }

    pub fn get_rotran(&self) -> ([[f32; 3]; 3], [f32; 3]) {
        (self.rot.unwrap(), self.tran.unwrap())
    }

    /// Point-to-plane score on the un-superposed coordinates.
    pub fn get_init_rms(&mut self) -> f32 {
        if self.init_rms.is_none() {
            let coords    = self.coords.clone().unwrap();
            let reference = self.reference_coords.clone().unwrap();
            let n_res     = self.natoms / 2;
            self.init_rms = Some(calc_ptp_score(&reference, &coords, n_res));
        }
        self.init_rms.unwrap()
    }

    pub fn get_rms(&self) -> f32 {
        self.rms.unwrap()
    }
}

// ─── scoring helper ────────────────────────────────────────────────────────

/// sqrt( mean symmetric-pt-to-plane residual² ) over `n_res` residue pairs.
/// Coords are interleaved: index 2i = Cα, 2i+1 = Cβ.
fn calc_ptp_score(reference: &[[f32; 3]], transformed: &[[f32; 3]], n_res: usize) -> f32 {
    let sum: f32 = (0..n_res)
        .map(|i| point_to_plane_distance(
            reference[2*i], reference[2*i+1],
            transformed[2*i], transformed[2*i+1],
        ))
        .sum();
    (sum / n_res as f32).sqrt()
}

// ─── point-to-plane residual ───────────────────────────────────────────────

/// Symmetric point-to-plane squared residual (Rusinkiewicz 2019):
///   [ (Cβ_ref − Cβ_target) · (np + nq) ]²
/// where np = Cβ_target − Cα_target, nq = Cβ_ref − Cα_ref.
///
/// Returns 0 if the combined normal is degenerate (anti-parallel CB vectors).
pub fn point_to_plane_distance(
    ca_ref:    [f32; 3],
    cb_ref:    [f32; 3],
    ca_target: [f32; 3],
    cb_target: [f32; 3],
) -> f32 {
    let np = [
        cb_target[0] - ca_target[0],
        cb_target[1] - ca_target[1],
        cb_target[2] - ca_target[2],
    ];
    let nq = [
        cb_ref[0] - ca_ref[0],
        cb_ref[1] - ca_ref[1],
        cb_ref[2] - ca_ref[2],
    ];
    let m = [np[0]+nq[0], np[1]+nq[1], np[2]+nq[2]];

    // epsilon guard for degenerate (anti-parallel) normals
    let norm_sq = m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
    if norm_sq < 1e-8_f32 { return 0.0; }

    let pq = [
        cb_ref[0] - cb_target[0],
        cb_ref[1] - cb_target[1],
        cb_ref[2] - cb_target[2],
    ];
    (pq[0]*m[0] + pq[1]*m[1] + pq[2]*m[2]).powi(2)
}

// ═══════════════════════════════════════════════════════════════════════════
// Internal math helpers (all f64 for numerical stability)
// ═══════════════════════════════════════════════════════════════════════════

#[inline] fn sub64(a: [f64;3], b: [f64;3]) -> [f64;3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }
#[inline] fn add64(a: [f64;3], b: [f64;3]) -> [f64;3] { [a[0]+b[0], a[1]+b[1], a[2]+b[2]] }
#[inline] fn dot64(a: [f64;3], b: [f64;3]) -> f64      { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }

#[inline]
fn cross64(a: [f64;3], b: [f64;3]) -> [f64;3] {
    [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    ]
}

#[inline]
fn matvec64(m: [[f64;3];3], v: [f64;3]) -> [f64;3] {
    [
        m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
        m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
        m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2],
    ]
}

#[inline]
fn matvec_f32(m: [[f32;3];3], v: [f32;3]) -> [f32;3] {
    [
        m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
        m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
        m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2],
    ]
}

fn matmul64(a: [[f64;3];3], b: [[f64;3];3]) -> [[f64;3];3] {
    let mut c = [[0.0f64;3];3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 { c[i][j] += a[i][k] * b[k][j]; }
        }
    }
    c
}

#[inline]
fn f32_to_f64_mat3(m: [[f32;3];3]) -> [[f64;3];3] {
    [
        [m[0][0] as f64, m[0][1] as f64, m[0][2] as f64],
        [m[1][0] as f64, m[1][1] as f64, m[1][2] as f64],
        [m[2][0] as f64, m[2][1] as f64, m[2][2] as f64],
    ]
}

#[inline]
fn f64_to_f32_mat3(m: [[f64;3];3]) -> [[f32;3];3] {
    [
        [m[0][0] as f32, m[0][1] as f32, m[0][2] as f32],
        [m[1][0] as f32, m[1][1] as f32, m[1][2] as f32],
        [m[2][0] as f32, m[2][1] as f32, m[2][2] as f32],
    ]
}

/// Rodrigues rotation matrix for rotation vector `v` (angle = ‖v‖, axis = v/‖v‖).
fn rodrigues64(v: [f64;3]) -> [[f64;3];3] {
    let angle = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt();
    if angle < 1e-14 {
        return [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    }
    let (s, c)  = angle.sin_cos();
    let t       = 1.0 - c;
    let ax      = v[0] / angle;
    let ay      = v[1] / angle;
    let az      = v[2] / angle;
    [
        [t*ax*ax + c,    t*ax*ay - s*az, t*ax*az + s*ay],
        [t*ax*ay + s*az, t*ay*ay + c,    t*ay*az - s*ax],
        [t*ax*az - s*ay, t*ay*az + s*ax, t*az*az + c   ],
    ]
}

/// Gaussian elimination with partial pivoting for a 6×6 system: A·x = b.
/// Returns `None` if the system is (near-)singular.
fn solve_6x6(mut a: [[f64;6];6], mut b: [f64;6]) -> Option<[f64;6]> {
    const N: usize = 6;
    let mut perm: [usize; N] = [0,1,2,3,4,5];

    for col in 0..N {
        // find pivot
        let mut max_val = a[col][col].abs();
        let mut max_row = col;
        for row in (col+1)..N {
            let v = a[row][col].abs();
            if v > max_val { max_val = v; max_row = row; }
        }
        if max_val < GAUSS_EPS { return None; }

        // swap rows
        if max_row != col {
            a.swap(col, max_row);
            b.swap(col, max_row);
            perm.swap(col, max_row);
        }

        let pivot = a[col][col];
        for row in (col+1)..N {
            let factor = a[row][col] / pivot;
            for k in col..N { a[row][k] -= factor * a[col][k]; }
            b[row] -= factor * b[col];
        }
    }

    // back-substitution
    let mut x = [0.0f64; N];
    for i in (0..N).rev() {
        let mut s = b[i];
        for j in (i+1)..N { s -= a[i][j] * x[j]; }
        x[i] = s / a[i][i];
    }
    Some(x)
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::kabsch::KabschSuperimposer;

    /// Build a trivial test: source = reference rotated 30° around Z + small translation.
    fn make_test_coords(n_res: usize) -> (Vec<[f32;3]>, Vec<[f32;3]>) {
        use std::f32::consts::PI;
        let theta = 30.0_f32 * PI / 180.0;
        let (s, c) = theta.sin_cos();
        let _rot = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0_f32]];
        let tran = [1.5_f32, -0.5, 0.3];

        let mut reference = Vec::with_capacity(n_res * 2);
        let mut mobile    = Vec::with_capacity(n_res * 2);

        for i in 0..n_res {
            let fi = i as f32;
            let ca_ref: [f32;3] = [fi,     fi * 0.5, fi * 0.1];
            let cb_ref: [f32;3] = [fi+0.5, fi * 0.5 + 1.0, fi * 0.1 + 0.5];

            // mobile = R⁻¹(ref − t)
            let derot = |p: [f32;3]| -> [f32;3] {
                let q = [p[0]-tran[0], p[1]-tran[1], p[2]-tran[2]];
                [c*q[0]+s*q[1], -s*q[0]+c*q[1], q[2]]
            };

            reference.push(ca_ref);
            reference.push(cb_ref);
            mobile.push(derot(ca_ref));
            mobile.push(derot(cb_ref));
        }
        (reference, mobile)
    }

    #[test]
    fn test_ptp_known_rotation() {
        let (reference, mobile) = make_test_coords(6);
        let mut sup = PointToPlaneSupimposer::new();
        sup.set(&reference, &mobile);
        sup.run();

        let rms = sup.get_rms();
        println!("point-to-plane rms after superposition: {rms:.6}");
        // After exact inverse rotation the score should be near zero.
        assert!(rms < 1e-3, "rms={rms} unexpectedly large");
    }

    #[test]
    fn test_get_rotran_and_transformed() {
        let (reference, mobile) = make_test_coords(8);
        let mut sup = PointToPlaneSupimposer::new();
        sup.set(&reference, &mobile);
        sup.run();

        let (rot, tran) = sup.get_rotran();
        let tf = sup.get_transformed();

        // Every transformed Cβ should be close to its reference Cβ
        let n_res = mobile.len() / 2;
        for i in 0..n_res {
            let cb_tf  = tf[2*i+1];
            let cb_ref = reference[2*i+1];
            let d: f32 = (0..3).map(|k| (cb_tf[k]-cb_ref[k]).powi(2)).sum::<f32>().sqrt();
            assert!(d < 0.1, "residue {i}: Cβ dist {d:.4} too large");
        }
        let _ = (rot, tran); // suppress unused warning
    }

    #[test]
    fn test_init_rms_gt_rms() {
        let (reference, mobile) = make_test_coords(10);
        let mut sup = PointToPlaneSupimposer::new();
        sup.set(&reference, &mobile);
        sup.run();

        let init = sup.get_init_rms();
        let fin  = sup.get_rms();
        println!("init_rms={init:.4}, rms={fin:.4}");
        assert!(init >= fin || fin < 1e-5,
            "init_rms should be ≥ rms after optimisation");
    }

    /// Mirror of the kabsch `test_kabsch_superimposer` test.
    /// Each source/target Coordinate is treated as a Cα; a synthetic Cβ is
    /// created by adding `cb_offset` so we get a valid CA/CB pair per residue.
    #[test]
    fn test_real_structure_coords() {
        // --- raw CA coordinates (same as kabsch test) ----------------------
        let source_ca: &[[f32; 3]] = &[
            [6.994,  8.354, 42.405],
            [9.429,  7.479, 48.266],
            [5.547,  0.158, 42.050],
        ];
        let target1_ca: &[[f32; 3]] = &[
            [-13.958, -1.741, -4.223],
            [-12.833,  3.134, -7.780],
            [ -5.720, -2.218, -3.368],
        ];
        let target2_ca: &[[f32; 3]] = &[
            [-4.924,  5.813, -9.485],
            [-0.499, 10.073, -8.059],
            [-0.792,  0.658, -4.430],
        ];

        // Synthetic Cβ = Cα + cb_offset (constant vector, same for all atoms).
        // Because the offset is rigid it behaves like a real CA→CB bond,
        // and the optimal rotation that aligns CAs also aligns CBs.
        let cb_offset: [f32; 3] = [1.0, 0.0, 0.0];
        let interleave = |cas: &[[f32; 3]]| -> Vec<[f32; 3]> {
            cas.iter()
                .flat_map(|&ca| {
                    let cb = [ca[0] + cb_offset[0], ca[1] + cb_offset[1], ca[2] + cb_offset[2]];
                    [ca, cb]
                })
                .collect()
        };

        let source_coords  = interleave(source_ca);
        let target1_coords = interleave(target1_ca);
        let target2_coords = interleave(target2_ca);

        // --- Kabsch RMSD (Cα only) for reference -------------------------
        let ca_to_f32 = |cas: &[[f32;3]]| -> Vec<[f32;3]> { cas.to_vec() };
        let mut kab = KabschSuperimposer::new();

        kab.set(&ca_to_f32(target1_ca), &ca_to_f32(source_ca));
        kab.run();
        let kab_rms1 = kab.get_rms();

        kab.set(&ca_to_f32(target2_ca), &ca_to_f32(source_ca));
        kab.run();
        let kab_rms2 = kab.get_rms();

        // --- superpose source (moving) onto target1 (fixed) ----------------
        let mut sup = PointToPlaneSupimposer::new();
        sup.set(&target1_coords, &source_coords);
        sup.run();
        let init1 = sup.get_init_rms();
        let ptp_rms1 = sup.get_rms();
        let (rot1, tran1) = sup.get_rotran();
        let tf1 = sup.get_transformed();
        println!("source → target1");
        println!("  Kabsch RMSD (Cα):        {kab_rms1:.6}");
        println!("  PtP init_rms:            {init1:.6}");
        println!("  PtP rms (after opt):     {ptp_rms1:.6}");
        println!("  rot:  {:?}", rot1);
        println!("  tran: {:?}", tran1);
        println!("  transformed[0]: {:?}", tf1[0]);

        // --- superpose source (moving) onto target2 (fixed) ----------------
        sup.set(&target2_coords, &source_coords);
        sup.run();
        let init2 = sup.get_init_rms();
        let ptp_rms2 = sup.get_rms();
        let (rot2, tran2) = sup.get_rotran();
        let tf2 = sup.get_transformed();
        println!("source → target2");
        println!("  Kabsch RMSD (Cα):        {kab_rms2:.6}");
        println!("  PtP init_rms:            {init2:.6}");
        println!("  PtP rms (after opt):     {ptp_rms2:.6}");
        println!("  rot:  {:?}", rot2);
        println!("  tran: {:?}", tran2);
        println!("  transformed[0]: {:?}", tf2[0]);

        // The point-to-plane score should be small for near-identical structures.
        assert!(ptp_rms1 < 1.0, "source→target1 ptp_rms={ptp_rms1:.6} unexpectedly large");
        assert!(ptp_rms2 < 1.0, "source→target2 ptp_rms={ptp_rms2:.6} unexpectedly large");
        // Optimisation must not increase the score.
        assert!(ptp_rms1 <= init1 + 1e-4, "ptp_rms1 > init1 (optimisation diverged)");
        assert!(ptp_rms2 <= init2 + 1e-4, "ptp_rms2 > init2 (optimisation diverged)");
    }

    #[test]
    fn test_point_to_plane_distance_degenerate() {
        // Anti-parallel CBs → norm = 0 → should return 0, not NaN
        let ca = [0.0_f32; 3];
        let cb_ref    = [1.0, 0.0, 0.0_f32];
        let ca_target = [0.0_f32; 3];
        let cb_target = [-1.0, 0.0, 0.0_f32];
        let v = point_to_plane_distance(ca, cb_ref, ca_target, cb_target);
        assert!(!v.is_nan(), "should not be NaN");
        assert_eq!(v, 0.0);
    }
}
