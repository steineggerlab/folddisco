// File: lms_qcp.rs
// Created: 2025-08-16
// Author: Hyunbin Kim (khb7840@gmail.com)
// Description:
//   Reference: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-29
//   - Find subset of residues with minimal RMSD instead of all matching residues
//   - Tunable quantiles: q_seed (default 0.5), q_report (default 0.75).

use crate::structure::coordinate::Coordinate;

#[derive(Debug, Clone, Copy)]
pub struct LmsQcpParams {
    /// Å cutoff to stop forward growth once min_core is reached.
    pub r_max: f32,
    /// Minimal core size before allowing stop (Default: N/2 set in run()).
    pub min_core: Option<usize>,
    /// Minimal subset for seeding (3 for rigid 3D).
    pub k: usize,
    /// Random seed trials for seeding.
    pub t_init: usize,
    /// RNG seed.
    pub seed: u64,
    /// Quantile for seed scoring (0..1). Default 0.5 (median).
    pub q_seed: f32,
    /// Quantile reported over ALL pairs under final transform (0..1). Default 0.75.
    pub q_report: f32,
}

impl Default for LmsQcpParams {
    fn default() -> Self {
        Self {
            r_max: 2.0,
            min_core: None,
            k: 3,
            t_init: 500,
            seed: 0xC0FF_EE00_5EEDu64,
            q_seed: 0.5,
            q_report: 0.75,
        }
    }
}

#[derive(Debug)]
pub struct LmsQcpSuperimposer {
    pub reference_coords: Option<Vec<[f32; 3]>>,
    pub coords: Option<Vec<[f32; 3]>>,
    pub rot: Option<[[f32; 3]; 3]>,
    pub tran: Option<[f32; 3]>,
    pub transformed_coords: Option<Vec<[f32; 3]>>,
    pub rms_inliers: Option<f32>,            // RMSD over final core
    pub q_report_residual_all: Option<f32>,  // q_report distance over ALL pairs
    pub core_indices: Vec<usize>,
    pub natoms: usize,
    pub params: LmsQcpParams,
}

impl LmsQcpSuperimposer {
    pub fn new() -> Self { Self::with_params(LmsQcpParams::default()) }
    pub fn with_params(params: LmsQcpParams) -> Self {
        Self {
            reference_coords: None,
            coords: None,
            rot: None,
            tran: None,
            transformed_coords: None,
            rms_inliers: None,
            q_report_residual_all: None,
            core_indices: Vec::new(),
            natoms: 0,
            params,
        }
    }

    pub fn set_atoms(&mut self, fixed: &[Coordinate], moving: &[Coordinate]) {
        assert_eq!(fixed.len(), moving.len(), "Fixed and moving atom lists differ");
        let fixed_coords: Vec<[f32; 3]> = fixed.iter().map(|c| c.to_array()).collect();
        let moving_coords: Vec<[f32; 3]> = moving.iter().map(|c| c.to_array()).collect();
        self.set(&fixed_coords, &moving_coords);
    }

    pub fn set(&mut self, reference_coords: &[[f32; 3]], coords: &[[f32; 3]]) {
        assert_eq!(reference_coords.len(), coords.len(), "Lengths differ");
        assert!(coords.len() >= 3, "Need ≥3 pairs");
        self.reference_coords = Some(reference_coords.to_vec());
        self.coords = Some(coords.to_vec());
        self.natoms = coords.len();
        self.rot = None; self.tran = None; self.transformed_coords = None;
        self.rms_inliers = None; self.q_report_residual_all = None; self.core_indices.clear();
    }

    pub fn run(&mut self) {
        let n = self.natoms;
        let mut rng = SmallRng::new(self.params.seed);
        let k = self.params.k.max(3).min(n);
        let t_init = self.params.t_init.max(1);
        let q_seed = self.params.q_seed.clamp(0.0, 1.0);

        // ---------- 1) Seed by minimizing q_seed quantile over residuals ----------
        // We implement a fast path for k=3 (recommended).
        let mut best_seed: [usize; 3] = [0, 1, 2];
        let mut best_qval = f32::INFINITY;

        {
            let movs = self.coords.as_ref().unwrap();
            let refs = self.reference_coords.as_ref().unwrap();
            let mut res_sq: Vec<f32> = Vec::with_capacity(n.saturating_sub(3));

            for _ in 0..t_init {
                let s = if k == 3 {
                    sample_three_non_collinear(movs, &mut rng, 64)
                } else {
                    // fallback: sample 3 anyway; growth will correct
                    sample_three_non_collinear(movs, &mut rng, 64)
                };
                if s.is_none() { continue; }
                let seed = s.unwrap();

                let mut st = RunningStats::new();
                st.add(seed[0], movs, refs);
                st.add(seed[1], movs, refs);
                st.add(seed[2], movs, refs);

                let (r, t) = qcp_from_stats(&st);

                res_sq.clear();
                for i in 0..n {
                    if i == seed[0] || i == seed[1] || i == seed[2] { continue; }
                    let p = apply(r, t, movs[i]);
                    res_sq.push(dist2(p, refs[i]));
                }
                let qv = select_quantile_squared(&mut res_sq, q_seed);
                if qv < best_qval {
                    best_qval = qv;
                    best_seed = seed;
                }
            }
        } // short-lived borrow ends here

        // ---------- 2) Forward-search LMS: grow core until r_max (once min_core reached) ----------
        let min_core = self.params.min_core.unwrap_or(n / 2).max(3);
        let r2_max = self.params.r_max * self.params.r_max;

        // init core + running stats
        let mut st = RunningStats::new();
        let mut in_core = vec![false; n];
        for &i in &best_seed {
            let movs = self.coords.as_ref().unwrap();
            let refs = self.reference_coords.as_ref().unwrap();
            st.add(i, movs, refs);
            in_core[i] = true;
        }
        self.core_indices = best_seed.to_vec();

        loop {
            let (r, t) = qcp_from_stats(&st);

            // find best candidate using short-lived borrows
            let (best_i, best_r2) = {
                let movs = self.coords.as_ref().unwrap();
                let refs = self.reference_coords.as_ref().unwrap();
                let mut best_i = None;
                let mut best_r2 = f32::INFINITY;
                for i in 0..n {
                    if in_core[i] { continue; }
                    let p = apply(r, t, movs[i]);
                    let d2 = dist2(p, refs[i]);
                    if d2 < best_r2 { best_r2 = d2; best_i = Some(i); }
                }
                (best_i, best_r2)
            }; // borrows end before we mutate self

            if let Some(i_add) = best_i {
                let allow_stop = self.core_indices.len() >= min_core;
                if allow_stop && best_r2 > r2_max {
                    self.finish(r, t);
                    return;
                }
                // update stats and core
                {
                    let movs = self.coords.as_ref().unwrap();
                    let refs = self.reference_coords.as_ref().unwrap();
                    st.add(i_add, movs, refs);
                }
                in_core[i_add] = true;
                self.core_indices.push(i_add);
                if self.core_indices.len() == n {
                    self.finish(r, t);
                    return;
                }
            } else {
                self.finish(r, t);
                return;
            }
        }
    }

    #[inline]
    pub fn get_rotran(&self) -> ([[f32; 3]; 3], [f32; 3]) {
        (self.rot.expect("run() not called"), self.tran.expect("run() not called"))
    }

    pub fn get_transformed(&mut self) -> Vec<[f32; 3]> {
        if self.transformed_coords.is_none() {
            let movs = self.coords.as_ref().unwrap();
            let (r, t) = self.get_rotran();
            self.transformed_coords = Some(movs.iter().copied().map(|p| apply(r, t, p)).collect());
        }
        self.transformed_coords.clone().unwrap()
    }

    #[inline]
    pub fn get_rms_inliers(&self) -> f32 { self.rms_inliers.expect("run() not called") }

    #[inline]
    pub fn get_q_report_all(&self) -> f32 { self.q_report_residual_all.expect("run() not called") }

    #[inline]
    pub fn inliers(&self) -> &[usize] { &self.core_indices }

    #[inline(always)]
    pub fn core_percent(&self) -> f32 {
        (self.core_indices.len() as f32) / (self.natoms as f32)
    }

    // finalize metrics & caches; reads coords/refs internally (no external borrows)
    fn finish(&mut self, r: [[f32; 3]; 3], t: [f32; 3]) {
        self.rot = Some(r);
        self.tran = Some(t);

        let movs = self.coords.as_ref().unwrap();
        let refs = self.reference_coords.as_ref().unwrap();

        // RMS over current core
        let mut sum = 0.0f32;
        for &i in &self.core_indices {
            let p = apply(r, t, movs[i]);
            sum += dist2(p, refs[i]);
        }
        self.rms_inliers = Some((sum / (self.core_indices.len() as f32)).sqrt());

        // q_report over ALL pairs (select quantile of squared residuals, then sqrt)
        let q = self.params.q_report.clamp(0.0, 1.0);
        let mut res_sq: Vec<f32> = (0..movs.len()).map(|i| dist2(apply(r, t, movs[i]), refs[i])).collect();

        let qsq = select_quantile_squared(&mut res_sq, q);
        self.q_report_residual_all = Some(qsq.sqrt());

        // cache transformed
        self.transformed_coords = Some(movs.iter().copied().map(|p| apply(r, t, p)).collect());
    }
}

// ====== Running stats for incremental QCP ======

#[derive(Clone, Debug)]
struct RunningStats {
    n: usize,
    sum_x: [f64; 3],
    sum_y: [f64; 3],
    // squared norms sums
    sxx: f64, // Σ ||x||^2
    syy: f64, // Σ ||y||^2
    // cross outer sum: Σ y x^T (row-major)
    syx: [[f64; 3]; 3],
}

impl RunningStats {
    #[inline(always)]
    fn new() -> Self {
        Self {
            n: 0,
            sum_x: [0.0; 3],
            sum_y: [0.0; 3],
            sxx: 0.0,
            syy: 0.0,
            syx: [[0.0; 3]; 3],
        }
    }
    #[inline(always)]
    fn add(&mut self, i: usize, movs: &[[f32; 3]], refs: &[[f32; 3]]) {
        // Convert to f64 for internal calculations
        let x = [movs[i][0] as f64, movs[i][1] as f64, movs[i][2] as f64];
        let y = [refs[i][0] as f64, refs[i][1] as f64, refs[i][2] as f64];
        self.n += 1;
        self.sum_x[0] += x[0]; self.sum_x[1] += x[1]; self.sum_x[2] += x[2];
        self.sum_y[0] += y[0]; self.sum_y[1] += y[1]; self.sum_y[2] += y[2];
        self.sxx += x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        self.syy += y[0]*y[0] + y[1]*y[1] + y[2]*y[2];
        // syx += y x^T
        self.syx[0][0] += y[0]*x[0]; self.syx[0][1] += y[0]*x[1]; self.syx[0][2] += y[0]*x[2];
        self.syx[1][0] += y[1]*x[0]; self.syx[1][1] += y[1]*x[1]; self.syx[1][2] += y[1]*x[2];
        self.syx[2][0] += y[2]*x[0]; self.syx[2][1] += y[2]*x[1]; self.syx[2][2] += y[2]*x[2];
    }
    #[inline(always)]
    fn means(&self) -> ([f64; 3], [f64; 3]) {
        let inv = 1.0f64 / self.n as f64;
        ([self.sum_x[0]*inv, self.sum_x[1]*inv, self.sum_x[2]*inv],
         [self.sum_y[0]*inv, self.sum_y[1]*inv, self.sum_y[2]*inv])
    }
}

// QCP from incremental stats:
// - A = Σ (y-μy)(x-μx)^T = syx - n * μy μx^T
// - e0 = (trace(Σ (y-μy)(y-μy)^T) + trace(Σ (x-μx)(x-μx)^T)) / 2
//      = ( (syy - n||μy||^2) + (sxx - n||μx||^2) ) / 2
#[inline]
fn qcp_from_stats(st: &RunningStats) -> ([[f32; 3]; 3], [f32; 3]) {
    debug_assert!(st.n >= 3);
    let (mux, muy) = st.means();
    // A = syx - n * μy μx^T (all in f64)
    let mut a = [[0.0f64; 3]; 3];
    a[0][0] = st.syx[0][0] - st.n as f64 * (muy[0]*mux[0]);
    a[0][1] = st.syx[0][1] - st.n as f64 * (muy[0]*mux[1]);
    a[0][2] = st.syx[0][2] - st.n as f64 * (muy[0]*mux[2]);
    a[1][0] = st.syx[1][0] - st.n as f64 * (muy[1]*mux[0]);
    a[1][1] = st.syx[1][1] - st.n as f64 * (muy[1]*mux[1]);
    a[1][2] = st.syx[1][2] - st.n as f64 * (muy[1]*mux[2]);
    a[2][0] = st.syx[2][0] - st.n as f64 * (muy[2]*mux[0]);
    a[2][1] = st.syx[2][1] - st.n as f64 * (muy[2]*mux[1]);
    a[2][2] = st.syx[2][2] - st.n as f64 * (muy[2]*mux[2]);

    let mu2x = mux[0]*mux[0] + mux[1]*mux[1] + mux[2]*mux[2];
    let mu2y = muy[0]*muy[0] + muy[1]*muy[1] + muy[2]*muy[2];
    let e0 = 0.5f64 * ((st.syy - st.n as f64 * mu2y) + (st.sxx - st.n as f64 * mu2x)).max(0.0);

    let rot = qcp_from_a_e0(a, e0, st.n);

    // translation: t = μy - R μx
    let rx = apply_rot_f64(rot, mux);
    let t = [muy[0] - rx[0], muy[1] - rx[1], muy[2] - rx[2]];
    
    // Convert back to f32 for output
    let rot_f32 = [
        [rot[0][0] as f32, rot[0][1] as f32, rot[0][2] as f32],
        [rot[1][0] as f32, rot[1][1] as f32, rot[1][2] as f32],
        [rot[2][0] as f32, rot[2][1] as f32, rot[2][2] as f32],
    ];
    let t_f32 = [t[0] as f32, t[1] as f32, t[2] as f32];
    (rot_f32, t_f32)
}

// ====== QCP kernel (A,e0 form) ======

#[inline]
fn qcp_from_a_e0(a: [[f64; 3]; 3], e0: f64, _n: usize) -> [[f64; 3]; 3] {
    let sxx = a[0][0]; let sxy = a[0][1]; let sxz = a[0][2];
    let syx = a[1][0]; let syy = a[1][1]; let syz = a[1][2];
    let szx = a[2][0]; let szy = a[2][1]; let szz = a[2][2];

    let sxx2 = sxx*sxx; let syy2 = syy*syy; let szz2 = szz*szz;
    let sxy2 = sxy*sxy; let syz2 = syz*syz; let sxz2 = sxz*sxz;
    let syx2 = syx*syx; let szy2 = szy*szy; let szx2 = szx*szx;

    let syz_szy_m_syy_szz2 = 2.0*(syz*szy - syy*szz);
    let sxx2_syy2_szz2_syz2_szy2 = syy2 + szz2 - sxx2 + syz2 + szy2;

    let c2 = -2.0*(sxx2 + syy2 + szz2 + sxy2 + syx2 + sxz2 + szx2 + syz2 + szy2);
    let c1 = 8.0*(sxx*syz*szy + syy*szx*sxz + szz*sxy*syx - sxx*syy*szz - syz*szx*sxy - szy*syx*sxz);

    let sxz_p_szx = sxz + szx; let syz_p_szy = syz + szy; let sxy_p_syx = sxy + syx;
    let syz_m_szy = syz - szy; let sxz_m_szx = sxz - szx; let sxy_m_syx = sxy - syx;
    let sxx_p_syy = sxx + syy; let sxx_m_syy = sxx - syy;
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

    // Newton on quartic for max eigenvalue
    let mut lam = e0.max(0.0);
    let eps = 1e-15f64;
    for _ in 0..10 {
        let x2 = lam*lam;
        let b = (x2 + c2)*lam;
        let a = b + c1;
        let f = a*lam + c0;
        let fp = 2.0*x2*lam + b + a;
        let delta = f / (fp + eps);
        let nlam = (lam - delta).abs();
        if (nlam - lam).abs() < eps * nlam { lam = nlam; break; }
        lam = nlam;
    }

    // Build quaternion (implicit 4×4) and normalize
    let a11 = sxx_p_syy + szz - lam;
    let a12 = syz_m_szy;
    let a13 = neg_sxz_m_szx;
    let a14 = sxy_m_syx;
    let a21 = a12;
    let a22 = sxx_m_syy - szz - lam;
    let a23 = sxy_p_syx;
    let a24 = sxz_p_szx;
    let a31 = a13;
    let a32 = a23;
    let a33 = syy - sxx - szz - lam;
    let a34 = syz_p_szy;
    let a41 = a14;
    let a42 = a24;
    let a43 = a34;
    let a44 = szz - sxx_p_syy - lam;

    let a3344_4334 = a33*a44 - a43*a34;
    let a3244_4234 = a32*a44 - a42*a34;
    let a3243_4233 = a32*a43 - a42*a33;
    let a3143_4133 = a31*a43 - a41*a33;
    let a3144_4134 = a31*a44 - a41*a34;
    let a3142_4132 = a31*a42 - a41*a32;

    let mut q1 = a22*a3344_4334 - a23*a3244_4234 + a24*a3243_4233;
    let mut q2 = -a21*a3344_4334 + a23*a3144_4134 - a24*a3143_4133;
    let mut q3 = a21*a3244_4234 - a22*a3144_4134 + a24*a3142_4132;
    let mut q4 = -a21*a3243_4233 + a22*a3143_4133 - a23*a3142_4132;

    let mut qsqr = q1*q1 + q2*q2 + q3*q3 + q4*q4;
    let evec_prec = 1e-12f64;

    if qsqr < evec_prec {
        // fallback 1
        q1 = a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
        q3 = a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
        q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
        qsqr = q1*q1 + q2*q2 + q3*q3 + q4*q4;

        if qsqr < evec_prec {
            // fallback 2 -> identity
            return [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
        }
    }

    let inv = qsqr.sqrt().recip();
    q1*=inv; q2*=inv; q3*=inv; q4*=inv;

    // quaternion -> rotation
    let a2 = q1*q1; let x2 = q2*q2; let y2 = q3*q3; let z2 = q4*q4;
    let xy = q2*q3; let az = q1*q4; let zx = q4*q2; let ay = q1*q3; let yz = q3*q4; let ax = q1*q2;

    let mut rot = [[0.0f64;3];3];
    rot[0][0] = a2 + x2 - y2 - z2;
    rot[0][1] = 2.0*(xy + az);
    rot[0][2] = 2.0*(zx - ay);
    rot[1][0] = 2.0*(xy - az);
    rot[1][1] = a2 - x2 + y2 - z2;
    rot[1][2] = 2.0*(yz + ax);
    rot[2][0] = 2.0*(zx + ay);
    rot[2][1] = 2.0*(yz - ax);
    rot[2][2] = a2 - x2 - y2 + z2;
    rot
}

// ====== Fast quantile selection over squared distances ======

#[inline]
fn select_quantile_squared(res_sq: &mut [f32], q: f32) -> f32 {
    if res_sq.is_empty() { return 0.0; }
    if res_sq.len() == 1 { return res_sq[0]; }
    let q = q.clamp(0.0, 1.0);
    let n = res_sq.len();
    let pos = (q * ((n - 1) as f32)).round() as usize;
    // let (_, nth, _) = res_sq.select_nth_unstable_by(pos, |a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let (_, nth, _) = res_sq.select_nth_unstable_by(pos, |a, b| a.partial_cmp(b).unwrap());

    *nth
}

// ====== Utils (inlined & branch-light) ======

#[inline(always)]
fn apply_rot(r: [[f32; 3]; 3], v: [f32; 3]) -> [f32; 3] {
    [
        r[0][0]*v[0] + r[0][1]*v[1] + r[0][2]*v[2],
        r[1][0]*v[0] + r[1][1]*v[1] + r[1][2]*v[2],
        r[2][0]*v[0] + r[2][1]*v[1] + r[2][2]*v[2],
    ]
}
#[inline(always)]
fn apply_rot_f64(r: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        r[0][0]*v[0] + r[0][1]*v[1] + r[0][2]*v[2],
        r[1][0]*v[0] + r[1][1]*v[1] + r[1][2]*v[2],
        r[2][0]*v[0] + r[2][1]*v[1] + r[2][2]*v[2],
    ]
}
#[inline(always)]
fn apply(r: [[f32; 3]; 3], t: [f32; 3], v: [f32; 3]) -> [f32; 3] {
    let w = apply_rot(r, v);
    [w[0] + t[0], w[1] + t[1], w[2] + t[2]]
}
#[inline(always)]
fn dist2(a: [f32; 3], b: [f32; 3]) -> f32 {
    let dx = a[0]-b[0]; let dy = a[1]-b[1]; let dz = a[2]-b[2];
    dx*dx + dy*dy + dz*dz
}

// Non-collinear 3-sample (rejects nearly collinear triplets)
#[inline]
fn sample_three_non_collinear(movs: &[[f32; 3]], rng: &mut SmallRng, max_tries: usize) -> Option<[usize; 3]> {
    let n = movs.len();
    if n < 3 { return None; }
    for _ in 0..max_tries {
        let i = rng.gen_range(n);
        let mut j = rng.gen_range(n); if j == i { j = (j + 1) % n; }
        let mut k = rng.gen_range(n); while k == i || k == j { k = (k + 1) % n; }
        let v1 = [movs[j][0]-movs[i][0], movs[j][1]-movs[i][1], movs[j][2]-movs[i][2]];
        let v2 = [movs[k][0]-movs[i][0], movs[k][1]-movs[i][1], movs[k][2]-movs[i][2]];
        let cx = v1[1]*v2[2] - v1[2]*v2[1];
        let cy = v1[2]*v2[0] - v1[0]*v2[2];
        let cz = v1[0]*v2[1] - v1[1]*v2[0];
        let area2 = cx*cx + cy*cy + cz*cz;
        if area2 > 1e-6 { return Some([i, j, k]); }
    }
    None
}

// Tiny RNG (xorshift64*)
struct SmallRng { state: u64 }
impl SmallRng {
    #[inline(always)]
    fn new(seed: u64) -> Self {
        // SplitMix-ish seed scramble
        let mut x = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        x = (x ^ (x>>30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        x = (x ^ (x>>27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        Self { state: x ^ (x>>31) }
    }
    #[inline(always)]
    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        self.state = x; x
    }
    #[inline(always)]
    fn gen_range(&mut self, end: usize) -> usize {
        (self.next_u64() as usize) % end
    }
}

// ====== tests ======

#[cfg(test)]
mod tests {
    use super::*;

    // TODO: Make this test runnable
    #[test]
    #[ignore]
    fn lms_qcp_excludes_outlier() {
        // 3 coherent + 1 outlier
        let fixed = vec![
            Coordinate::new(0.0, 0.0, 0.0),
            Coordinate::new(1.0, 0.0, 0.0),
            Coordinate::new(0.0, 1.0, 0.0),
            Coordinate::new(50.0, -20.0, 7.0), // outlier
        ];
        let moving = vec![
            Coordinate::new(1.0, 0.0, 0.0),
            Coordinate::new(2.0, 0.0, 0.0),
            Coordinate::new(1.0, 1.0, 0.0),
            Coordinate::new(-200.0, 111.0, -50.0), // mismatched
        ];

        let mut sup = LmsQcpSuperimposer::with_params(LmsQcpParams {
            r_max: 2.0, min_core: Some(3), k: 3, t_init: 400,
            seed: 123, q_seed: 0.5, q_report: 0.75
        });
        sup.set_atoms(&fixed, &moving);
        sup.run();

        assert!(sup.inliers().len() >= 3);
        assert!(sup.get_rms_inliers() < 1e-4);
        // Also ensure we can get rotran
        let _ = sup.get_rotran();
    }
}
