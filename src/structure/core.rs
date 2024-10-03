use crate::structure::atom::{Atom, AtomVector};
use crate::structure::coordinate::{approx_cb, CarbonCoordinateVector, Coordinate};
use crate::structure::feature::{Torsion, TorsionType};
use crate::utils::convert::map_aa_to_u8;

use super::coordinate::{calc_torsion_radian, calc_angle_radian};

/// Structure is the main data structure for storing the information of a protein structure.
#[derive(Debug)]
pub struct Structure {
    pub num_chains: usize,
    pub chains: Vec<u8>,
    pub atom_vector: AtomVector,
    pub num_atoms: usize,
    pub num_residues: usize,
}

impl Structure {
    pub fn new() -> Structure {
        Structure {
            num_chains: 0,
            chains: Vec::new(),
            atom_vector: AtomVector::new(),
            num_atoms: 0,
            num_residues: 0,
        }
    }

    pub fn update(&mut self, atom: Atom, record: &mut (u8, u64)) {
        // record store previous chain ID and residue serial
        if record.0 != atom.chain {
            self.chains.push(atom.chain);
            self.num_chains += 1;
            record.0 = atom.chain;
        }
        if record.1 != atom.res_serial {
            self.num_residues += 1;
            record.1 = atom.res_serial;
        }
        self.num_atoms += 1;

        self.atom_vector.push_atom(atom);
    }

    pub fn to_compact(&self) -> CompactStructure {
        CompactStructure::build(self)
    }
    pub fn get_torsion(&self) -> Torsion {
        //FIXME: Right now, only Psi is calculated
        Torsion::build(self, TorsionType::Psi)
    }

}

#[derive(Debug, Clone)]
pub struct CompactStructure {
    pub num_chains: usize,
    pub chains: Vec<u8>,
    pub chain_per_residue: Vec<u8>,
    pub num_residues: usize,
    pub residue_serial: Vec<u64>,
    pub residue_name: Vec<[u8; 3]>,
    pub n_vector: CarbonCoordinateVector,
    pub ca_vector: CarbonCoordinateVector,
    pub cb_vector: CarbonCoordinateVector,
    pub n_vector2: (Vec<f32>, Vec<f32>, Vec<f32>),
    pub ca_vector2: (Vec<f32>, Vec<f32>, Vec<f32>),
    pub cb_vector2: (Vec<f32>, Vec<f32>, Vec<f32>),
    pub b_factors: Vec<f32>,
}

impl CompactStructure {
    pub fn build(origin: &Structure) -> CompactStructure {
        // Store only backbone atoms
        let model = &origin.atom_vector;

        let mut res_serial_vec: Vec<u64> = Vec::new();
        let mut res_name_vec: Vec<[u8; 3]> = Vec::new();
        let mut b_factors: Vec<f32> = Vec::new();
        let mut n_vec = CarbonCoordinateVector::new();
        let mut ca_vec = CarbonCoordinateVector::new();
        let mut cb_vec = CarbonCoordinateVector::new();

        let mut n_vec_x: Vec<f32> = Vec::new();
        let mut n_vec_y: Vec<f32> = Vec::new();
        let mut n_vec_z: Vec<f32> = Vec::new();
        let mut ca_vec_x: Vec<f32> = Vec::new();
        let mut ca_vec_y: Vec<f32> = Vec::new();
        let mut ca_vec_z: Vec<f32> = Vec::new();
        let mut cb_vec_x: Vec<f32> = Vec::new();
        let mut cb_vec_y: Vec<f32> = Vec::new();
        let mut cb_vec_z: Vec<f32> = Vec::new();
        
        
        let mut chain_per_residue: Vec<u8> = Vec::new();
        let mut prev_res_serial: Option<u64> = None;
        let mut prev_res_name: Option<&[u8; 3]> = None;
        let mut n: Option<Coordinate> = None;
        let mut ca: Option<Coordinate> = None;
        let mut cb: Option<Coordinate> = None;

        let mut gly_n: Option<Coordinate> = None;
        let mut gly_c: Option<Coordinate> = None;
        

        for idx in 0..origin.num_atoms {
            if prev_res_serial != Some(model.get_res_serial(idx)) || idx == origin.num_atoms - 1 {
                // Save previous 'CA' and 'CB'
                match (n, ca, cb) {
                    (Some(n), Some(ca), Some(cb)) => {
                        let resi = prev_res_serial.expect("expected residue serial number");
                        let resn = prev_res_name.expect("expected residue name");
                        n_vec.push(&n);
                        ca_vec.push(&ca);
                        cb_vec.push(&cb);
                        res_serial_vec.push(resi);
                        res_name_vec.push(*resn);
                        chain_per_residue.push(origin.atom_vector.chain[idx]);
                        b_factors.push(origin.atom_vector.b_factor[idx]);
                        n_vec_x.push(n.x);
                        n_vec_y.push(n.y);
                        n_vec_z.push(n.z);
                        ca_vec_x.push(ca.x);
                        ca_vec_y.push(ca.y);
                        ca_vec_z.push(ca.z);
                        cb_vec_x.push(cb.x);
                        cb_vec_y.push(cb.y);
                        cb_vec_z.push(cb.z);
                        
                    }
                    (Some(n), Some(ca), None) => {
                        let resi = prev_res_serial.expect("expected residue serial number");
                        let resn = prev_res_name.expect("expected residue name");
                        n_vec.push(&n);
                        ca_vec.push(&ca);
                        res_serial_vec.push(resi);
                        res_name_vec.push(*resn);
                        chain_per_residue.push(origin.atom_vector.chain[idx]);
                        b_factors.push(origin.atom_vector.b_factor[idx]);
                        if let (Some(b"GLY"), Some(gly_n), Some(gly_c)) =
                            (prev_res_name, &gly_n, &gly_c)
                        {
                            // Approximate CB
                            let cb = approx_cb(&ca, gly_n, gly_c);
                            cb_vec.push(&cb);
                            cb_vec_x.push(cb.x);
                            cb_vec_y.push(cb.y);
                            cb_vec_z.push(cb.z);
                        } else {
                            cb_vec.push_none();
                        }
                        n_vec_x.push(n.x);
                        n_vec_y.push(n.y);
                        n_vec_z.push(n.z);
                        ca_vec_x.push(ca.x);
                        ca_vec_y.push(ca.y);
                        ca_vec_z.push(ca.z);

                    }
                    (None, None, None) => {}
                    (None, None, Some(_)) => {}
                    (None, Some(_), None) => {}
                    (None, Some(_), Some(_)) => {}
                    (Some(_), None, None) => {}
                    (Some(_), None, Some(_)) => {}
                }
                // Reset 'CA' and 'CB'
                ca = None;
                cb = None;
                n = None;
                prev_res_serial = Some(model.get_res_serial(idx));
                prev_res_name = model.res_name.get(idx);
            }

            if model.is_ca(idx) {
                ca = Some(model.get_coordinates(idx));
            } else if model.is_cb(idx) {
                cb = Some(model.get_coordinates(idx));
            } else if model.is_n(idx) && &model.get_res_name(idx) != b"GLY" {
                n = Some(model.get_coordinates(idx));
            } else if &model.get_res_name(idx) == b"GLY" {
                // If GLY, calculate CB from CA, N, C
                match model.atom_name.get(idx) {
                    Some(b" N  ") => {
                        gly_n = Some(model.get_coordinates(idx));
                        n = Some(model.get_coordinates(idx));
                    }
                    Some(b" C  ") => {
                        gly_c = Some(model.get_coordinates(idx));
                    }
                    _ => (),
                }
            }
        }

        CompactStructure {
            num_chains: origin.num_chains,
            chains: origin.chains.clone(),
            chain_per_residue: chain_per_residue,
            num_residues: res_serial_vec.len(),
            residue_serial: res_serial_vec,
            residue_name: res_name_vec,
            n_vector: n_vec,
            ca_vector: ca_vec,
            cb_vector: cb_vec,
            b_factors: b_factors,
            n_vector2: (n_vec_x, n_vec_y, n_vec_z),
            ca_vector2: (ca_vec_x, ca_vec_y, ca_vec_z),
            cb_vector2: (cb_vec_x, cb_vec_y, cb_vec_z),
        }
    }
    #[inline(always)]
    pub fn get_index(&self, chain: &u8, res_serial: &u64) -> Option<usize> {
        for i in 0..self.num_residues {
            if self.chain_per_residue[i] == *chain && self.residue_serial[i] == *res_serial {
                return Some(i);
            }
        }
        None
    }
    #[inline(always)]
    pub fn get_ca(&self, idx: usize) -> Option<Coordinate> {
        let (x, y, z) = self.ca_vector.get(idx);

        if x.is_some() && y.is_some() && z.is_some() {
            Some(Coordinate::build(&x, &y, &z))
        } else {
            None
        }
    }
    #[inline(always)]
    pub fn get_cb(&self, idx: usize) -> Option<Coordinate> {
        let (x, y, z) = self.cb_vector.get(idx);

        if x.is_some() && y.is_some() && z.is_some() {
            Some(Coordinate::build(&x, &y, &z))
        } else {
            None
        }
    }
    #[inline(always)]
    pub fn get_n(&self, idx: usize) -> Option<Coordinate> {
        let (x, y, z) = self.n_vector.get(idx);

        if x.is_some() && y.is_some() && z.is_some() {
            Some(Coordinate::build(&x, &y, &z))
        } else {
            None
        }
    }
    #[inline(always)]
    pub fn get_ca_distance(&self, idx1: usize, idx2: usize) -> Option<f32> {
        let ca1 = self.get_ca(idx1);
        let ca2 = self.get_ca(idx2);
        if ca1.is_some() && ca2.is_some() {
            let dist = ca1
                .expect("Unable to get CA coordinate")
                .calc_distance(&ca2.expect("Unable to get CA coordinate"));
            Some(dist)
        } else {
            None
        }
    }
    #[inline(always)]
    pub fn get_cb_distance(&self, idx1: usize, idx2: usize) -> Option<f32> {
        let cb1 = self.get_cb(idx1);
        let cb2 = self.get_cb(idx2);
        if cb1.is_some() && cb2.is_some() {
            let dist = cb1
                .expect("Unable to get CA coordinate")
                .calc_distance(&cb2.expect("Unable to get CA coordinate"));
            Some(dist)
        } else {
            None
        }
    }
    #[inline(always)]
    pub fn get_ca_cb_angle(&self, idx1: usize, idx2: usize, return_radian: bool) -> Option<f32> {
        let ca1 = self.get_ca(idx1);
        let cb1 = self.get_cb(idx1);
        let ca2 = self.get_ca(idx2);
        let cb2 = self.get_cb(idx2);
        if ca1.is_some() && cb1.is_some() && ca2.is_some() && cb2.is_some() {
            let angle = ca1.expect("Unable to get CA1 coordinate").calc_angle(
                &cb1.expect("Unable to get CB1 coordinate"),
                &ca2.expect("Unable to get CA2 coordinate"),
                &cb2.expect("Unable to get CB2 coordinate"),
                return_radian,
            );
            Some(angle)
        } else {
            None
        }
    }
    #[inline(always)]
    pub fn get_res_name(&self, idx: usize) -> &[u8; 3] {
        if idx >= self.num_residues {
            return b"UNK";
        }
        &self.residue_name[idx]
    }
    #[inline(always)]
    pub fn get_res_serial(&self, idx1: usize, idx2: usize) -> (u64, u64) {
        (self.residue_serial[idx1], self.residue_serial[idx2])
    }
    #[inline(always)]
    pub fn get_ppf(&self, idx1: usize, idx2: usize, dist_cutoff: f32) -> Option<[f32; 4]> {
        let ca1 = self.get_ca(idx1);
        let cb1 = self.get_cb(idx1);
        let cb2 = self.get_cb(idx2);
        if ca1.is_some() && cb1.is_some() && cb2.is_some() {
            let rel_cb1 = cb1.unwrap().sub(&ca1.unwrap());
            let rel_cb2 = cb2.unwrap().sub(&ca1.unwrap());
            let ppf = rel_cb1.get_ppf(&rel_cb2);
            if ppf[0] > dist_cutoff {
                return None;
            }
            Some(ppf)
        } else {
            None
        }
    }
    #[inline(always)]
    pub fn get_trrosetta_feature(&self, idx1: usize, idx2: usize, dist_cutoff: f32) -> Option<(f32, f32, f32, f32, f32, f32)> {
        let ca1 = self.get_ca(idx1);
        let ca2 = self.get_ca(idx2);
        let cb1 = self.get_cb(idx1);
        let cb2 = self.get_cb(idx2);
        let n1 = self.get_n(idx1);
        let n2 = self.get_n(idx2);
        if let (Some(ca1), Some(ca2), Some(cb1), Some(cb2), Some(n1), Some(n2)) =
            (ca1, ca2, cb1, cb2, n1, n2)
        {
            let cb_dist = cb1.calc_distance(&cb2);
            if cb_dist > dist_cutoff {
                return None;
            }
            let omega = calc_torsion_radian(&ca1, &cb1, &cb2, &ca2);
            let theta1 = calc_torsion_radian(&n1, &ca1, &cb1, &cb2);
            let theta2 = calc_torsion_radian(&cb1, &cb2, &ca2, &n2);
            let phi1 = calc_angle_radian(&ca1, &cb1, &cb2);
            let phi2 = calc_angle_radian(&cb1, &cb2, &ca2);
            Some((cb_dist, omega, theta1, theta2, phi1, phi2))
        } else {
            None
        }
    }
    #[inline(always)]
    pub fn get_trrosetta_feature2(&self, idx1: usize, idx2: usize) -> Option<[f32; 7]> {
        let ca1 = self.get_ca(idx1);
        let ca2 = self.get_ca(idx2);
        let cb1 = self.get_cb(idx1);
        let cb2 = self.get_cb(idx2);
        let n1 = self.get_n(idx1);
        let n2 = self.get_n(idx2);
        if let (Some(ca1), Some(ca2), Some(cb1), Some(cb2), Some(n1), Some(n2)) =
            (ca1, ca2, cb1, cb2, n1, n2)
        {
            let log_dist = idx2 as f32 - idx1 as f32;
            let dist_sign = if idx2 > idx1 { 1.0 } else { -1.0 };
            let log_dist = log_dist.abs().ln() * dist_sign;
            let cb_dist = cb1.calc_distance(&cb2);
            let omega = calc_torsion_radian(&ca1, &cb1, &cb2, &ca2);
            let theta1 = calc_torsion_radian(&n1, &ca1, &cb1, &cb2);
            let theta2 = calc_torsion_radian(&cb1, &cb2, &ca2, &n2);
            let phi1 = calc_angle_radian(&ca1, &cb1, &cb2);
            let phi2 = calc_angle_radian(&cb1, &cb2, &ca2);
            let feature = [omega, theta1, theta2, phi1, phi2, cb_dist, log_dist];
            Some(feature)
        } else {
            None
        }
    }

    #[inline(always)]
    pub fn get_pdb_tr_feature(&self, idx1: usize, idx2: usize, dist_cutoff: f32) -> Option<(f32, f32, f32, f32, f32)> {
        let ca1 = self.get_ca(idx1);
        let ca2 = self.get_ca(idx2);
        let cb1 = self.get_cb(idx1);
        let cb2 = self.get_cb(idx2);
        let n1 = self.get_n(idx1);
        let n2 = self.get_n(idx2);
        if let (Some(ca1), Some(ca2), Some(cb1), Some(cb2), Some(n1), Some(n2)) =
            (ca1, ca2, cb1, cb2, n1, n2)
        {

            let ca_dist = ca1.calc_distance(&ca2);
            if ca_dist > dist_cutoff {
                return None;
            }
            let cb_dist = cb1.calc_distance(&cb2);
            // let ca_cb_angle = self.get_ca_cb_angle(idx1, idx2, true).unwrap();
            let ca_cb_angle = ca1.calc_angle(&cb1, &ca2, &cb2, true);
            let theta1 = calc_torsion_radian(&n1, &ca1, &cb1, &cb2);
            let theta2 = calc_torsion_radian(&cb1, &cb2, &ca2, &n2);
            Some((ca_dist, cb_dist, ca_cb_angle, theta1, theta2))
            // Some(feature)
        } else {
            None
        }
    }
    
    
    
    #[inline(always)]
    pub fn get_bfactor(&self, idx: usize) -> f32 {
        self.b_factors[idx]
    }
    #[inline(always)]
    pub fn get_plddt(&self, idx: usize) -> f32 {
        self.get_bfactor(idx)
    }
    #[inline(always)]
    pub fn get_avg_bfactor(&self) -> f32 {
        let mut sum = 0.0;
        for i in 0..self.num_residues {
            sum += self.b_factors[i];
        }
        sum / self.num_residues as f32
    }
    #[inline(always)]
    pub fn get_avg_plddt(&self) -> f32 {
        self.get_avg_bfactor()
    }
    #[inline(always)]
    pub fn get_list_amino_acids_and_distances(&self, i: usize, j: usize) -> Option<(u8, u8, f32)> {
        // Return i, j, aa_i, aa_j, distance
        let aa_i = map_aa_to_u8(self.get_res_name(i));
        let aa_j = map_aa_to_u8(self.get_res_name(j));
        let distance = self.get_ca_distance(i, j);
        if distance.is_none() {
            None
        } else {
            let distance = distance.unwrap();
            if distance <= 20.0 {
                Some((aa_i, aa_j, distance))
            } else {
                None
            }
        }
    }
    
}

#[cfg(test)]
mod structure_tests {
    #[test]
    fn test_gly_integration() {
        let data = crate::structure::io::pdb::Reader::from_file("data/homeobox/1akha-.pdb")
            .expect("Unable to read test file");
        let structure = &data.read_structure().expect("Unable to read structure");
        let compact = &structure.to_compact();

        for idx in 0..compact.num_residues {
            let res_name = compact
                .residue_name
                .get(idx)
                .expect("expected residue name");
            if res_name == b"GLY" {
                let res_str = std::str::from_utf8(res_name).expect("expected residue name");
                let ca = compact.get_ca(idx);
                let cb = compact.get_cb(idx);
                println!(
                    "residue {} ({}) is {} CA: {:?} CB: {:?}",
                    idx, compact.residue_serial[idx], res_str, ca, cb
                );
                assert!(cb.is_some());
            }
        }
        println!("{:?}", compact.get_res_name(compact.num_residues - 1));
        assert_eq!(compact.num_residues, structure.num_residues);
    }
    
    #[test]
    fn test_avg_bfactor() {
        let data = crate::structure::io::pdb::Reader::from_file("data/homeobox/1akha-.pdb")
            .expect("Unable to read test file");
        let structure = &data.read_structure().expect("Unable to read structure");
        let compact = &structure.to_compact();
        let avg_bfactor = compact.get_avg_bfactor();
        println!("Average B-factor: {}", avg_bfactor);
        assert!(avg_bfactor > 0.0);
    }
}
