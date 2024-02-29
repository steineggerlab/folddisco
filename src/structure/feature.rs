use crate::structure::coordinate::Coordinate;
use crate::structure::core::Structure;

#[derive(Debug, Clone)]
pub enum TorsionType {
    Psi,
    Phi,
    None,
}

#[derive(Debug, Clone)]
pub struct Torsion {
    pub torsion_type: TorsionType,
    pub residue: Vec<u64>,
    pub torsion: Vec<f32>,
}

impl Torsion {
    pub fn new() -> Self {
        Torsion {
            torsion_type: TorsionType::None,
            residue: Vec::new(),
            torsion: Vec::new(),
        }
    }

    pub fn build(structure: &Structure, torsiontype: TorsionType) -> Self {
        let mut torsion = Torsion::new();
        torsion.set_torsion_type(torsiontype);

        for n in 0..structure.num_residues {
            let res_serial = structure.atom_vector.get_nth_n(n).get_res_serial();
            let psi = Torsion::calc_psi_angle(structure, n);
            torsion.push(res_serial, psi);
        }

        torsion
    }

    pub fn set_torsion_type(&mut self, torsion_type: TorsionType) {
        self.torsion_type = torsion_type;
    }

    pub fn get_torsion_type(&self) -> &TorsionType {
        &self.torsion_type
    }

    pub fn push(&mut self, res_serial: u64, torsion: f32) {
        self.residue.push(res_serial);
        self.torsion.push(torsion);
    }

    pub fn calc_torsion_angle(
        atom1: &Coordinate,
        atom2: &Coordinate,
        atom3: &Coordinate,
        atom4: &Coordinate,
    ) -> f32 {
        let (a, b, c, d) = (atom1, atom2, atom3, atom4);
        /* If psi angle, b: N(i), ca: CA(i), c: C(i), d: N(i+1)
        If phi angle, b: C(i-1), c: N(i), ca: CA(i), d: C(i) */

        // Form vectors
        let v1 = b.sub(a);
        let v2 = c.sub(b);
        let v3 = d.sub(c);

        // Form normal vectors via cross products
        let r = v1.cross(&v2).normalize();
        let s = v2.cross(&v3).normalize();
        let t = r.cross(&v2.normalize()).normalize();
        let x = r.dot(&s);
        let y = s.dot(&t);
        -y.atan2(x).to_degrees()
    }

    pub fn calc_psi_angle(structure: &Structure, nth: usize) -> f32 {
        let psi_vec: Vec<f32> = Vec::new();

        let n = structure.atom_vector.get_nth_n(nth).get_coordinate();
        let ca = structure.atom_vector.get_nth_ca(nth).get_coordinate();
        let c = structure.atom_vector.get_nth_c(nth).get_coordinate();
        let n2 = structure.atom_vector.get_nth_n(nth + 1).get_coordinate();

        let psi = Torsion::calc_torsion_angle(&n, &ca, &c, &n2);
        psi
    }
    // TODO: implement get_torsion_angle (phi)
}

#[cfg(tests)]
mod tests {
    use super::*;
    #[test]
    fn test_calc_torsion_angle() {
        let a = Coordinate {
            x: 24.969,
            y: 13.428,
            z: 30.692,
        };
        let b = Coordinate {
            x: 24.044,
            y: 12.661,
            z: 29.808,
        };
        let c = Coordinate {
            x: 22.785,
            y: 13.482,
            z: 29.543,
        };
        let d = Coordinate {
            x: 21.951,
            y: 13.670,
            z: 30.431,
        };

        let actual_torsion = -71.21515;
        let test_torsion = calc_torsion_angle(&a, &b, &c, &d);

        println!("actual_torsion: {:?}", actual_torsion);
        println!("test_torsion: {:?}", test_torsion);
    }
}
