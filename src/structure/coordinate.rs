use crate::utils::calculator::Calculate;
use crate::utils::constants::CA_CB_DIST;
use crate::structure::core::Structure;

use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Coordinate {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Coordinate {
    pub fn new(x: f32, y: f32, z: f32) -> Coordinate {
        Coordinate { x, y, z }
    }
    pub fn build(x: &Option<f32>, y: &Option<f32>, z: &Option<f32>) -> Self {
        Coordinate {
            x: x.expect("Unable to get coordinate"),
            y: y.expect("Unable to get coordinate"),
            z: z.expect("Unable to get coordinate"),
        }
    }
    pub fn add(&self, other: &Coordinate) -> Coordinate {
        Coordinate {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
    pub fn sub(&self, other: &Coordinate) -> Coordinate {
        Coordinate {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
    pub fn dot(&self, other: &Coordinate) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    pub fn cross(&self, other: &Coordinate) -> Coordinate {
        Coordinate {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
    pub fn norm(&self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
    pub fn normalize(&self) -> Coordinate {
        let norm = self.norm();
        Coordinate {
            x: self.x / norm,
            y: self.y / norm,
            z: self.z / norm,
        }
    }
    pub fn scale(&self, factor: f32) -> Coordinate {
        Coordinate {
            x: self.x * factor,
            y: self.y * factor,
            z: self.z * factor,
        }
    }
    pub fn distance(&self, other: &Coordinate) -> f32 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

impl Calculate for Coordinate {
    fn calc_distance(&self, other: &Coordinate) -> f32 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        dist
    }

    fn calc_angle(&self, atom2: &Coordinate, atom3: &Coordinate, atom4: &Coordinate) -> f32 {
        let (a, b, c, d) = (self, atom2, atom3, atom4);
        // Form vectors
        let v1 = b.sub(a); // vector 1
        let v2 = d.sub(c); // vector 2
        let dot = v1.dot(&v2); // dot product
        let v1_len = v1.norm(); // length of vector 1
        let v2_len = v2.norm(); // length of vector 2
        let cos = dot / (v1_len * v2_len); // cos of angle
        let radian = cos.acos(); // angle in radians
        let degree = radian.to_degrees(); // angle in degrees
        degree
    }

    fn calc_torsion_angle(&self, atom2: &Coordinate, atom3: &Coordinate, atom4: &Coordinate) -> f32 {
        let (ca, b, c, d) = (self, atom2, atom3, atom4);
        // If psi angle, b: N(i), ca: CA(i), c: C(i), d: N(i+1)
        // If phi angle, b: C(i-1), c: N(i), ca: CA(i), d: C(i)

        // Form vectors
        // FIXME: order of atoms do not match
        let v1 = b.sub(ca);
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

}

// Originally from foldseek StructureTo3DiBase::approxCBetaPosition
// link: https://github.com/steineggerlab/foldseek/blob/master/lib/3di/structureto3di.cpp
pub fn approx_cb(ca: &Coordinate, n: &Coordinate, c: &Coordinate) -> Coordinate {
    // Assumption: CA forms with its four ligands a tetrahedral.
    let v1 = c.sub(ca).normalize();
    let v2 = n.sub(ca).normalize();

    let b1 = v2.add(&v1.scale(1.0 / 3.0));
    let b2 = v1.cross(&b1);

    let u1 = b1.normalize();
    let u2 = b2.normalize();

    //direction from c_alpha to c_beta
    let v4 = u1.scale(-1.0 / 2.0).sub(&u2.scale(3.0f32.sqrt() / 2.0));
    let v4 = v4.scale(8.0f32.sqrt() / 3.0);
    let v4 = v4.add(&v1.scale(-1.0 / 3.0));

    let cb = ca.add(&v4.scale(CA_CB_DIST));

    cb
}

#[derive(Debug, Clone)]
pub struct CoordinateVector {
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Vec<f32>,
    pub size: usize,
}

impl CoordinateVector {
    pub fn new() -> CoordinateVector {
        CoordinateVector {
            x: Vec::new(),
            y: Vec::new(),
            z: Vec::new(),
            size: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CarbonCoordinateVector {
    x: Vec<Option<f32>>,
    y: Vec<Option<f32>>,
    z: Vec<Option<f32>>,
}

impl CarbonCoordinateVector {
    pub fn new() -> Self {
        CarbonCoordinateVector {
            x: Vec::new(),
            y: Vec::new(),
            z: Vec::new(),
        }
    }
    pub fn get(&self, idx: usize) -> (Option<f32>, Option<f32>, Option<f32>) {
        (self.x[idx], self.y[idx], self.z[idx])
    }

    pub fn push(&mut self, coordinate: &Coordinate) {
        self.x.push(Some(coordinate.x));
        self.y.push(Some(coordinate.y));
        self.z.push(Some(coordinate.z));
    }

    pub fn push_none(&mut self) {
        self.x.push(None);
        self.y.push(None);
        self.z.push(None);
    }

}

#[derive(Debug, Clone)]
pub enum TorsionType {
    Psi,
    Phi,
    None,
}

#[derive(Debug, Clone)]
pub struct Torsion {
    pub torsion_type: TorsionType,
    pub torsion_map: HashMap<u64,f32>,
}

impl Torsion {
    pub fn new() -> Self {
        Torsion {
            torsion_type: TorsionType::None,
            torsion_map: HashMap::new(),
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
        self.torsion_map.insert(res_serial, torsion);
    }

    pub fn get(&self, res_serial: u64) -> Option<&f32> {
        self.torsion_map.get(&res_serial)
    }

    pub fn calc_torsion_angle(atom1: &Coordinate, atom2: &Coordinate, atom3: &Coordinate, atom4: &Coordinate) -> f32 {
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
        let n2 = structure.atom_vector.get_nth_n(nth+1).get_coordinate();
    
        let psi = Torsion::calc_torsion_angle(&n, &ca,&c, &n2);
        psi
    }
    // TODO: implement get_torsion_angle (phi)
}

#[cfg(test)]
mod coordinate_tests {
    use super::*;
    #[test]
    fn test_approx_cb() {
        let n = Coordinate {
            x: 199.117,
            y: 245.223,
            z: 222.801,
        };
        let ca = Coordinate {
            x: 200.327,
            y: 244.427,
            z: 222.675,
        };
        let c = Coordinate {
            x: 200.278,
            y: 243.297,
            z: 223.701,
        };
        let actual_cb = Coordinate {
            x: 201.553,
            y: 245.296,
            z: 222.934,
        };

        let test_cb = approx_cb(&ca, &n, &c);

        println!("actual_cb: {:?}", actual_cb);
        println!("test_cb: {:?}", test_cb);
        println!("distance: {:?}", actual_cb.distance(&test_cb));
    }

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
        let test_torsion = a.calc_torsion_angle(&b, &c, &d);

        println!("actual_torsion: {:?}", actual_torsion);
        println!("test_torsion: {:?}", test_torsion);
    }
}
