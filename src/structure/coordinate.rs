use crate::utils::calculator::Calculate;
use crate::utils::constants::CA_CB_DIST;
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

    pub fn get_ppf(&self, other: &Coordinate) -> [f32; 4] {
        let n1 = self.normalize();
        let n2=other.normalize();
        let d=other.sub(self);
        let nd=d.normalize();
        let n1_nd = n1.dot(&nd).acos().to_degrees();
        let n2_nd = n2.dot(&nd).acos().to_degrees();
        let n1_n2 = n1.dot(&n2).acos().to_degrees();
        [d.norm(), n1_nd, n2_nd, n1_n2]
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
        let v1 = (b.x - a.x, b.y - a.y, b.z - a.z); // vector 1
        let v2 = (d.x - c.x, d.y - c.y, d.z - c.z); // vector 2
        let dot = v1.0 * v2.0 + v1.1 * v2.1 + v1.2 * v2.2; // dot product
        let v1_len = (v1.0.powf(2.0) + v1.1.powf(2.0) + v1.2.powf(2.0)).sqrt(); // length of vector 1
        let v2_len = (v2.0.powf(2.0) + v2.1.powf(2.0) + v2.2.powf(2.0)).sqrt(); // length of vector 2
        let cos = dot / (v1_len * v2_len); // cos of angle
        let radian = cos.acos(); // angle in radians
        let degree = radian.to_degrees(); // angle in degrees
        degree
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

pub fn calc_torsion_angle(
    a: &Coordinate, b: &Coordinate, c: &Coordinate, d: &Coordinate,
) -> f32 {
    let v1 = b.sub(a);
    let v2 = c.sub(b);
    let v3 = d.sub(c);

    let r = v1.cross(&v2).normalize();
    let s = v2.cross(&v3).normalize();
    let t = r.cross(&v2.normalize()).normalize();
    let x = r.dot(&s);
    let y = s.dot(&t);
    -y.atan2(x).to_degrees()
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
    pub fn get(&self, idx: usize) -> (f32, f32, f32) {
        (self.x[idx], self.y[idx], self.z[idx])
    }
    pub fn calc_torsion_angle(
        &self,
        a: usize,
        b: usize,
        c: usize,
        d: usize,
    ) -> f32 {
        let (a_x, a_y, a_z) = self.get(a);
        let (b_x, b_y, b_z) = self.get(b);
        let (c_x, c_y, c_z) = self.get(c);
        let (d_x, d_y, d_z) = self.get(d);

        let a = Coordinate {x: a_x, y: a_y, z: a_z};
        let b = Coordinate {x: b_x, y: b_y, z: b_z};
        let c = Coordinate {x: c_x, y: c_y, z: c_z};
        let d = Coordinate {x: d_x, y: d_y, z: d_z};
        calc_torsion_angle(&a, &b, &c, &d)
    }

    pub fn calc_all_torsion_angles(&self) -> Vec<f32> {
        let mut torsion_angles = Vec::new();
        for i in 0..self.x.len() - 3 {
            let a = i;
            let b = i + 1;
            let c = i + 2;
            let d = i + 3;
            torsion_angles.push(self.calc_torsion_angle(a, b, c, d));
        }
        torsion_angles
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

    pub fn calc_torsion_angle(
        &self,
        a: usize,
        b: usize,
        c: usize,
        d: usize,
    ) -> Option<f32> {
        let (a_x, a_y, a_z) = self.get(a);
        let (b_x, b_y, b_z) = self.get(b);
        let (c_x, c_y, c_z) = self.get(c);
        let (d_x, d_y, d_z) = self.get(d);

        if a_x.is_none() || a_y.is_none() || a_z.is_none()
        || b_x.is_none() || b_y.is_none() || b_z.is_none()
        || c_x.is_none() || c_y.is_none() || c_z.is_none()
        || d_x.is_none() || d_y.is_none() || d_z.is_none()
        {
            return None;
        }

        let a = Coordinate {
            x: a_x.unwrap(),
            y: a_y.unwrap(),
            z: a_z.unwrap(),
        };
        let b = Coordinate {
            x: b_x.unwrap(),
            y: b_y.unwrap(),
            z: b_z.unwrap(),
        };
        let c = Coordinate {
            x: c_x.unwrap(),
            y: c_y.unwrap(),
            z: c_z.unwrap(),
        };
        let d = Coordinate {
            x: d_x.unwrap(),
            y: d_y.unwrap(),
            z: d_z.unwrap(),
        };

        Some(calc_torsion_angle(&a, &b, &c, &d))
    }

    pub fn calc_all_torsion_angles(&self) -> Vec<Option<f32>> {
        let mut torsion_angles = Vec::new();
        for i in 0..self.x.len() - 3 {
            let a = i;
            let b = i + 1;
            let c = i + 2;
            let d = i + 3;
            torsion_angles.push(self.calc_torsion_angle(a, b, c, d));
        }
        torsion_angles
    }

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

        let actual_phi = -71.21515;
        let test_phi = calc_torsion_angle(&a, &b, &c, &d);

        println!("actual_phi: {:?}", actual_phi);
        println!("test_phi: {:?}", test_phi);
    }

    #[test]
    fn test_calc_torsion_angle_vec() {
        let coord_vec = CoordinateVector {
            x: vec![24.969, 24.044, 22.785, 21.951],
            y: vec![13.428, 12.661, 13.482, 13.670],
            z: vec![30.692, 29.808, 29.543, 30.431],
            size: 4,
        };
        let torsion_vec = coord_vec.calc_all_torsion_angles();
        println!("torsion_vec: {}", torsion_vec[0]);
    }


}
