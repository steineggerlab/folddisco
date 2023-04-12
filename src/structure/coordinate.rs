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
}
