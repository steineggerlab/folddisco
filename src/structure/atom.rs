use crate::utils::calculator::Calculate;
#[derive(Debug, Clone)]
pub struct Coordinate {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Coordinate {
    pub fn build(x:&Option<f32>, y:&Option<f32>, z:&Option<f32>) -> Self {
        Coordinate { x: x.unwrap() , y: y.unwrap() , z: z.unwrap()  }
    }
}

impl Calculate for Coordinate {
    fn calc_distance(&self, other: &Coordinate) -> f32{
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        dist
    }
    
    fn calc_dihedral(&self, atom2: &Coordinate, atom3: &Coordinate, atom4: &Coordinate) -> f32 {
        // code referred from pdbtbx
        // https://github.com/douweschulte/pdbtbx

        let (a,b,c,d) = (self, atom2, atom3, atom4);
        // Form vectors
        let ba = [a.x - b.x, a.y - b.y, a.z - b.z];
        let bc = [c.x - b.x, c.y - b.y, c.z - b.z];
        let cb = [b.x - c.x, b.y - c.y, b.z - c.z];
        let cd = [d.x - c.x, d.y - c.y, d.z - c.z];

        // Form two normal vectors via cross products
        let n1 = [
            ba[1] * bc[2] - ba[2] * bc[1],
            ba[2] * bc[0] - ba[0] * bc[2],
            ba[0] * bc[1] - ba[1] * bc[0],
        ];
        let n2 = [
            cb[1] * cd[2] - cb[2] * cd[1],
            cb[2] * cd[0] - cb[0] * cd[2],
            cb[0] * cd[1] - cb[1] * cd[0],
        ];

        // calculate abs of vecs
        let abs_n1 = n1.iter().fold(0.0, |acc, x| acc + (x * x)).sqrt();
        let abs_n2 = n2.iter().fold(0.0, |acc, x| acc + (x * x)).sqrt();

        let dot = n1
            .iter()
            .zip(n2.iter())
            .fold(0.0, |acc, (a, b)| acc + (a * b));
        (dot / (abs_n1 * abs_n2)).acos().to_degrees()
    }
}


#[derive(Debug)]
pub struct Atom {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub atom_name: [u8;4],
    pub atom_serial: u64,
    pub res_name: [u8;3],     
    pub res_serial: u64,
    pub chain: u8,
    pub b_factor: f32,
}

impl Atom {
    pub fn build(
        x: f32, y: f32, z: f32, atom_name: [u8;4], atom_serial: u64,
        res_name: [u8;3], res_serial: u64, chain: u8, b_factor: f32
    ) -> Atom {
        Atom {
            x, y,z, 
            atom_name, atom_serial,
            res_name, res_serial,
            chain, b_factor,
        }
    }
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
        CoordinateVector { x: Vec::new(), y: Vec::new(), z: Vec::new(), size: 0 }
    }
}

/// AtomVector
#[derive(Debug, Clone)]
pub struct AtomVector {
    pub coordinates: CoordinateVector,
    pub atom_name: Vec<[u8; 4]>,
    pub atom_serial: Vec<u64>,
    pub res_name: Vec<[u8; 3]>,
    pub res_serial: Vec<u64>,
    pub chain: Vec<u8>,
    pub b_factor: Vec<f32>,
}

impl AtomVector {
    pub fn new() -> AtomVector {
        AtomVector {
            atom_name: Vec::new(),
            coordinates: CoordinateVector::new(),
            atom_serial: Vec::new(),
            res_name: Vec::new(),
            res_serial: Vec::new(),
            chain: Vec::new(),
            b_factor: Vec::new(),
        }
    }

    pub fn push(
        &mut self, atom_name: [u8; 4], x: f32, y: f32, z: f32,
        atom_serial: u64, res_name: [u8; 3], res_serial: u64,
        chain: u8, b_factor: f32
    ) {
        self.atom_name.push(atom_name);
        self.coordinates.x.push(x);
        self.coordinates.y.push(y);
        self.coordinates.z.push(z);
        self.coordinates.size += 1;
        self.atom_serial.push(atom_serial);
        self.res_name.push(res_name);
        self.res_serial.push(res_serial);
        self.chain.push(chain);
        self.b_factor.push(b_factor);
    }

    pub fn push_atom(&mut self, atom: Atom) {
        self.atom_name.push(atom.atom_name);
        self.coordinates.x.push(atom.x);
        self.coordinates.y.push(atom.y);
        self.coordinates.z.push(atom.z);
        self.coordinates.size += 1;
        self.atom_serial.push(atom.atom_serial);
        self.res_name.push(atom.res_name);
        self.res_serial.push(atom.res_serial);
        self.chain.push(atom.chain);
        self.b_factor.push(atom.b_factor);
    }

    pub fn get(&self, index: usize) -> Atom {
        // ??? index should be 1--n ???
        Atom {
            // atom_name: self.atom_name[index].clone(),
            atom_name: self.atom_name[index],
            x: self.coordinates.x[index],
            y: self.coordinates.y[index],
            z: self.coordinates.z[index],
            atom_serial: self.atom_serial[index],
            // res_name: self.res_name[index].clone(),
            res_name: self.res_name[index],
            res_serial: self.res_serial[index],
            chain: self.chain[index],
            b_factor: self.b_factor[index],
        }
    }

    pub fn is_CA(&self, index: usize) -> bool {
        &self.atom_name[index] == b" CA "
    }

    pub fn is_CB(&self, index: usize) -> bool {
        &self.atom_name[index] == b" CB "
    }

    pub fn get_coordinates(&self, index: usize) -> Coordinate {
        Coordinate {
            x: self.coordinates.x[index],
            y: self.coordinates.y[index],
            z: self.coordinates.z[index],
        }
    }

    pub fn get_res_serial(&self, index:usize) -> u64 {
        self.res_serial[index]
    }

}
