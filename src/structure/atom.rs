pub struct Atom {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub atom_name: [u8; 4],
    pub atom_serial: u64,
    pub res_name: [u8; 4],
    pub res_serial: u64,
    pub chain: u8,
    pub b_factor: f64,
}

impl Atom {
    pub fn new(
        x: f64, y: f64, z: f64, atom_name: [u8; 4], atom_serial: u64,
        res_name: [u8; 4], res_serial: u64, chain: u8, b_factor: f64
    ) -> Atom {
        Atom {
            x: x, y: y, z: z,
            atom_name: atom_name, atom_serial: atom_serial,
            res_name: res_name, res_serial: res_serial,
            chain: chain, b_factor: b_factor,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CoordinateVector {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub size: usize,
}

/// AtomVector
#[derive(Debug, Clone)]
pub struct AtomVector {
    pub coordinates: CoordinateVector,
    pub atom_name: Vec<[u8; 4]>,
    pub atom_serial: Vec<u64>,
    pub res_name: Vec<[u8; 4]>,
    pub res_serial: Vec<u64>,
    pub chain: Vec<u8>,
    pub b_factor: Vec<f64>,
}

impl AtomVector {
    pub fn new() -> AtomVector {
        AtomVector {
            atom_name: Vec::new(),
            coordinates: CoordinateVector {
                x: Vec::new(),
                y: Vec::new(),
                z: Vec::new(),
                size: 0,
            },
            atom_serial: Vec::new(),
            res_name: Vec::new(),
            res_serial: Vec::new(),
            chain: Vec::new(),
            b_factor: Vec::new(),
        }
    }

    pub fn push(
        &mut self, atom_name: [u8; 4], x: f64, y: f64, z: f64,
        atom_serial: u64, res_name: [u8; 4], res_serial: u64,
        chain: u8, b_factor: f64
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
        self.b_factor.push(atom.b_factor);
    }

    pub fn get(&self, index: usize) -> Atom {
        Atom {
            atom_name: self.atom_name[index],
            x: self.coordinates.x[index],
            y: self.coordinates.y[index],
            z: self.coordinates.z[index],
            atom_serial: self.atom_serial[index],
            res_name: self.res_name[index],
            res_serial: self.res_serial[index],
            chain: self.chain[index],
            b_factor: self.b_factor[index],
        }
    }

}
