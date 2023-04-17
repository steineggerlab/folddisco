use crate::structure::coordinate::{Coordinate, CoordinateVector};

#[derive(Debug)]
pub struct Atom {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub atom_name: [u8; 4],
    pub atom_serial: u64,
    pub res_name: [u8; 3],
    pub res_serial: u64,
    pub chain: u8,
    pub b_factor: f32,
}

impl Atom {
    // Constructor
    pub fn new(
        x: f32,
        y: f32,
        z: f32,
        atom_name: [u8; 4],
        atom_serial: u64,
        res_name: [u8; 3],
        res_serial: u64,
        chain: u8,
        b_factor: f32,
    ) -> Atom {
        Atom {
            x,
            y,
            z,
            atom_name,
            atom_serial,
            res_name,
            res_serial,
            chain,
            b_factor,
        }
    }
    pub fn new_empty() -> Atom {
        Atom {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            atom_name: [0; 4],
            atom_serial: 0,
            res_name: [0; 3],
            res_serial: 0,
            chain: 0,
            b_factor: 0.0,
        }
    }
    pub fn is_empty(&self) -> bool {
        self.atom_serial == 0 && self.res_serial == 0
    }
    pub fn get_coordinate(&self) -> Coordinate {
        Coordinate::new(self.x, self.y, self.z)
    }
    pub fn get_res_name(&self) -> [u8;3] {
        self.res_name
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
        &mut self,
        atom_name: [u8; 4],
        x: f32,
        y: f32,
        z: f32,
        atom_serial: u64,
        res_name: [u8; 3],
        res_serial: u64,
        chain: u8,
        b_factor: f32,
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

    pub fn is_ca(&self, index: usize) -> bool {
        &self.atom_name[index] == b" CA "
    }

    pub fn is_cb(&self, index: usize) -> bool {
        &self.atom_name[index] == b" CB "
    }

    pub fn is_c(&self, index: usize) -> bool {
        &self.atom_name[index] == b" C  "
    }

    pub fn is_n(&self, index: usize) -> bool {
        &self.atom_name[index] == b" N  "
    }

    pub fn is_backbone(&self, index: usize) -> bool {
        self.is_ca(index) || self.is_c(index) || self.is_n(index)
    }

    pub fn get_coordinates(&self, index: usize) -> Coordinate {
        Coordinate {
            x: self.coordinates.x[index],
            y: self.coordinates.y[index],
            z: self.coordinates.z[index],
        }
    }

    pub fn get_res_serial(&self, index: usize) -> u64 {
        self.res_serial[index]
    }

    pub fn get_res_name(&self, index: usize) -> [u8; 3] {
        self.res_name[index]
    }

    pub fn get_atom_name(&self, index: usize) -> [u8; 4] {
        self.atom_name[index]
    }

    pub fn len(&self) -> usize {
        self.coordinates.size
    }

    pub fn _print_residue(&self, index: usize) -> String {
        let res_byte = self.get_res_name(index);
        let res_str = std::str::from_utf8(&res_byte).unwrap();
        String::from(res_str)
    }

    pub fn _print_atom(&self, index: usize) {
        let atom_name = self.get_atom_name(index);
        println!(
            "{}{}{}{}",
            atom_name[0] as char, atom_name[1] as char, atom_name[2] as char, atom_name[3] as char
        );
    }
}
