use crate::utils::calculator::Calculate;
use crate::structure::coordinate::{Coordinate, approx_cb};
use crate::structure::atom::{Atom, AtomVector};

/// Structure is the main data structure for storing the information of a protein structure.
#[derive(Debug)]
pub struct Structure {
    pub num_chains: usize,
    pub chains: Vec<u8>,
    pub atom_vector : AtomVector,
    // pub atom_vectors: Vec<AtomVector>, // ???왜 AtomVector이 아닌 Vec<AtomVector>???
    pub num_atoms: usize,
    pub num_residues: usize,
}

impl Structure {
    pub fn new() -> Structure {
        Structure {
            num_chains : 0,
            chains : Vec::new(),
            atom_vector : AtomVector::new(),
            num_atoms : 0,
            num_residues : 0,
        }
    }

    pub fn update(&mut self, atom : Atom, record: &mut (u8, u64) ) {
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

    fn _fill_gly_cbeta(&mut self) {
        // Iterate over self.atom_vector

        let mut curr_res_serial: u64 = 0;

        let mut gly_n = Atom::new_empty();
        let mut gly_ca = Atom::new_empty();
        let mut gly_c = Atom::new_empty();
        let mut gly_o = Atom::new_empty();

        for i in 0..self.atom_vector.len() {
            if self.atom_vector.res_name.get(i) == Some(b"GLY") {
                match self.atom_vector.res_serial.get(i) {
                    Some(res_serial) => {
                        if curr_res_serial != *res_serial {
                            curr_res_serial = *res_serial;
                            gly_n = Atom::new_empty();
                            gly_ca = Atom::new_empty();
                            gly_c = Atom::new_empty();
                            gly_o = Atom::new_empty();
                        }
                    },
                    None => (),
                }
                match self.atom_vector.atom_name.get(i) {
                    Some(b" N  ") => gly_n = self.atom_vector.get(i),
                    Some(b" CA ") => gly_ca = self.atom_vector.get(i),
                    Some(b" C  ") => gly_c = self.atom_vector.get(i),
                    Some(b" O  ") => gly_o = self.atom_vector.get(i),
                    _ => (),
                }
                if gly_n.is_empty() || gly_ca.is_empty() || gly_c.is_empty() || gly_o.is_empty() {
                    continue;
                } else {
                    let cbeta_coord = approx_cb(
                        &gly_ca.get_coordinate(), &gly_n.get_coordinate(), &gly_c.get_coordinate()
                    );
                    // // Make new cbeta atom
                    // let cbeta_atom = Atom::new(cbeta_coord.x, cbeta_coord.y, cbeta_coord.z, b"CB  ", b"GLY", gly_ca.chain, gly_ca.res_serial, gly_ca.res_name, );

                }
            }

        }
    }



    pub fn to_compact(&self) -> CompactStructure {
        CompactStructure::build(self)
    }

    // pub fn count_chains() {}
    // pub fn count_atoms() {}
    // pub fn count_residues() {}

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
            z: Vec::new() }
    }
    pub fn get(&self, idx: usize) -> (Option<f32>, Option<f32>, Option<f32>) {
        (self.x[idx], self.y[idx], self.z[idx])
    }
}

#[derive(Debug, Clone)]
pub struct CompactStructure {
    pub num_chains: usize,
    pub chains: Vec<u8>,
    pub num_residues : usize,
    pub residues: Vec<u64>,
    pub ca_vector: CarbonCoordinateVector,
    pub cb_vector: CarbonCoordinateVector,
}

impl CompactStructure {
    pub fn build(origin: &Structure) -> CompactStructure {

        let model = &origin.atom_vector;

        let mut res_vec: Vec<u64> = Vec::new();
        let mut ca_vec = CarbonCoordinateVector::new();
        let mut cb_vec = CarbonCoordinateVector::new();
        let mut prev_res_serial: Option<u64> = None;

        for idx in 0..origin.num_atoms {
            // TODO : implement iterator for Structure.atom_vector
            // TODO : revise using match, .... etc
            let atom = model.get(idx);

            if &atom.atom_name == b" CA " {
                if prev_res_serial == Some(atom.res_serial) {
                    ca_vec.x.last_mut().unwrap().replace(atom.x);
                    ca_vec.y.last_mut().unwrap().replace(atom.y);
                    ca_vec.z.last_mut().unwrap().replace(atom.z);
                } else {
                    res_vec.push(atom.res_serial);
                    ca_vec.x.push(Some(atom.x));
                    ca_vec.y.push(Some(atom.y));
                    ca_vec.z.push(Some(atom.z));
                    cb_vec.x.push(None);
                    cb_vec.y.push(None);
                    cb_vec.z.push(None);
                    prev_res_serial = Some(atom.res_serial);
                }
            } else if &atom.atom_name == b" CB " {
                if prev_res_serial == Some(atom.res_serial) {
                    cb_vec.x.last_mut().unwrap().replace(atom.x);
                    cb_vec.y.last_mut().unwrap().replace(atom.y);
                    cb_vec.z.last_mut().unwrap().replace(atom.z);
                } else {
                    res_vec.push(atom.res_serial);
                    cb_vec.x.push(Some(atom.x));
                    cb_vec.y.push(Some(atom.y));
                    cb_vec.z.push(Some(atom.z));
                    ca_vec.x.push(None);
                    ca_vec.y.push(None);
                    ca_vec.z.push(None);
                    prev_res_serial = Some(atom.res_serial);
                }
            }
        }

        CompactStructure {
            num_chains: origin.num_chains,
            chains: origin.chains.clone(),
            num_residues: origin.num_residues,
            residues: res_vec,
            ca_vector: ca_vec,
            cb_vector: cb_vec,
        }
    }

    pub fn get_ca(&self, idx: usize) -> Option<Coordinate> {
        let x =  self.ca_vector.x.get(idx).unwrap_or(&None);
        let y =  self.ca_vector.y.get(idx).unwrap_or(&None);
        let z =  self.ca_vector.z.get(idx).unwrap_or(&None);
        if x.is_some() && y.is_some() && z.is_some() {
            Some(Coordinate::build(x,y,z))
        } else { None }
    }

    pub fn get_cb(&self, idx: usize) -> Option<Coordinate> {
        let x =  self.cb_vector.x.get(idx).unwrap_or(&None);
        let y =  self.cb_vector.y.get(idx).unwrap_or(&None);
        let z =  self.cb_vector.z.get(idx).unwrap_or(&None);
        if x.is_some() && y.is_some() && z.is_some() {
            Some(Coordinate::build(x,y,z))
        } else { None }
    }

    pub fn get_distance(&self, idx1: usize, idx2:usize) -> Option<f32> {
        let ca1 = self.get_ca(idx1) ;
        let ca2 = self.get_ca(idx2);
        if ca1.is_some() && ca2.is_some() {
            let dist = ca1.unwrap().calc_distance(&ca2.unwrap());
            Some(dist)
        } else {None}
    }
    pub fn get_angle(&self, idx1: usize, idx2:usize) -> Option<f32> {
        let ca1 = self.get_ca(idx1) ;
        let cb1 = self.get_cb(idx1) ;
        let ca2 = self.get_ca(idx2) ;
        let cb2 = self.get_cb(idx2) ;
        if ca1.is_some() && cb1.is_some() && ca2.is_some() &&  cb2.is_some() {
            // let angle = ca1.unwrap().calc_dihedral(&ca2.unwrap(), &cb1.unwrap(), &cb2.unwrap());
            let angle = ca1.unwrap().calc_angle(&cb1.unwrap(), &ca2.unwrap(), &cb2.unwrap());
            Some(angle)
        } else {None}
    }
}
