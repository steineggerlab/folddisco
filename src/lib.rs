//! # About project
//!
//! Motifsearch is a tool for finding discontinuous motifs in protein structures.

pub mod geometry;
pub mod index;
pub mod structure;
pub mod utils;
pub mod controller;

pub struct MotifSearch {
    pub pdb_files: Vec<String>,
    // pub pdb_table: HashMap<index, PDB>,
    pub controller: controller::Controller,
}





pub fn run() {
    let f = structure::io::pdb::Reader::from_file("data/111l_alpha.pdb").unwrap();
    let structure = &f.read_structure().unwrap();
    let compact = &structure.to_compact();

    let _a1_ca = compact.get_ca(161).unwrap();
    let _a1_cb = compact.get_ca(5).unwrap();

    let dist = compact.get_distance(161, 5).unwrap();
    let angle = compact.get_angle(161,5).unwrap();
    let hashvalue = geometry::hash::HashValue::perfect_hash(dist, angle);
    println!("original dist:{}, angle:{}", dist, angle);
    println!("{:?}", hashvalue);
    // println!("reverse hash(dist,angle): {:?}", hashvalue.reverse_hash());
}