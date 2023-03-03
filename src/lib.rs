//! # About project
//!
//! Motifsearch is a tool for finding discontinuous motifs in protein structures.

pub mod geometry;
pub mod index;
pub mod structure;
pub mod utils;
pub mod controller;

pub fn Run() {
    use crate::utils::calculator::Calculate;
    let f = structure::io::pdb::Reader::from_file("data/111l_alpha.pdb").unwrap();
    let structure = &f.read_structure().unwrap();
    let compact = &structure.to_compact();

    let a1_ca = compact.get_CA(161).unwrap();
    let a1_cb = compact.get_CA(5).unwrap();

    let dist = compact.get_distance(161, 5).unwrap();
    let angle = compact.get_angle(161,5).unwrap();
    let hashvalue = geometry::hash::HashValue::perfect_hash(dist, angle);
    println!("dist:{}, angle:{}", dist, angle);
    println!("{:?}", hashvalue);
    println!("reverse hash(dist,angle): {:?}", hashvalue.reverse_hash());
}