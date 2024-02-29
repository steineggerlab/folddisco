use super::atom;

#[derive(Debug)]
pub struct Residue {
    pub chain: u8,
    pub res_serial: u64,
    pub res_name: [u8; 3],
}

#[derive(Debut,Clone)]
pub struct ResidueVector {

}