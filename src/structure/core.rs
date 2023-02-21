use super::atom::AtomVector;

/// Structure is the main data structure for storing the information of a protein structure.
pub struct Structure {
    pub num_chains: usize,
    pub chains: Vec<u8>,
    pub atom_vectors: Vec<AtomVector>,
    pub num_atoms: usize,
    pub num_residues: usize,
}
