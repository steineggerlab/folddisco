use std::fmt;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u16);

impl HashValue {
    // Constructor
    // pub fn new(dist: u16, angle: u16) -> Self {
    //     let hashvalue = dist << 8 | angle ;
    //     HashValue(hashvalue)
    // }

    pub fn perfect_hash(dist: f32, angle: f32) -> Self {
        let dist = dist.round() as u16;
        let angle = angle.round() as u16;
        let hashvalue = dist << 8 | angle ;
        HashValue(hashvalue)
    }
    pub fn reverse_hash(&self) -> (u16, u16) {
        let dist = self.0 >> 8;
        let angle = self.0 & 0b11111111;
        (dist, angle)
    }
}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (dist, angle) = self.reverse_hash();
        write!(f, "HashValue({}), dist={}, angle={}", self.0, dist, angle)
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (dist, angle) = self.reverse_hash();
        write!(f, "{}\t{}\t{}", self.0, dist, angle)
        // write!(f, "{}", self.0)
    }
}

pub type HashCollection = Vec<HashValue>;


// #[derive(Debug)]
// pub struct GeometricHasher {
//     angles: Vec<f32>,
//     distances: Vec<f32>,
//     hashvalue: Vec<u64>, // TODO: Figure out the data type for hash
// }

// impl GeometricHasher {
//     pub fn new(angles: Vec<f32>, distances: Vec<f32>) -> GeometricHasher {
//         GeometricHasher {
//             angles: angles,
//             distances: distances,
//             hashvalue: Vec::new(),
//         }
//     }
// }