#[derive(Debug)]
pub struct HashValue(u16);

impl HashValue {
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

#[derive(Debug)]
pub struct HashCollection(pub Vec<HashValue>);

impl HashCollection {
    pub fn new() -> Self{
        HashCollection(Vec::new())
    }

    pub fn push(&mut self, hashvalue : HashValue) {
        self.push(hashvalue);
    }
}

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