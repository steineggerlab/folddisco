#[derive(Debug)]
pub struct HashValue {
    // revise HashValue
    // TODO: Figure out the data type for hash
    dist: f32,
    angle: f32,
}

impl HashValue {
    pub fn hash(dist: &f32, angle: &f32) -> Self {
        HashValue { dist: *dist, angle: *angle }
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