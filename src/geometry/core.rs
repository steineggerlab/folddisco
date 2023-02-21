

#[derive(Debug)]
pub struct GeometricHasher {
    angles: Vec<f32>,
    distances: Vec<f32>,
    hash: Vec<u64>, // TODO: Figure out the data type for hash
}

impl GeometricHasher {
    pub fn new(angles: Vec<f32>, distances: Vec<f32>) -> GeometricHasher {
        GeometricHasher {
            angles: angles,
            distances: distances,
            hash: Vec::new(),
        }
    }
}

