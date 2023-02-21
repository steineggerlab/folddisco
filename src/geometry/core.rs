

#[derive(Debug)]
pub struct GeometricHasher {
    angles: Vec<f64>,
    distances: Vec<f64>,
    hash: Vec<u64>, // TODO: Figure out the data type for hash
}

impl GeometricHasher {
    pub fn new(angles: Vec<f64>, distances: Vec<f64>) -> GeometricHasher {
        GeometricHasher {
            angles: angles,
            distances: distances,
            hash: Vec::new(),
        }
    }
}

