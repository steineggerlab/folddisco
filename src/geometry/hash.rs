use std::fmt;
// use std::hash::Hasher;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u16);

impl HashValue {
    // Constructor
    // pub fn new(dist: u16, angle: u16) -> Self {
    //     let hashvalue = dist << 8 | angle ;
    //     HashValue(hashvalue)
    // }

    pub fn from_u16(hashvalue: u16) -> Self {
        HashValue(hashvalue)
    }

    pub fn perfect_hash(dist: f32, angle: f32) -> Self {
        let dist = dist.round() as u16;
        let angle = angle.round() as u16;
        let hashvalue = dist << 8 | angle;
        HashValue(hashvalue)
    }
    pub fn reverse_hash(&self) -> (u16, u16) {
        let dist = self.0 >> 8;
        let angle = self.0 & 0b11111111;
        (dist, angle)
    }
    pub fn loose_hash(&self) -> Self {
        let (dist, angle) = self.reverse_hash();
        let dist = dist / 2; // 20 seems to be too loose
        let angle = angle / 5; // 20 seems to be too loose
        let hashvalue = dist << 8 | angle;
        HashValue(hashvalue)
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

#[cfg(test)]
mod tests {
    // use super::*;
    use crate::{geometry, structure};
    use crate::utils::calculator::Calculate;
    #[test]
    fn test_hash1() {
        let f = structure::io::pdb::Reader::from_file("data/111l_alpha.pdb").unwrap();
        let structure = &f.read_structure().unwrap();
        let compact = &structure.to_compact();
        
        for i in 0..compact.num_residues {
            for j in i+1..compact.num_residues {
                let ci = compact.get_ca(i).expect("compact failed to get CA");
                let cj = compact.get_ca(j).expect("compact failed to get CA");

                let dist = compact.get_distance(i,j).expect("compact failed to get distance");
                let angle = compact.get_angle(i,j).expect("compact failed to get angle");
                let torsion = ci.calc_torsion_angle(&cj, &cj, &cj);

                let hashvalue = geometry::hash::HashValue::perfect_hash(dist, angle);
                let reverse = hashvalue.reverse_hash();

                // println!("residue1: {}, residue2: {}, dist: {}=={}, angle: {}=={}", i, j, dist.round(), reverse.0, angle.round(), reverse.1);
                assert_eq!(dist.round() as u16, reverse.0);
                assert_eq!(angle.round() as u16, reverse.1);
            }
        }
    }

    #[test]
    fn test_hash2() {
        let f = structure::io::pdb::Reader::from_file("data/AF-Q8W3K0-F1-model_v4.pdb").unwrap();
        let structure = &f.read_structure().unwrap();
        let compact = &structure.to_compact();
        
        for i in 0..compact.num_residues {
            for j in i+1..compact.num_residues {
                let ci = compact.get_ca(i).expect("compact failed to get CA");
                let cj = compact.get_ca(j).expect("compact failed to get CA");

                let dist = compact.get_distance(i,j).expect("compact failed to get distance");
                let angle = compact.get_angle(i,j).expect("compact failed to get angle");

                let hashvalue = geometry::hash::HashValue::perfect_hash(dist, angle);
                let reverse = hashvalue.reverse_hash();

                println!("residue1: {}, residue2: {}, dist: {}=={}, angle: {}=={}", i, j, dist.round(), reverse.0, angle.round(), reverse.1);
                assert_eq!(dist.round() as u16, reverse.0);
                assert_eq!(angle.round() as u16, reverse.1);
            }
        } 
    }

}
