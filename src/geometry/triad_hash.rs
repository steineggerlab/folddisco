use std::fmt;

#[derive(Ord, PartialOrd, Eq, PartialEq, Clone, Copy, Hash)]
pub struct HashValue(u16);

impl HashValue {

    pub fn perfect_hash(edge1: f32, edge2: f32, edge3: f32) -> Self {
        let mut hashvalue = 0;
        let mut edges = vec![edge1, edge2, edge3];
        // If edges are not sorted, panic
        if edges[0] > edges[1] {
            panic!("01 Edges are not sorted");
        } else if edges[1] > edges[2] {
            panic!("02 Edges are not sorted");
        }
        // If edges are not in range, panic
        if edges[0] < 3.5 || edges[0] > 19.5 {
            panic!("03 Edges are not in range");
        } else if edges[1] < 3.5 || edges[1] > 19.5 {
            panic!("04 Edges are not in range");
        } else if edges[2] < 3.5 || edges[2] > 19.5 {
            println!("{}", edges[2]);
            panic!("05 Edges are not in range");
        }

        let edge1 = (edge1.round() - 4.0_f32) as u16;
        // Assert that higher bits are empty
        let edge2 = (edge2.round() - 4.0_f32) as u16;
        let edge3 = (edge3.round() - 4.0_f32) as u16;
        hashvalue += edge1 << 8;
        hashvalue += edge2 << 4;
        hashvalue += edge3;
        HashValue(hashvalue)
    }

    pub fn reverse_hash(&self) -> (u16, u16, u16) {
        let edge1 = ((self.0 >> 8) & 0b1111) + 4;
        let edge2 = ((self.0 >> 4) & 0b1111) + 4;
        let edge3 = (self.0 & 0b1111) + 4;
        if edge1 < 4 || edge1 > 19 {
            println!("{}", edge1);
            panic!("06 Edges are not in range");
        } else if edge2 < 4 || edge2 > 19 {
            panic!("07 Edges are not in range");
        } else if edge3 < 4 || edge3 > 19 {
            println!("{}", edge3);
            panic!("08 Edges are not in range");
        }
        (edge1, edge2, edge3)
    }

}

#[cfg(test)]
mod temp_test {
    #[test]
    fn test_one_hash() {
        let hash = super::HashValue::perfect_hash(5.0, 8.0, 11.0);
        assert_eq!(hash.0, 0b0000000101000111);
        let reverse = hash.reverse_hash();
        assert_eq!(reverse, (5, 8, 11));
    }
}

impl fmt::Debug for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (edge1, edge2, edge3) = self.reverse_hash();
        write!(f, "HashValue({}), edge1={}, edge2={}, edge3={}", self.0, edge1, edge2, edge3)
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (edge1, edge2, edge3) = self.reverse_hash();
        write!(f, "{}\t{}\t{}\t{}", self.0, edge1, edge2, edge3)
    }
}

pub type HashCollection = Vec<HashValue>;


pub struct PossibleTriangles {
    pub dist_vec: Vec<u8>,
    pub triangles: Vec<(u8, u8, u8)>,
}

impl PossibleTriangles {
    pub fn new(dist_vec: Vec<u8>) -> Self {
        let mut triangles = Vec::new();
        for i in 0..dist_vec.len() {
            for j in i+1..dist_vec.len() {
                for k in j+1..dist_vec.len() {
                    if dist_vec[i] + dist_vec[j] > dist_vec[k] &&
                       dist_vec[i] + dist_vec[k] > dist_vec[j] &&
                       dist_vec[j] + dist_vec[k] > dist_vec[i] {
                        triangles.push((dist_vec[i], dist_vec[j], dist_vec[k]));
                    }
                }
            }
        }
        PossibleTriangles {
            dist_vec,
            triangles,
        }
    }
    pub fn get_triangles(&self) -> &Vec<(u8, u8, u8)> {
        &self.triangles
    }
    pub fn get_dist_vec(&self) -> &Vec<u8> {
        &self.dist_vec
    }

    pub fn all_triangles(&self) -> Vec<(u8, u8, u8)> {
        let mut triangles = Vec::new();
        for i in 0..self.dist_vec.len() {
            for j in i+1..self.dist_vec.len() {
                for k in j+1..self.dist_vec.len() {
                    if self.dist_vec[i] + self.dist_vec[j] > self.dist_vec[k] &&
                       self.dist_vec[i] + self.dist_vec[k] > self.dist_vec[j] &&
                       self.dist_vec[j] + self.dist_vec[k] > self.dist_vec[i] {
                        triangles.push((self.dist_vec[i], self.dist_vec[j], self.dist_vec[k]));
                    }
                }
            }
        }
        triangles
    }

}