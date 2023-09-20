// WARNING: Leave this until the prototyping is done 2023-06-08 12:53:26
use std::collections::{HashSet, HashMap};
use crate::structure::core::CompactStructure;
use crate::structure::coordinate::{approx_cb, CarbonCoordinateVector, Coordinate};

pub struct Grid {
    pub grid_size: usize,
    pub spacing: f32,
    pub grid: HashMap<(usize, usize, usize), HashSet<usize>>,
}

impl Grid {
    pub fn new(structure: &CompactStructure, spacing: f32) -> Grid {
        // Get max & min coordinates for each axis
        let (mut min_x, mut min_y, mut min_z) = (std::f32::MAX, std::f32::MAX, std::f32::MAX);
        let (mut max_x, mut max_y, mut max_z) = (std::f32::MIN, std::f32::MIN, std::f32::MIN);
        for i in 0..structure.num_residues {
            let ca_i = structure.get_ca(i).unwrap();
            if ca_i.x < min_x {
                min_x = ca_i.x;
            }
            if ca_i.y < min_y {
                min_y = ca_i.y;
            }
            if ca_i.z < min_z {
                min_z = ca_i.z;
            }
            if ca_i.x > max_x {
                max_x = ca_i.x;
            }
            if ca_i.y > max_y {
                max_y = ca_i.y;
            }
            if ca_i.z > max_z {
                max_z = ca_i.z;
            }
        }

        let min_coord = Coordinate::new(min_x, min_y, min_z);
        let max_coord = Coordinate::new(max_x, max_y, max_z);
        let max_dist = max_coord.distance(&min_coord);

        // Calculate the grid size
        let grid_size = (max_dist / spacing).ceil() as usize;

        // Initialize the grid
        let mut grid = HashMap::new();
        // Iterate over each residue in the structure
        for i in 0..structure.num_residues {
            // Get the coordinates of the residue's CA atom
            let ca_i = structure.get_ca(i);
            // Calculate the grid indices of the residue's CA atom
            if let Some(ca_i) = ca_i {
                let (x_i, y_i, z_i) = (
                    (ca_i.x / spacing) as usize,
                    (ca_i.y / spacing) as usize,
                    (ca_i.z / spacing) as usize,
                );
                // Add the residue index to the grid cell
                let cell = grid.entry((x_i, y_i, z_i)).or_insert(HashSet::new());
                cell.insert(i);
            }
        }

        Grid {
            grid_size,
            spacing,
            grid,
        }
    }
    
    pub fn calc_pairwise_dist(&self, structure: &CompactStructure) -> Vec<f32> {
        let mut dist_vec = Vec::new();
        for i in 0..structure.num_residues {
            let ca_i = structure.get_ca(i);
            if let Some(ca_i) = ca_i {
                for j in 0..structure.num_residues {
                    let ca_j = structure.get_ca(j);
                    if let Some(ca_j) = ca_j {
                        let dist = ca_i.distance(&ca_j);
                        dist_vec.push(dist);
                    }
                }
            }
        }
        dist_vec
    }
    
}

// Test
#[cfg(test)]
mod grid_tests {
    use super::*;
    use crate::PDBReader;
    use crate::structure::core::CompactStructure;
    use crate::structure::coordinate::CarbonCoordinateVector;
    
    fn pairwise_distance(structure: &CompactStructure) -> Vec<f32> {
        let mut dist_vec = Vec::new();
        for i in 0..structure.num_residues {
            let ca_i = structure.get_ca(i);
            if let Some(ca_i) = ca_i {
                for j in 0..structure.num_residues {
                    let ca_j = structure.get_ca(j);
                    if let Some(ca_j) = ca_j {
                        let dist = ca_i.distance(&ca_j);
                        dist_vec.push(dist);
                    }
                }
            }
        }
        dist_vec
    }

    #[test]
    fn test_grid() {
        let path = String::from("data/serine_peptidases_filtered/1qfm.pdb");
        let reader = PDBReader::from_file(&path).expect("Failed to read PDB file");
        let structure: CompactStructure = reader.read_structure().expect("Failed to convert").to_compact();
        let grid = Grid::new(&structure, 20.0);
        let start_time = std::time::Instant::now();
        let dist_vec = grid.calc_pairwise_dist(&structure);
        let end_time = std::time::Instant::now();
        println!("Time elapsed - GRID: {:?}", end_time - start_time);
        let start_time = std::time::Instant::now();
        let dist_vec2 = pairwise_distance(&structure);
        let end_time = std::time::Instant::now();
        println!("Time elapsed - NO GRID: {:?}", end_time - start_time);
    }
    
}