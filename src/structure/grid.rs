// Function to split the whole structure into grids
// Max grid is 6x6x6

use crate::structure::coordinate::Coordinate;
use crate::structure::core::{Structure, CompactStructure};

pub const MAX_GRID_SIZE: u8 = 6;
pub const DEFAULT_GRID_WIDTH: f32 = 50.0;

pub fn get_grid_count(compact: &CompactStructure, grid_width: f32) -> (Coordinate, Coordinate, Coordinate, Coordinate) {
    // Return min, max, bins
    let min = compact.ca_vector.min_coord();
    let max = compact.ca_vector.max_coord();
    let diff = max.sub(&min);
    let mut bins = diff.clone();
    bins.x = (bins.x / grid_width).ceil().min(MAX_GRID_SIZE as f32);
    bins.y = (bins.y / grid_width).ceil().min(MAX_GRID_SIZE as f32);
    bins.z = (bins.z / grid_width).ceil().min(MAX_GRID_SIZE as f32);
    let widths = Coordinate::new(diff.x / bins.x, diff.y / bins.y, diff.z / bins.z);
    (min, max, bins, widths)
}

pub fn get_grid_index_as_tuple(coord: &Coordinate, min: &Coordinate, widths: &Coordinate, bins: &Coordinate) -> (u8, u8, u8) {
    let diff = coord.sub(min);
    let x = (diff.x / widths.x).min(bins.x - 1.0).max(0.0)as u8;
    let y = (diff.y / widths.y).min(bins.y - 1.0).max(0.0)as u8;
    let z = (diff.z / widths.z).min(bins.y - 1.0).max(0.0)as u8;
    (x, y, z)
}

pub fn get_grid_index(coord: &Coordinate, min: &Coordinate, widths: &Coordinate, bins: &Coordinate) -> u8 {
    let (x, y, z) = get_grid_index_as_tuple(coord, min, widths, bins);
    let index: u8 = MAX_GRID_SIZE.pow(2) * x + MAX_GRID_SIZE * y + z;
    index
}

pub fn grid_index_to_tuple(index: u8) -> (u8, u8, u8) {
    let x = index / (MAX_GRID_SIZE * MAX_GRID_SIZE);
    let y = (index - x * MAX_GRID_SIZE * MAX_GRID_SIZE) / MAX_GRID_SIZE;
    let z = index - x * MAX_GRID_SIZE * MAX_GRID_SIZE - y * MAX_GRID_SIZE;
    (x, y, z)
}
pub fn tuple_to_grid_index(tuple: (u8, u8, u8)) -> u8 {
    MAX_GRID_SIZE.pow(2) * tuple.0 + MAX_GRID_SIZE * tuple.1 + tuple.2
}

pub fn get_grid_index_vector_from_compact(compact: &CompactStructure, grid_width: f32) -> Vec<u8> {
    let (min, _, bins, widths) = get_grid_count(&compact, grid_width);
    let mut grid_index_vector: Vec<u8> = vec![0; compact.num_residues as usize];
    for i in 0..compact.num_residues {
        let ca = compact.get_ca(i);
        match ca {
            Some(coord) => {
                grid_index_vector[i] = get_grid_index(&coord, &min, &widths, &bins);
            },
            None => {
                grid_index_vector[i] = 0;
            }
        }
    }
    grid_index_vector
}

pub fn grid_distance_square(index1: u8, index2: u8) -> u8 {
    let (x1, y1, z1) = grid_index_to_tuple(index1);
    let (x2, y2, z2) = grid_index_to_tuple(index2);
    let x = (x1 as i8 - x2 as i8).abs();
    let y = (y1 as i8 - y2 as i8).abs();
    let z = (z1 as i8 - z2 as i8).abs();
    (x.pow(2) + y.pow(2) + z.pow(2)) as u8
}

pub fn grid_near_than_distance(index1: u8, index2: u8, distance: u8) -> bool {
    let dist = grid_distance_square(index1, index2);
    dist <= distance
}

pub fn nearby(index1: u8, index2: u8) -> bool {
    grid_near_than_distance(index1, index2, 3)
}

pub fn merge_id_with_grid(id: usize, grid_index: u8) -> usize {
    id << 8 | grid_index as usize
}

pub fn split_id_grid_from_usize(id_grid: usize) -> (usize, u8) {
    let id = (id_grid >> 8) as usize;
    let grid_index = (id_grid & 0xFF) as u8;
    (id, grid_index)
}

pub fn split_id_grid_from_bytes(bytes: &[u8]) -> (usize, u8) {
    let id = (bytes[0] as usize) << 8 | (bytes[1] as usize);
    let grid_index = bytes[2];
    (id, grid_index)
}

pub fn convert_to_id_grid_vector(bytes: &[u8]) -> Vec<(usize, u8)> {
    let mut id_grid_vector = Vec::new();
    for i in 0..bytes.len() / 3 {
        let (id, grid_index) = split_id_grid_from_bytes(&bytes[i * 3..i * 3 + 3]);
        id_grid_vector.push((id, grid_index));
    }
    id_grid_vector
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::io::pdb::Reader;
    use std::fs::File;
    use std::path::Path;
    #[test]
    fn test_get_grid_index() {
        let path = Path::new("data/homeobox/1akha-.pdb");
        // let path = Path::new("analysis/h_sapiens_pdb/AF-Q02817-F3-model_v4.pdb");
        // let path = Path::new("analysis/query/4cha.pdb");
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure().unwrap();
        let compact = structure.to_compact();
        println!("{:?}", compact.num_residues);
        let (min, max, bins, widths) = get_grid_count(&compact, 20.0);
        println!("Min: {:?}, Max: {:?}, Bins: {:?}, Widths: {:?}", min, max, bins, widths);
        let grid_index_vector = get_grid_index_vector_from_compact(&compact, 20.0);
        // let grid_index_tuple_vector: Vec<(u8, u8, u8)> = grid_index_vector.iter().map(|&index| grid_index_to_tuple(index)).collect();
        println!("{:?}", grid_index_vector);
        let grid_distance_vector: Vec<u8> = grid_index_vector.iter().map(|&index| grid_distance_square(index, grid_index_vector[150])).collect();
        println!("{:?}", grid_distance_vector);
        // Count the elements in grid_distance_vector less or equal than root3
        let count = grid_distance_vector.iter().filter(|&&dist| dist <= 3).count();
        println!("Count: {}", count);
        // println!("{:?}", grid_index_tuple_vector);
    }
}