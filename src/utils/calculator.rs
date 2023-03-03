pub trait Calculate {
    // self(f32,f32,f32)
    fn calc_distance(&self, other:&Self) -> f32;
    fn calc_angle(&self, atom2:&Self, atom3: &Self, atom4: &Self) -> f32;
}

use crate::structure::atom::Coordinate;
