// feature.rs 
pub fn get_geometric_hash_with_grid(structure: &CompactStructure, id: usize, hash_type: HashType, nbin_dist: usize, nbin_angle: usize, grid_width: f32) -> Vec<(GeometricHash, usize)> {
    let grid_indices = get_grid_index_vector_from_compact(structure, grid_width);
    let res_bound = CombinationIterator::new(structure.num_residues);
    let mut hash_vec = Vec::new();

    res_bound.for_each(|(i, j)| {
        if !nearby(grid_indices[i], grid_indices[j]) {
            return;
        }
        let feature = get_single_feature(i, j, structure, hash_type);
        if let Some(feature) = feature {
            let hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(feature, hash_type)
            } else {
                GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
            };
            let id_grid = merge_id_with_grid(id, grid_indices[i]);
            hash_vec.push((hash, id_grid));
        }
    });

    // Reduce memory usage
    hash_vec.shrink_to_fit();
    hash_vec
}
