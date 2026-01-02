pub fn parse_distance_angle_pairs(input: &str) -> Vec<(usize, usize)> {
    input
        .split(',')
        .filter_map(|pair_str| {
            let parts: Vec<_> = pair_str.trim().split('-').collect();
            if parts.len() == 2 {
                match (parts[0].parse(), parts[1].parse()) {
                    (Ok(a), Ok(b)) => Some((a, b)),
                    _ => None,
                }
            } else {
                None
            }
        })
        .collect()
}