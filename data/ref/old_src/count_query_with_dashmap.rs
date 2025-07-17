
// Fallback to DashMap approach for sparse data
pub fn _count_query_idmode<'a>(
    queries: &[GeometricHash],
    query_map: &HashMap<GeometricHash, ((usize, usize), bool)>,
    offset_table: &SimpleHashMap,
    value_vec: &[u16],
    lookup: &'a [(String, usize, usize, f32, usize)],
    sampling_ratio: Option<f32>,
    sampling_count: Option<usize>,
    freq_filter: Option<f32>,
    length_penalty_power: Option<f32>,
) -> Vec<(usize, StructureResult<'a>)> {
    let sampled = sample_query_idmode(queries, offset_table, sampling_ratio, sampling_count);
    let total_hits = lookup.len() as f32;
    let lp = length_penalty_power.unwrap_or(0.5);
    let query_count_map = DashMap::new();
    sampled.par_iter().for_each(|query| {
        if let Some((off, hc)) = offset_table.get(query) {
            if freq_filter.map_or(false, |f| (*hc as f32 / total_hits) > f) {
                return;
            }
            
            let edge = query_map.get(query).unwrap().0;

            for &val in get_values_with_offset_u16(value_vec, *off, *hc).iter() {
                let (ref id, nid, nres, plddt, db_key) = &lookup[val as usize];
                let idf = (total_hits / (*hc as f32)).log2();

                let mut is_new = false;
                let entry = query_count_map.entry(*nid);
                let mut ref_mut = entry.or_insert_with(|| {
                    is_new = true;
                    StructureResult::new(id, *nid, 1, 2, 1, idf, *nres, *plddt, &edge, *db_key)
                });
                
                let result = ref_mut.value_mut();
                if !is_new {
                    // Use temporary HashSets for deduplication but don't store them
                    let mut temp_nodes = HashSet::new();
                    temp_nodes.insert(edge.0);
                    temp_nodes.insert(edge.1);
                    result.node_count += temp_nodes.len();
                    
                    let mut temp_edges = HashSet::new();
                    temp_edges.insert(edge);
                    result.edge_count += temp_edges.len();
                    
                    result.total_match_count += 1;
                    result.idf += idf;
                }
            }
        }
    });
    
    // Apply length penalty and convert to Vec
    query_count_map.into_par_iter().map(|(nid, mut sr)| {
        sr.idf *= (sr.nres as f32).powf(-lp);
        (nid, sr)
    }).collect()
}
