// File: query.rs
// Created: 2023-12-22 17:00:50
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2024 Hyunbin Kim, All rights reserved

use rustc_hash::FxHashMap as HashMap;
use crate::geometry::core::{GeometricHash, HashType};
use crate::index::indextable::FolddiscoIndex;
use crate::structure::chain_id::{split_chain_and_rest, ChainId};
use crate::utils::convert::{is_aa_group_char, map_one_letter_to_u8_vec};
use crate::utils::combination::CombinationIterator;
use crate::utils::log::{log_msg, print_log_msg, FAIL, WARN};
use crate::utils::convert::map_aa_to_u8;
use super::expand::{for_each_expanded_feature, FeatureExpander, ToleranceConfig};
use super::feature::get_single_feature;
use super::io::read_compact_structure;
use super::substitution::{index_matching_substitution, resolve_substitution, substitution_variants, SubstitutionScheme, SCHEME_MARKER};
use crate::structure::core::CompactStructure;

/// Ceiling on hashes generated for one residue pair. Substitutions multiply with the
/// geometric neighbourhood; enumeration is nearest-first, so the cut drops the farthest.
pub const MAX_HASHES_PER_PAIR: usize = 4096;

/// Score factor per substituted residue in a hash, so exact matches rank first.
/// Chosen in docs/feature_evaluation.md §14.
pub const SUBSTITUTION_WEIGHT: f32 = 0.75;

/// What a query hash stands for.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QueryHash {
    /// Query residue indices of the pair.
    pub pair: (usize, usize),
    /// 1 for the query's own amino acids, `SUBSTITUTION_WEIGHT` per substituted side.
    pub weight: f32,
    /// The observed hash's IDF; a substituted hash keeps its own if that is lower.
    pub idf: f32,
    /// Observed geometry, not a neighbouring bin.
    pub observed: bool,
}

impl QueryHash {
    /// Exact residues at the observed geometry.
    pub fn is_exact(&self) -> bool {
        self.weight >= 1.0 && self.observed
    }
}

pub type QueryHashMap = HashMap<GeometricHash, QueryHash>;

/// IDF of `hash` in `index`: log2(total structures / structures with the hash); 0 if absent.
pub fn calculate_idf_for_hash(
    hash: &GeometricHash,
    index: &Option<&FolddiscoIndex>,
    total_structures: f32,
) -> f32 {
    if let Some(ref idx) = index {
        let entries = idx.get_entries(hash.as_u32());
        let hash_count = entries.len();
        if hash_count > 0 {
            return (total_structures / (hash_count as f32)).log2();
        }
    }
    0.0
}

pub fn parse_threshold_string(threshold_string: Option<String>) -> Vec<f32> {
    if threshold_string.is_none() {
        return Vec::new();
    }
    let threshold_string = threshold_string.unwrap();
    // Remove whitespace
    let threshold_string = threshold_string.replace(" ", "");
    let mut thresholds: Vec<f32> = Vec::new();
    for threshold in threshold_string.split(',') {
        let threshold = threshold.parse::<f32>().expect(
            &log_msg(FAIL, "Failed to parse threshold")
        );
        thresholds.push(threshold);
    }
    thresholds
}

/// Insert `feature`'s hash(es). An existing entry is kept unless the new one has a higher
/// weight, or the same weight at observed geometry, so an exact hash is never shadowed by
/// another pair's substitution or neighbouring bin.
fn insert_binned_hash(
    hash_collection: &mut QueryHashMap,
    feature: &Vec<f32>, hash_type: HashType,
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
    value: QueryHash,
) {
    if let Some(multiple_bin) = multiple_bin {
        for (nbin_dist, nbin_angle) in multiple_bin.iter() {
            let hash_value = if *nbin_dist == 0 || *nbin_angle == 0 {
                GeometricHash::perfect_hash_default(feature, hash_type)
            } else {
                GeometricHash::perfect_hash(feature, hash_type, *nbin_dist, *nbin_angle)
            };
            insert_if_heavier(hash_collection, hash_value, value);
        }
    } else {
        let hash_value = if nbin_dist == 0 || nbin_angle == 0 {
            GeometricHash::perfect_hash_default(feature, hash_type)
        } else {
            GeometricHash::perfect_hash(feature, hash_type, nbin_dist, nbin_angle)
        };
        insert_if_heavier(hash_collection, hash_value, value);
    }
}

fn insert_if_heavier(hash_collection: &mut QueryHashMap, hash: GeometricHash, value: QueryHash) {
    match hash_collection.get(&hash) {
        Some(old) if (old.weight, old.observed) >= (value.weight, value.observed) => {}
        _ => { hash_collection.insert(hash, value); }
    }
}

/// Resolve `:*` markers, and with `apply_to_all` every residue lacking an explicit
/// `:ALT`, into `scheme`'s alternatives for the residue observed in `compact`.
/// Residues missing from the structure are left unchanged.
pub fn resolve_query_substitutions(
    compact: &CompactStructure, query_residues: &[(ChainId, u64)],
    substitutions: &[Option<Vec<u8>>], scheme: SubstitutionScheme, apply_to_all: bool,
    serial_query: bool,
) -> Vec<Option<Vec<u8>>> {
    map_observed_residues(compact, query_residues, substitutions, serial_query, |substitution, observed| {
        resolve_substitution(substitution, observed, scheme, apply_to_all)
    })
}

/// Substitutions residue matching needs on an index expanded with `scheme`; see
/// `index_matching_substitution`. `substitutions` should already be resolved.
pub fn index_matching_substitutions(
    compact: &CompactStructure, query_residues: &[(ChainId, u64)],
    substitutions: &[Option<Vec<u8>>], scheme: SubstitutionScheme, serial_query: bool,
) -> Vec<Option<Vec<u8>>> {
    map_observed_residues(compact, query_residues, substitutions, serial_query, |substitution, observed| {
        index_matching_substitution(substitution, observed, scheme)
    })
}

/// Apply `f` to each residue's substitution list and observed amino acid code.
/// Residues missing from the structure keep their list unchanged.
fn map_observed_residues<F>(
    compact: &CompactStructure, query_residues: &[(ChainId, u64)],
    substitutions: &[Option<Vec<u8>>], serial_query: bool, f: F,
) -> Vec<Option<Vec<u8>>>
where
    F: Fn(&Option<Vec<u8>>, u8) -> Option<Vec<u8>>,
{
    query_residues.iter().zip(substitutions.iter()).map(|((chain, ri), substitution)| {
        let index = if serial_query {
            Some(*ri as usize).filter(|&i| i < compact.num_residues)
        } else {
            compact.get_index(chain, ri)
        };
        match index {
            Some(index) => f(substitution, map_aa_to_u8(compact.get_res_name(index))),
            None => substitution.clone(),
        }
    }).collect()
}

/// Build the query hash map for `query_residues` of the structure at `path`.
///
/// Returns (`QueryHashMap`, residue indices, and amino acid pair ->
/// [(Ca distance, first residue index)]), the last used by matching.
pub fn make_query_map(
    path: &String, query_residues: &Vec<(ChainId, u64)>, hash_type: HashType, 
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
    tolerance: &ToleranceConfig,
    amino_acid_substitutions: &Vec<Option<Vec<u8>>>, distance_cutoff: f32, serial_query: bool,
    index: &Option<&FolddiscoIndex>,
    total_structures: f32,
) -> (QueryHashMap, Vec<usize>, HashMap<(u8, u8), Vec<(f32, usize)>>) {
    let (compact, _) = read_compact_structure(path).expect("Failed to read compact structure");
    make_query_map_from_structure(
        &compact, query_residues, hash_type, nbin_dist, nbin_angle, multiple_bin,
        tolerance, amino_acid_substitutions, distance_cutoff, serial_query,
        index, total_structures,
    )
}

/// Same as `make_query_map` for a structure already in memory.
fn make_query_map_from_structure(
    compact: &CompactStructure, query_residues: &Vec<(ChainId, u64)>, hash_type: HashType,
    nbin_dist: usize, nbin_angle: usize, multiple_bin: &Option<Vec<(usize, usize)>>,
    tolerance: &ToleranceConfig,
    amino_acid_substitutions: &Vec<Option<Vec<u8>>>, distance_cutoff: f32, serial_query: bool,
    index: &Option<&FolddiscoIndex>,
    total_structures: f32,
) -> (QueryHashMap, Vec<usize>, HashMap<(u8, u8), Vec<(f32, usize)>>) {
    let mut hash_collection = HashMap::default();
    let mut observed_distance_map: HashMap<(u8, u8), Vec<(f32, usize)>> = HashMap::default();
    
    // Convert residue indices to vector indices
    let mut indices = Vec::new();
    let mut query_residues = query_residues.clone();
    let mut amino_acid_substitutions = amino_acid_substitutions.clone();

    if query_residues.is_empty() {
        // Iterate over all residues and set to query_residues
        for i in 0..compact.num_residues {
            let chain = compact.chain_per_residue[i];
            let residue_index = compact.residue_serial[i];
            query_residues.push((chain, residue_index));
            amino_acid_substitutions.push(None);
        }
    }

    let mut substitution_map: HashMap<usize, Vec<u8>> = HashMap::default();
    
    for (i, (chain, ri)) in query_residues.iter().enumerate() {
        let index = if serial_query { Some(*ri as usize) } else { compact.get_index(&chain, &ri) };
        if let Some(index) = index {
            let _residue: String = compact.get_res_name(index).iter().map(|&c| c as char).collect();
            indices.push(index);
            if let Some(substitution) = amino_acid_substitutions[i].clone() {
                substitution_map.insert(index, substitution);
            }
        }
    }
    // Amino acid positions inside the feature vector, if this hash type encodes them
    let aa_indices = hash_type.amino_acid_index().map(|idx| (idx[0], idx[1]));
    let mut expander = FeatureExpander::new(hash_type, nbin_dist, nbin_angle, tolerance);
    let comb_iter = CombinationIterator::new(indices.len());
    let mut feature = vec![0.0; 9];
    let mut variant = vec![0.0; 9];
    comb_iter.for_each(|(i, j)| {
        if i == j {
            return;
        }
        let is_feature = get_single_feature(
            indices[i], indices[j], &compact, hash_type, distance_cutoff, &mut feature
        );

        if is_feature {
            let observed_hash = if nbin_dist == 0 || nbin_angle == 0 {
                GeometricHash::perfect_hash_default(&feature, hash_type)
            } else {
                GeometricHash::perfect_hash(&feature, hash_type, nbin_dist, nbin_angle)
            };
            let idf = calculate_idf_for_hash(&observed_hash, index, total_structures);
            let pair = (indices[i], indices[j]);
            insert_binned_hash(
                &mut hash_collection, &feature, hash_type, nbin_dist, nbin_angle, multiple_bin,
                QueryHash { pair, weight: 1.0, idf, observed: true },
            );

            // Substitutions apply to every geometric neighbour too, so `164:H` is not
            // limited to the bins the observed residue fell into.
            let aa_variants = match aa_indices {
                Some((a0, a1)) => substitution_variants(
                    (feature[a0], feature[a1]),
                    substitution_map.get(&indices[i]).map(|v| v.as_slice()),
                    substitution_map.get(&indices[j]).map(|v| v.as_slice()),
                ),
                None => Vec::new(),
            };

            // Matching accepts a target residue pair only if its amino acids appear
            // here, so substituted pairs are registered alongside the observed one.
            if let Some((aa_i, aa_j, ca_dist)) = compact.get_list_amino_acids_and_distances(pair.0, pair.1) {
                observed_distance_map.entry((aa_i, aa_j)).or_default().push((ca_dist, pair.0));
                for &(sub_i, sub_j) in &aa_variants {
                    let entry = observed_distance_map.entry((sub_i as u8, sub_j as u8)).or_default();
                    if !entry.contains(&(ca_dist, pair.0)) {
                        entry.push((ca_dist, pair.0));
                    }
                }
            }

            let observed_aa = aa_indices.map(|(a0, a1)| (feature[a0], feature[a1]));
            for_each_expanded_feature(
                &mut expander, &feature, aa_indices, &aa_variants, MAX_HASHES_PER_PAIR,
                &mut variant, |variant| {
                    let (weight, idf) = match (aa_indices, observed_aa) {
                        (Some((a0, a1)), Some((o0, o1))) if (variant[a0], variant[a1]) != (o0, o1) => {
                            let substituted = (variant[a0] != o0) as i32 + (variant[a1] != o1) as i32;
                            (SUBSTITUTION_WEIGHT.powi(substituted),
                             substituted_idf(variant, hash_type, nbin_dist, nbin_angle, idf, index, total_structures))
                        }
                        _ => (1.0, idf),
                    };
                    let observed = (0..variant.len()).all(|k| {
                        aa_indices.map_or(false, |(a0, a1)| k == a0 || k == a1) || variant[k] == feature[k]
                    });
                    insert_binned_hash(
                        &mut hash_collection, variant, hash_type, nbin_dist, nbin_angle, multiple_bin,
                        QueryHash { pair, weight, idf, observed },
                    )
                },
            );
        }
    });
    (hash_collection, indices, observed_distance_map)
}

/// IDF of a substituted hash, capped by the observed hash's so a rare substitution
/// does not outscore the query's own residues. Uncapped when the observed hash is absent.
fn substituted_idf(
    variant: &Vec<f32>, hash_type: HashType, nbin_dist: usize, nbin_angle: usize,
    observed_idf: f32, index: &Option<&FolddiscoIndex>, total_structures: f32,
) -> f32 {
    let hash = if nbin_dist == 0 || nbin_angle == 0 {
        GeometricHash::perfect_hash_default(variant, hash_type)
    } else {
        GeometricHash::perfect_hash(variant, hash_type, nbin_dist, nbin_angle)
    };
    let own = calculate_idf_for_hash(&hash, index, total_structures);
    if observed_idf > 0.0 { own.min(observed_idf) } else { own }
}

/// Widest residue span one `-q` range may expand to. Above the default 50,000-residue
/// index limit; stops a typo like `A204-2150000000` from allocating ~80 GB.
const MAX_RESIDUE_RANGE_SPAN: u64 = 100_000;

/// Residue number at one end of a range, or a single position.
/// The token may repeat the segment's chain (`F204-F215`), but not name another one.
fn parse_residue_number(token: &str, chain: ChainId, segment: &str) -> Result<u64, String> {
    let digits = match split_chain_and_rest(token) {
        (Some(token_chain), rest) if token_chain == chain => rest,
        (Some(token_chain), _) => return Err(format!(
            "Query '{}' names chain '{}' and chain '{}' in '{}'. A range stays in one chain, and a \
             multi-character or numeric chain needs the separator, as in 'AA_250'",
            segment, chain, token_chain, token
        )),
        (None, rest) => rest,
    };
    digits.parse::<u64>().map_err(|_| format!(
        "Query '{}' has '{}' where a residue number was expected", segment, token
    ))
}

/// Parse a query string for the CLI; **exits the process** on malformed input.
/// Library callers should use `parse_query_string_checked`.
pub fn parse_query_string(query_string: &str, default_chain: ChainId) -> (Vec<(ChainId, u64)>, Vec<Option<Vec<u8>>>) {
    match parse_query_string_checked(query_string, default_chain) {
        Ok((query_residues, amino_acid_substitutions)) => {
            warn_on_duplicate_residues(query_string, &query_residues);
            (query_residues, amino_acid_substitutions)
        }
        Err(err) => {
            print_log_msg(FAIL, &err);
            std::process::exit(1);
        }
    }
}

/// Warn when a query names a residue more than once (`B57,B57`, `A1-A5,A3-A7`).
/// Repeats are kept because the list length is the denominator of the coverage ratios.
fn warn_on_duplicate_residues(query_string: &str, query_residues: &[(ChainId, u64)]) {
    let mut distinct = query_residues.to_vec();
    distinct.sort_unstable();
    distinct.dedup();
    if distinct.len() < query_residues.len() {
        print_log_msg(WARN, &format!(
            "Query '{}' names {} residues but only {} distinct ones. The repeats still \
             count toward the query length, so coverage ratios and the novelty verdict are \
             computed against {} and read low",
            query_string, query_residues.len(), distinct.len(), query_residues.len()
        ));
    }
}

/// Parse a `-q` string into `(chain, residue)` pairs and per-residue substitutions,
/// returning an error message for malformed input.
///
/// Accepts `A250`, `A_250`, `AA_250`, `10_250` (separator required for multi-character or
/// numeric chains) and bare `250` (`default_chain`). Ranges (`A250-252`) and
/// substitutions (`A250:R`, `A250:*`) work with every spelling.
pub fn parse_query_string_checked(
    query_string: &str, default_chain: ChainId
)-> Result<(Vec<(ChainId, u64)>, Vec<Option<Vec<u8>>>) , String> {
    let mut query_residues = Vec::new();
    let mut amino_acid_substitutions = Vec::new();

    if query_string.is_empty() {
        return Ok((query_residues, amino_acid_substitutions));
    }
    // A blank or non-alphanumeric single-byte default falls back to chain A
    let default_chain = if default_chain.is_empty()
        || (default_chain.len() == 1 && !default_chain.first_byte().is_ascii_alphanumeric())
    {
        ChainId::from_byte(b'A')
    } else {
        default_chain
    };
    // Remove whitespace
    let query_string = query_string.replace(" ", "");
    for segment in query_string.split(',') {
        let (chain, rest) = match split_chain_and_rest(segment) {
            (Some(chain), rest) => (chain, rest),
            (None, rest) => (default_chain, rest),
        };

        let (range_part, subst_part) = match rest.split_once(':') {
            Some((r, s)) => {
                let mut sub_vec = s
                    .chars()
                    .filter(|c| is_aa_group_char(*c))
                    .flat_map(|c| map_one_letter_to_u8_vec(c))
                    .collect::<Vec<_>>();
                // `*` asks for the substitution scheme; resolved once the residue is known
                if s.contains('*') {
                    sub_vec.push(SCHEME_MARKER);
                }
                (r, Some(sub_vec))
            }
            None => (rest, None),
        };

        if let Some((start_str, end_str)) = range_part.split_once('-') {
            let start = parse_residue_number(start_str, chain, segment)?;
            let end = parse_residue_number(end_str, chain, segment)?;
            if end < start {
                return Err(format!(
                    "Query '{}' ends before it starts; a range runs from the lower residue",
                    segment
                ));
            }
            if end - start >= MAX_RESIDUE_RANGE_SPAN {
                return Err(format!(
                    "Query '{}' spans {} residues, past the {} a range may expand to. \
                     A range is expanded one residue at a time, so check for a mistyped \
                     digit before this becomes an out-of-memory kill",
                    segment, (end - start).saturating_add(1), MAX_RESIDUE_RANGE_SPAN
                ));
            }
            for r in start..=end {
                query_residues.push((chain, r));
                amino_acid_substitutions.push(subst_part.clone());
            }
        } else {
            let residue_num = parse_residue_number(range_part, chain, segment)?;
            query_residues.push((chain, residue_num));
            amino_acid_substitutions.push(subst_part);
        }
    }

    Ok((query_residues, amino_acid_substitutions))
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::chain_id::chain;
    
    fn zinc_finger_query_map(
        tolerance: &ToleranceConfig, substitutions: Vec<Option<Vec<u8>>>,
    ) -> QueryHashMap {
        let path = String::from("query/1G2F.pdb");
        let query_residues = vec![(chain("F"), 207), (chain("F"), 212), (chain("F"), 225)];
        let (hash_collection, _indices, _observed_dist_map) = make_query_map(
            &path, &query_residues, HashType::PDBTrRosetta, 16, 4, &None,
            tolerance, &substitutions, 20.0, false, &None, 1000.0
        );
        hash_collection
    }

    #[test]
    fn test_make_query_map() {
        let path= String::from("query/1G2F.pdb");
        let query_residues = vec![
            (ChainId::from_byte(b'F'), 207), (ChainId::from_byte(b'F'), 212),
            (ChainId::from_byte(b'F'), 225)
        ];
        let amino_acid_substitutions = vec![None; query_residues.len()];
        let hash_type = HashType::PDBTrRosetta;
        let (hash_collection, _index_found, _observed_dist_map) = make_query_map(
            &path, &query_residues, hash_type, 16, 4, &None,
            &ToleranceConfig::default_query(), &amino_acid_substitutions, 20.0, false,
            &None, 1000.0
        );
        // Without substitutions every hash carries the query's own residues
        assert!(hash_collection.values().all(|h| h.weight == 1.0));
        // Observed geometry for at most one hash per ordered pair, neighbours for the rest
        let observed = hash_collection.values().filter(|h| h.observed).count();
        assert!(observed > 0 && observed <= 6, "observed hashes: {}", observed);
        assert!(hash_collection.len() > 6);
    }

    #[test]
    fn wider_tolerance_only_adds_hashes() {
        let no_substitution = vec![None; 3];
        let tight = zinc_finger_query_map(
            &ToleranceConfig::new(vec![0.5], vec![5.0], 1), no_substitution.clone()
        );
        let loose = zinc_finger_query_map(
            &ToleranceConfig::new(vec![0.5], vec![5.0], 2), no_substitution
        );
        // A wider radius is a pure addition: every hash of the tight query survives
        for hash in tight.keys() {
            assert!(loose.contains_key(hash), "radius 2 lost a hash of radius 1");
        }
        assert!(loose.len() > tight.len());
    }

    #[test]
    fn substitutions_compose_with_geometric_tolerance() {
        // His variants must appear across the tolerance neighbourhood, not only at the observed geometry
        let tolerance = ToleranceConfig::new(vec![0.5], vec![5.0], 1);
        let plain = zinc_finger_query_map(&tolerance, vec![None; 3]);
        let substituted = zinc_finger_query_map(
            &tolerance, vec![Some(vec![8]), None, None] // 8 = HIS
        );
        for hash in plain.keys() {
            assert!(substituted.contains_key(hash));
        }
        // Pairs with residue 207 gain one His variant per geometric variant
        assert!(
            substituted.len() > plain.len() + 2,
            "substitution did not compose with tolerance: {} vs {}",
            substituted.len(), plain.len()
        );
    }

    #[test]
    fn scheme_marker_resolves_against_the_observed_residue() {
        let (residues, subs) = parse_query_string("F207,F212:*,F225:Q*", chain("A"));
        assert_eq!(subs, vec![None, Some(vec![SCHEME_MARKER]), Some(vec![5, SCHEME_MARKER])]);
        let (compact, _) = read_compact_structure(&String::from("query/1G2F.pdb")).unwrap();
        let scheme = SubstitutionScheme::Group;
        let observed = |k: usize| map_aa_to_u8(compact.get_res_name(compact.get_index(&residues[k].0, &residues[k].1).unwrap()));
        let resolved = resolve_query_substitutions(&compact, &residues, &subs, scheme, false, false);
        assert_eq!(resolved[0], None);
        assert_eq!(resolved[1], Some(scheme.alternatives(observed(1))));
        let mut expected = vec![5];
        expected.extend(scheme.alternatives(observed(2)).into_iter().filter(|&aa| aa != 5));
        assert_eq!(resolved[2], Some(expected));
        // Global application fills residues without an explicit list
        let all = resolve_query_substitutions(&compact, &residues, &subs, scheme, true, false);
        assert_eq!(all[0], Some(scheme.alternatives(observed(0))));
    }

    #[test]
    fn substituted_pairs_are_registered_for_matching() {
        let path = String::from("query/1G2F.pdb");
        let residues = vec![(chain("F"), 207), (chain("F"), 212), (chain("F"), 225)];
        let (_, _, plain) = make_query_map(
            &path, &residues, HashType::PDBTrRosetta, 16, 4, &None,
            &ToleranceConfig::default_query(), &vec![None; 3], 20.0, false, &None, 1000.0,
        );
        let (_, _, substituted) = make_query_map(
            &path, &residues, HashType::PDBTrRosetta, 16, 4, &None,
            &ToleranceConfig::default_query(), &vec![Some(vec![17]), None, None], 20.0, false, &None, 1000.0,
        );
        // Every observed pair survives, and Trp (17) now appears on the first residue's side
        for key in plain.keys() {
            assert!(substituted.contains_key(key));
        }
        assert!(substituted.keys().any(|&(aa_i, _)| aa_i == 17));
        assert!(!plain.keys().any(|&(aa_i, aa_j)| aa_i == 17 || aa_j == 17));
    }

    #[test]
    fn substituted_hashes_are_weighted_and_never_shadow_exact_ones() {
        let tolerance = ToleranceConfig::default_query();
        let plain = zinc_finger_query_map(&tolerance, vec![None; 3]);
        // 207 and 212 are both Cys; substituting both yields one- and two-sided variants
        let substituted = zinc_finger_query_map(
            &tolerance, vec![Some(vec![8]), Some(vec![8]), None] // 8 = HIS
        );
        for (hash, h) in plain.iter() {
            assert_eq!(h.weight, 1.0);
            assert_eq!(substituted[hash].weight, 1.0, "an exact hash lost its weight");
            assert_eq!(substituted[hash].observed, h.observed);
        }
        let weight = |w: f32| substituted.values().filter(|h| h.weight == w).count();
        assert_eq!(weight(1.0), plain.len());
        assert!(weight(SUBSTITUTION_WEIGHT) > 0);
        assert!(weight(SUBSTITUTION_WEIGHT * SUBSTITUTION_WEIGHT) > 0);
        assert_eq!(weight(1.0) + weight(SUBSTITUTION_WEIGHT) + weight(SUBSTITUTION_WEIGHT * SUBSTITUTION_WEIGHT), substituted.len());
    }

    #[test]
    fn unknown_substitution_code_is_ignored() {
        // 255 (unknown one-letter code) would overflow the residue field
        let tolerance = ToleranceConfig::default_query();
        let plain = zinc_finger_query_map(&tolerance, vec![None; 3]);
        let bogus = zinc_finger_query_map(&tolerance, vec![Some(vec![255]), None, None]);
        assert_eq!(plain.len(), bogus.len());
    }

    #[test]
    fn test_parse_query_string() {
        let query_string = "A250,B232,C269";
        let query_residues = parse_query_string(query_string, chain("A"));
        assert_eq!(query_residues, (vec![(chain("A"), 250), (chain("B"), 232), (chain("C"), 269)], vec![None, None, None]));
    }
    #[test]
    fn test_parse_query_string_with_space() {
        let query_string = "A250, A232, A269";
        let query_residues = parse_query_string(query_string, chain("A"));
        assert_eq!(query_residues, (vec![(chain("A"), 250), (chain("A"), 232), (chain("A"), 269)], vec![None, None, None]));
    }
    
    #[test]
    fn test_parse_query_string_with_space_and_no_chain() {
        let query_string = "250, 232, 269";
        let query_residues = parse_query_string(query_string, chain("A"));
        assert_eq!(query_residues, (vec![(chain("A"), 250), (chain("A"), 232), (chain("A"), 269)], vec![None, None, None]));
    }

    #[test]
    fn test_parse_query_string_with_aa_substitution() {
        let query_string = "A250:R,B232:K,C269:QK";
        let query_residues = parse_query_string(query_string, chain("A"));
        // R = 1, K = 11, Q = 5
        assert_eq!(query_residues, (vec![(chain("A"), 250), (chain("B"), 232), (chain("C"), 269)], vec![Some(vec![1]), Some(vec![11]), Some(vec![5, 11])]));
        let query_string = "250:R,232:K,269:QK";
        let query_residues = parse_query_string(query_string, chain("A"));
        // R = 1, K = 11, Q = 5
        assert_eq!(query_residues, (vec![(chain("A"), 250), (chain("A"), 232), (chain("A"), 269)], vec![Some(vec![1]), Some(vec![11]), Some(vec![5, 11])]));
    }
    #[test]
    fn range_end_uses_declared_chain_only() {
        let parsed = parse_query_string_checked("F204-215,F222-232", chain("A")).unwrap();
        assert_eq!(parsed.0.len(), 23);
        assert_eq!(parsed.0[0], (chain("F"), 204));
        assert_eq!(*parsed.0.last().unwrap(), (chain("F"), 232));

        // Repeating the segment's chain on the end is the same range; another chain is not.
        assert_eq!(parse_query_string_checked("F204-F215", chain("A")).unwrap().0.len(), 12);
        assert_eq!(parse_query_string_checked("AA_1-AA_3", chain("A")).unwrap().0.len(), 3);
        let err = parse_query_string_checked("F204-G215", chain("A")).unwrap_err();
        assert!(err.contains("G215"), "{}", err);
    }

    #[test]
    fn an_unbounded_range_is_a_diagnostic_not_an_allocation() {
        // One mistyped digit would ask for 2.15e9 residues
        let err = parse_query_string_checked("A204-2150000000", chain("A")).unwrap_err();
        assert!(err.contains("2150000000") || err.contains("spans"), "{}", err);
        assert!(err.contains("100000"), "the message should name the limit: {}", err);
        // The default 50,000-residue index limit stays below the cap
        let (residues, _) = parse_query_string_checked("A1-50000", chain("A")).unwrap();
        assert_eq!(residues.len(), 50_000);
        // Exactly at the cap is allowed, one past it is not
        assert_eq!(parse_query_string_checked("A1-100000", chain("A")).unwrap().0.len(), 100_000);
        assert!(parse_query_string_checked("A1-100001", chain("A")).is_err());
    }

    #[test]
    fn duplicate_residues_are_kept_but_countable() {
        // Overlapping ranges repeat residues 3-5; repeats are kept (and warned about)
        let (residues, _) = parse_query_string_checked("A1-A5,A3-A7", chain("A")).unwrap();
        assert_eq!(residues.len(), 10);
        let mut distinct = residues.clone();
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(distinct.len(), 7);
        // The plain repeated-residue form behaves the same way
        assert_eq!(parse_query_string_checked("B57,B57,B57", chain("A")).unwrap().0.len(), 3);
    }

    #[test]
    fn malformed_queries_are_diagnosed_not_panicked() {
        // A range cannot span two chains
        let err = parse_query_string_checked("F204-G215", chain("A")).unwrap_err();
        assert!(err.contains("G215"), "{}", err);
        // Reversed ranges are rejected
        let err = parse_query_string_checked("F215-F204", chain("A")).unwrap_err();
        assert!(err.contains("ends before it starts"), "{}", err);
        // Non-numeric tokens name themselves
        let err = parse_query_string_checked("F20x", chain("A")).unwrap_err();
        assert!(err.contains("F20x") && err.contains("'20x'"), "{}", err);
        assert!(parse_query_string_checked("F204-", chain("A")).is_err());
        // The documented good forms keep working
        assert!(parse_query_string_checked("B57,B102,C195", chain("A")).is_ok());
        assert!(parse_query_string_checked("1-10,11:X", chain("A")).is_ok());
        assert!(parse_query_string_checked("164:H,195,221,247:ND", chain("A")).is_ok());
    }

    #[test]
    fn test_parse_query_string_with_range() {
        let query_string = "A250-252,B232-234,C269:Q";
        let query_residues = parse_query_string(query_string, chain("A"));
        assert_eq!(query_residues, (vec![
            (chain("A"), 250), (chain("A"), 251), (chain("A"), 252), 
            (chain("B"), 232), (chain("B"), 233), (chain("B"), 234), 
            (chain("C"), 269),
        ], vec![None, None, None, None, None, None, Some(vec![5])]));
    }

    #[test]
    fn test_parse_query_string_separator_is_optional_for_single_char_chains() {
        // Both spellings parse to the same motif
        assert_eq!(
            parse_query_string("A250,B232,C269", chain("A")),
            parse_query_string("A_250,B_232,C_269", chain("A"))
        );
        assert_eq!(
            parse_query_string("A250-252,B232-234,C269:Q", chain("A")),
            parse_query_string("A_250-252,B_232-234,C_269:Q", chain("A"))
        );
        assert_eq!(
            parse_query_string("A250:R,B232:K", chain("A")),
            parse_query_string("A_250:R,B_232:K", chain("A"))
        );
    }

    #[test]
    fn test_parse_query_string_multi_char_chain() {
        assert_eq!(
            parse_query_string("AA_250,AB_232,AC_269", chain("A")),
            (vec![(chain("AA"), 250), (chain("AB"), 232), (chain("AC"), 269)], vec![None, None, None])
        );
        // Ranges and substitutions too
        assert_eq!(
            parse_query_string("AA_250-252,AB_232:K", chain("A")),
            (
                vec![(chain("AA"), 250), (chain("AA"), 251), (chain("AA"), 252), (chain("AB"), 232)],
                vec![None, None, None, Some(vec![11])]
            )
        );
    }

    #[test]
    fn test_parse_query_string_numeric_chain() {
        // A numeric chain needs `_`; `10_250` is chain 10, not residue 10250
        assert_eq!(
            parse_query_string("10_250,10_252", chain("A")),
            (vec![(chain("10"), 250), (chain("10"), 252)], vec![None, None])
        );
        assert_eq!(
            parse_query_string("1_250", chain("A")),
            (vec![(chain("1"), 250)], vec![None])
        );
    }

    #[test]
    fn test_parse_query_string_default_chain_may_be_multi_char() {
        // A bare residue uses the default chain, which may be multi-character or numeric
        assert_eq!(
            parse_query_string("250,252", chain("AA")),
            (vec![(chain("AA"), 250), (chain("AA"), 252)], vec![None, None])
        );
        assert_eq!(
            parse_query_string("250", chain("10")),
            (vec![(chain("10"), 250)], vec![None])
        );
        // A blank default falls back to chain A
        assert_eq!(
            parse_query_string("250", ChainId::from_byte(b' ')),
            (vec![(chain("A"), 250)], vec![None])
        );
        assert_eq!(
            parse_query_string("250", ChainId::empty()),
            (vec![(chain("A"), 250)], vec![None])
        );
    }
}
