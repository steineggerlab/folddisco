// File: expand.rs
// Author: Hyunbin Kim (khb7840@gmail.com)
// Copyright © 2026 Hyunbin Kim, All rights reserved
//
// Tolerance-driven expansion of a feature vector into neighbouring hash bins.
//
// - Up to `radius` dimensions deviate at once (features often cross boundaries together).
// - Tolerances wider than one bin are sub-stepped so no bin in between is skipped.
// - Perturbed angles wrap (torsions) or reflect (bounded angles); distances are left as is.
//
// Neighbours come in radius order (all 1-dim before 2-dim), nearer offsets first within a
// dimension, so early stopping keeps whole radius levels.

use crate::geometry::core::HashType;
use super::substitution::SubstitutionScheme;

/// Cap on sub-steps per dimension, bounding the cost of a very large threshold.
const MAX_SUBSTEPS: usize = 8;

/// How far a query hash is expanded around the observed geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct ToleranceConfig {
    /// Distance offsets in Angstroms.
    pub dist_thresholds: Vec<f32>,
    /// Angle offsets in degrees, as given on the command line.
    pub angle_thresholds: Vec<f32>,
    /// Feature dimensions allowed to deviate from the observed bin at once.
    pub radius: usize,
}

impl ToleranceConfig {
    pub fn new(
        dist_thresholds: Vec<f32>, angle_thresholds: Vec<f32>, radius: usize,
    ) -> Self {
        Self { dist_thresholds, angle_thresholds, radius }
    }

    /// Tolerances used by `folddisco query` when nothing is given.
    pub fn default_query() -> Self {
        Self::new(vec![0.5], vec![5.0], 1)
    }

    /// No expansion at all: only the observed geometry is queried.
    pub fn none() -> Self {
        Self::new(Vec::new(), Vec::new(), 0)
    }
}

/// Expansion stored into an index at build time: each target pair is also indexed
/// under the hashes a query with this tolerance and substitution scheme would reach.
#[derive(Debug, Clone, PartialEq)]
pub struct IndexExpansion {
    pub tolerance: ToleranceConfig,
    pub scheme: Option<SubstitutionScheme>,
}

impl IndexExpansion {
    /// `None` when neither geometry nor residue identity would be expanded.
    pub fn new(radius: usize, distance: f32, angle: f32, scheme: Option<SubstitutionScheme>) -> Option<Self> {
        if radius == 0 && scheme.is_none() {
            return None;
        }
        Some(Self { tolerance: ToleranceConfig::new(vec![distance], vec![angle], radius), scheme })
    }

    /// Tolerance matching must use so every candidate the expanded index returned for a
    /// query expanded by `query` can still be matched: offsets and radii add up.
    pub fn matching_tolerance(&self, query: &ToleranceConfig) -> ToleranceConfig {
        if self.tolerance.radius == 0 {
            return query.clone();
        }
        let add = |a: &[f32], b: &[f32]| vec![widest_tolerance(a) + widest_tolerance(b)];
        ToleranceConfig::new(
            add(&query.dist_thresholds, &self.tolerance.dist_thresholds),
            add(&query.angle_thresholds, &self.tolerance.angle_thresholds),
            query.radius + self.tolerance.radius,
        )
    }
}

#[derive(Debug, Clone)]
struct TolerantDim {
    /// Position of this dimension inside the feature vector
    index: usize,
    /// Offsets covering the widest tolerance, nearest first.
    offsets: Vec<f32>,
}

/// Enumerates the tolerance neighbourhood of a feature vector for one hash type.
/// Reuses internal buffers, so it does not allocate after construction.
pub struct FeatureExpander {
    hash_type: HashType,
    dims: Vec<TolerantDim>,
    radius: usize,
    /// Feature vector with offsets applied, before domain fix-up
    raw: Vec<f32>,
    /// `raw` pulled back into the valid domain; this is what callers see
    fixed: Vec<f32>,
}

impl FeatureExpander {
    pub fn new(
        hash_type: HashType, nbin_dist: usize, nbin_angle: usize, tolerance: &ToleranceConfig,
    ) -> Self {
        let dist_bin_width = hash_type.dist_bin_width(nbin_dist);
        let angle_bin_width = hash_type.angle_bin_width(nbin_angle);
        let mut dims: Vec<TolerantDim> = Vec::new();

        if let Some(dist_indices) = hash_type.dist_index() {
            let widest = widest_tolerance(&tolerance.dist_thresholds);
            if widest > 0.0 {
                let offsets = substep_offsets(widest, dist_bin_width);
                for index in dist_indices {
                    dims.push(TolerantDim { index, offsets: offsets.clone() });
                }
            }
        }

        if let Some(angle_indices) = hash_type.angle_index() {
            // CLI thresholds are degrees; only PDBMotif stores degrees
            let widest = widest_tolerance(&tolerance.angle_thresholds);
            let widest = if hash_type.angle_in_degrees() { widest } else { widest.to_radians() };
            if widest > 0.0 {
                let offsets = substep_offsets(widest, angle_bin_width);
                for index in angle_indices {
                    dims.push(TolerantDim { index, offsets: offsets.clone() });
                }
            }
        }

        Self {
            hash_type,
            dims,
            radius: tolerance.radius,
            raw: vec![0.0; 9],
            fixed: vec![0.0; 9],
        }
    }

    /// True when this expander would never produce a neighbour.
    pub fn is_noop(&self) -> bool {
        self.radius == 0
            || self.dims.is_empty()
            || self.dims.iter().all(|d| d.offsets.is_empty())
    }

    /// Number of feature dimensions that carry a tolerance.
    pub fn dim_count(&self) -> usize {
        self.dims.len()
    }

    /// Visit every neighbour of `feature` (not `feature` itself), nearest first.
    /// `visit` returns `false` to stop early.
    pub fn for_each_neighbor<F>(&mut self, feature: &[f32], mut visit: F)
    where
        F: FnMut(&[f32]) -> bool,
    {
        if self.is_noop() {
            return;
        }
        let radius = self.radius.min(self.dims.len());
        self.raw[..feature.len()].copy_from_slice(feature);
        for level in 1..=radius {
            if !self.visit_level(level, 0, feature, &mut visit) {
                return;
            }
        }
    }

    /// Choose `level` distinct dimensions from `dims[start..]` and offset each one.
    /// Returns false once `visit` asked to stop.
    fn visit_level<F>(&mut self, level: usize, start: usize, feature: &[f32], visit: &mut F) -> bool
    where
        F: FnMut(&[f32]) -> bool,
    {
        if level == 0 {
            let len = feature.len();
            self.fixed[..len].copy_from_slice(&self.raw[..len]);
            self.hash_type.sanitize_perturbed_feature(&mut self.fixed[..len]);
            return visit(&self.fixed[..len]);
        }
        // Leave room for the remaining levels
        let last = self.dims.len() - level;
        for dim in start..=last {
            let index = self.dims[dim].index;
            let original = self.raw[index];
            for offset_pos in 0..self.dims[dim].offsets.len() {
                self.raw[index] = original + self.dims[dim].offsets[offset_pos];
                if !self.visit_level(level - 1, dim + 1, feature, visit) {
                    self.raw[index] = original;
                    return false;
                }
            }
            self.raw[index] = original;
        }
        true
    }
}

/// Visit the expansion of one residue pair's `feature`: its amino acid variants, then
/// each geometric neighbour with the observed and every variant residue pair.
///
/// The observed feature itself is not visited. Stops before the total, observed hash
/// included, would exceed `cap`. `variant` is a scratch buffer as long as `feature`.
/// Query and index expansion share this so both sides reach the same hashes.
pub fn for_each_expanded_feature<F: FnMut(&Vec<f32>)>(
    expander: &mut FeatureExpander, feature: &Vec<f32>, aa_indices: Option<(usize, usize)>,
    aa_variants: &[(f32, f32)], cap: usize, variant: &mut Vec<f32>, mut visit: F,
) {
    let aa_variants = if aa_indices.is_some() { aa_variants } else { &[] };
    variant.copy_from_slice(feature);
    if let Some((a0, a1)) = aa_indices {
        for &(aa_i, aa_j) in aa_variants {
            variant[a0] = aa_i;
            variant[a1] = aa_j;
            visit(variant);
        }
    }
    let per_neighbor = 1 + aa_variants.len();
    let mut generated = per_neighbor;
    expander.for_each_neighbor(feature, |neighbor| {
        variant[..neighbor.len()].copy_from_slice(neighbor);
        visit(variant);
        if let Some((a0, a1)) = aa_indices {
            for &(aa_i, aa_j) in aa_variants {
                variant[a0] = aa_i;
                variant[a1] = aa_j;
                visit(variant);
            }
        }
        generated += per_neighbor;
        generated + per_neighbor <= cap
    });
}

/// Widest magnitude in a threshold list; smaller ones reach no extra bin.
fn widest_tolerance(thresholds: &[f32]) -> f32 {
    thresholds.iter().fold(0.0f32, |widest, t| widest.max(t.abs()))
}

/// Signed offsets covering `[-tolerance, +tolerance]`, nearest first, at most one bin
/// apart so no bin in between is skipped.
fn substep_offsets(tolerance: f32, bin_width: f32) -> Vec<f32> {
    if !(tolerance > 0.0) {
        return Vec::new();
    }
    let steps = if bin_width.is_finite() && bin_width > 0.0 {
        ((tolerance / bin_width).ceil() as usize).clamp(1, MAX_SUBSTEPS)
    } else {
        1
    };
    let mut offsets = Vec::with_capacity(steps * 2);
    for step in 1..=steps {
        let magnitude = tolerance * step as f32 / steps as f32;
        offsets.push(-magnitude);
        offsets.push(magnitude);
    }
    offsets
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::convert::{reflect_into_range, wrap_to_pi};
    use crate::geometry::core::GeometricHash;
    use rustc_hash::FxHashSet as HashSet;

    const PI: f32 = std::f32::consts::PI;

    /// PDBTrRosetta feature: aa1, aa2, ca_dist, cb_dist, ca_cb_angle, theta1, theta2
    fn pdbtr_feature() -> Vec<f32> {
        let mut feature = vec![0.0; 9];
        feature[0] = 8.0;   // HIS
        feature[1] = 3.0;   // ASP
        feature[2] = 9.4;
        feature[3] = 8.1;
        feature[4] = 1.2;
        feature[5] = -0.7;
        feature[6] = 2.3;
        feature
    }

    fn hash_of(variant: &[f32], hash_type: HashType, nbin_dist: usize, nbin_angle: usize) -> u32 {
        GeometricHash::perfect_hash_as_u32(&variant.to_vec(), hash_type, nbin_dist, nbin_angle)
    }

    #[test]
    fn radius_one_moves_one_dimension_at_a_time() {
        let tolerance = ToleranceConfig::new(vec![0.5], vec![5.0], 1);
        let mut expander = FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance);
        // 2 distance + 3 angle dimensions, two signs each
        assert_eq!(expander.dim_count(), 5);
        let mut count = 0;
        expander.for_each_neighbor(&pdbtr_feature(), |_| {
            count += 1;
            true
        });
        assert_eq!(count, 10);
    }

    #[test]
    fn radius_two_adds_joint_neighbors() {
        let feature = pdbtr_feature();
        let mut hashes = Vec::new();
        for radius in [1usize, 2] {
            let tolerance = ToleranceConfig::new(vec![0.5], vec![5.0], radius);
            let mut found = HashSet::default();
            FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance)
                .for_each_neighbor(&feature, |variant| {
                    found.insert(hash_of(variant, HashType::PDBTrRosetta, 16, 4));
                    true
                });
            hashes.push(found);
        }
        // Radius 2 is a strict superset: it still walks every single-dimension
        // neighbour, then adds the joint ones.
        assert!(hashes[0].is_subset(&hashes[1]));
        assert!(hashes[1].len() > hashes[0].len());
    }

    #[test]
    fn wide_tolerance_is_substepped_and_leaves_no_gap() {
        // 1.2 A bins: a 3.0 A tolerance must reach the bins in between
        let tolerance = ToleranceConfig::new(vec![3.0], vec![], 1);
        let mut expander = FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance);
        let feature = pdbtr_feature();
        let mut values: Vec<f32> = vec![feature[2]];
        expander.for_each_neighbor(&feature, |variant| {
            values.push(variant[2]);
            true
        });
        values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        values.dedup_by(|a, b| (*a - *b).abs() < 1e-4);
        assert!(values.len() >= 7, "expected sub-steps, got {:?}", values);
        assert!((values[0] - (feature[2] - 3.0)).abs() < 1e-3);
        assert!((values[values.len() - 1] - (feature[2] + 3.0)).abs() < 1e-3);
        // No gap wider than one bin anywhere in the covered interval
        for pair in values.windows(2) {
            assert!(pair[1] - pair[0] <= 1.2001, "gap {:?} skips a bin", pair);
        }
    }

    #[test]
    fn smaller_thresholds_reach_no_extra_bin() {
        // `-d 0.5,1.0` and `-d 1.0` reach the same bins
        let feature = pdbtr_feature();
        let mut sets = Vec::new();
        for thresholds in [vec![0.5, 1.0], vec![1.0]] {
            let tolerance = ToleranceConfig::new(thresholds, vec![], 2);
            let mut found = HashSet::default();
            FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance)
                .for_each_neighbor(&feature, |variant| {
                    found.insert(hash_of(variant, HashType::PDBTrRosetta, 16, 4));
                    true
                });
            sets.push(found);
        }
        assert_eq!(sets[0], sets[1]);
    }

    #[test]
    fn torsion_offsets_wrap_around_pi() {
        // 179 degrees must reach -179 degrees (2 degrees away)
        let mut feature = pdbtr_feature();
        feature[5] = 179.0_f32.to_radians();
        let tolerance = ToleranceConfig::new(vec![], vec![5.0], 1);
        let mut expander = FeatureExpander::new(HashType::FolddiscoDist, 32, 16, &tolerance);
        let mut wrapped = false;
        expander.for_each_neighbor(&feature, |variant| {
            // Every visited value stays inside the encodable domain
            assert!(variant[5] >= -PI - 1e-4 && variant[5] <= PI + 1e-4, "{}", variant[5]);
            if variant[5] < -3.0 {
                wrapped = true;
            }
            true
        });
        assert!(wrapped, "torsion offset past PI did not wrap");
    }

    #[test]
    fn perturbed_angles_stay_inside_their_bit_field() {
        // FolddiscoDist packs theta1/theta2 in 4 bits and Ca-Cb angle in 3; a wide
        // tolerance must not overflow into neighbouring fields.
        let tolerance = ToleranceConfig::new(vec![], vec![60.0], 3);
        let mut expander = FeatureExpander::new(HashType::FolddiscoDist, 32, 16, &tolerance);
        for &(ca_cb, theta1, theta2) in &[
            (0.0f32, PI, -PI), (PI, -PI, PI), (0.05, 3.10, -3.10),
        ] {
            let mut feature = pdbtr_feature();
            feature[4] = ca_cb;
            feature[5] = theta1;
            feature[6] = theta2;
            expander.for_each_neighbor(&feature, |variant| {
                assert!((0.0..=PI).contains(&variant[4]), "Ca-Cb angle {}", variant[4]);
                assert!((-PI..=PI).contains(&variant[5]), "theta1 {}", variant[5]);
                assert!((-PI..=PI).contains(&variant[6]), "theta2 {}", variant[6]);
                let hash = hash_of(variant, HashType::FolddiscoDist, 32, 16);
                // Residue pair: bits 21..29; distances: 11..20
                assert_eq!((hash >> 21) & 0x1FF, 8 * 20 + 3,
                    "angle overflowed into the residue field");
                assert_eq!((hash >> 11) & 0x3FF, expected_distance_bits(&feature),
                    "angle overflowed into the distance fields");
                true
            });
        }
    }

    /// Distance bits of an unperturbed FolddiscoDist hash, for overflow checks.
    fn expected_distance_bits(feature: &[f32]) -> u32 {
        (hash_of(feature, HashType::FolddiscoDist, 32, 16) >> 11) & 0x3FF
    }

    /// Hashes a query reaches for one residue pair, observed geometry included.
    fn reachable_hashes(
        feature: &[f32], hash_type: HashType, nbin_dist: usize, nbin_angle: usize,
        tolerance: &ToleranceConfig,
    ) -> HashSet<u32> {
        let mut found = HashSet::default();
        found.insert(hash_of(feature, hash_type, nbin_dist, nbin_angle));
        FeatureExpander::new(hash_type, nbin_dist, nbin_angle, tolerance)
            .for_each_neighbor(feature, |variant| {
                found.insert(hash_of(variant, hash_type, nbin_dist, nbin_angle));
                true
            });
        found
    }

    #[test]
    fn radius_two_finds_a_pair_that_crosses_two_boundaries_at_once() {
        // Boundary at 7.4 A; query and target straddle it in both Ca and Cb distance
        let mut query = pdbtr_feature();
        query[2] = 7.35;
        query[3] = 7.35;
        let mut target = query.clone();
        target[2] = 7.45;
        target[3] = 7.45;
        let target_hash = hash_of(&target, HashType::PDBTrRosetta, 16, 4);
        assert_ne!(target_hash, hash_of(&query, HashType::PDBTrRosetta, 16, 4));

        let radius_one = reachable_hashes(
            &query, HashType::PDBTrRosetta, 16, 4,
            &ToleranceConfig::new(vec![0.5], vec![5.0], 1),
        );
        let radius_two = reachable_hashes(
            &query, HashType::PDBTrRosetta, 16, 4,
            &ToleranceConfig::new(vec![0.5], vec![5.0], 2),
        );
        assert!(!radius_one.contains(&target_hash), "radius 1 should not reach a joint shift");
        assert!(radius_two.contains(&target_hash), "radius 2 missed the joint shift");
    }

    #[test]
    fn wrapping_finds_a_torsion_across_the_pi_boundary() {
        // 179 vs -179 degrees: opposite ends of the encoded range
        let mut query = pdbtr_feature();
        query[5] = 179.0_f32.to_radians();
        let mut target = query.clone();
        target[5] = -179.0_f32.to_radians();
        let target_hash = hash_of(&target, HashType::FolddiscoDist, 32, 16);
        assert_ne!(target_hash, hash_of(&query, HashType::FolddiscoDist, 32, 16));

        let reachable = reachable_hashes(
            &query, HashType::FolddiscoDist, 32, 16,
            &ToleranceConfig::new(vec![0.5], vec![5.0], 1),
        );
        assert!(reachable.contains(&target_hash), "wrapped torsion neighbour not reached");
    }

    #[test]
    fn visit_can_stop_early() {
        let tolerance = ToleranceConfig::new(vec![0.5], vec![5.0], 3);
        let mut expander = FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance);
        let mut count = 0;
        expander.for_each_neighbor(&pdbtr_feature(), |_| {
            count += 1;
            count < 4
        });
        assert_eq!(count, 4);
    }

    #[test]
    fn wrap_and_reflect_helpers() {
        assert!((wrap_to_pi(PI + 0.2) - (-PI + 0.2)).abs() < 1e-5);
        assert!((wrap_to_pi(-PI - 0.2) - (PI - 0.2)).abs() < 1e-5);
        assert!((wrap_to_pi(0.3) - 0.3).abs() < 1e-6);
        assert!((reflect_into_range(-0.2, 0.0, PI) - 0.2).abs() < 1e-5);
        assert!((reflect_into_range(PI + 0.2, 0.0, PI) - (PI - 0.2)).abs() < 1e-5);
        assert!((reflect_into_range(1.0, 0.0, PI) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn matching_tolerance_adds_only_geometric_index_expansion() {
        let query = ToleranceConfig::new(vec![0.5], vec![5.0], 1);
        let geometric = IndexExpansion::new(1, 0.5, 5.0, None).unwrap();
        assert_eq!(geometric.matching_tolerance(&query), ToleranceConfig::new(vec![1.0], vec![10.0], 2));
        assert_eq!(geometric.matching_tolerance(&ToleranceConfig::none()), ToleranceConfig::new(vec![0.5], vec![5.0], 1));
        let residues_only = IndexExpansion::new(0, 0.5, 5.0, Some(SubstitutionScheme::Group)).unwrap();
        assert_eq!(residues_only.matching_tolerance(&query), query);
        assert!(IndexExpansion::new(0, 0.5, 5.0, None).is_none());
    }

    #[test]
    fn empty_tolerance_produces_nothing() {
        let mut expander = FeatureExpander::new(
            HashType::PDBTrRosetta, 16, 4, &ToleranceConfig::none()
        );
        assert!(expander.is_noop());
        let mut count = 0;
        expander.for_each_neighbor(&pdbtr_feature(), |_| { count += 1; true });
        assert_eq!(count, 0);
    }

    #[test]
    fn every_hash_type_keeps_wide_angle_tolerance_inside_its_encoding() {
        // A wide angle tolerance keeps every hash type inside its declared bit width.
        // Distances stay mid-range; out-of-window distances are not bounded here.
        let tolerance = ToleranceConfig::new(vec![1.0], vec![60.0], 2);
        for hash_type in [
            HashType::PDBMotif, HashType::PDBMotifSinCos, HashType::TrRosetta,
            HashType::PDBTrRosetta, HashType::PointPairFeature,
            HashType::TertiaryInteraction, HashType::Hybrid,
            HashType::FolddiscoAngle, HashType::FolddiscoDist,
        ] {
            let nbin_dist = hash_type.default_dist_bin();
            let nbin_angle = hash_type.default_angle_bin();
            let mut expander = FeatureExpander::new(hash_type, nbin_dist, nbin_angle, &tolerance);
            let mut feature = pdbtr_feature();
            if hash_type == HashType::PDBMotif {
                feature[4] = 120.0; // this one keeps angles in degrees
            }
            if hash_type == HashType::TertiaryInteraction {
                feature[7] = 12.0;  // Ca distance
                feature[8] = 3.0;   // sequence separation
            }
            let bits = hash_type.encoding_bits() as u32;
            expander.for_each_neighbor(&feature, |variant| {
                let hash = hash_of(variant, hash_type, nbin_dist, nbin_angle);
                if bits < 32 {
                    assert_eq!(hash >> bits, 0, "{:?} exceeded {} bits", hash_type, bits);
                }
                true
            });
        }
    }

    /// Prints every expanded feature vector for one query, and times how long
    /// generating them takes. Run with `cargo test print_expanded_feature_sets --
    /// --nocapture` to see the output.
    #[test]
    fn print_expanded_feature_sets() {
        let tolerance = ToleranceConfig::new(vec![0.5], vec![5.0], 2);
        let mut expander = FeatureExpander::new(HashType::PDBTrRosetta, 16, 4, &tolerance);
        let feature = pdbtr_feature();

        println!("query feature: {:?}", feature);

        let mut variants = Vec::new();
        let start = std::time::Instant::now();
        expander.for_each_neighbor(&feature, |variant| {
            variants.push(variant.to_vec());
            true
        });
        let elapsed = start.elapsed();

        for (i, variant) in variants.iter().enumerate() {
            println!("variant {:4}: {:?}", i, variant);
        }
        println!(
            "generated {} additional feature sets in {:?} ({:?}/variant)",
            variants.len(),
            elapsed,
            variants.len().checked_sub(0).filter(|&n| n > 0)
                .map(|n| elapsed / n as u32)
                .unwrap_or_default(),
        );

        assert!(!variants.is_empty());
    }
}
