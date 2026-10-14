//! Amino acid substitution schemes for query and index expansion.
//!
//! Residue codes follow `map_aa_to_u8`: A R N D C Q E G H I L K M F P S T W Y V (0-19),
//! which is also the BLOSUM62 row order.

use std::fmt;

/// Marker for `:*` in a query residue; resolved to the active scheme's alternatives.
pub const SCHEME_MARKER: u8 = 254;

/// How alternatives for an observed amino acid are chosen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubstitutionScheme {
    /// Positive BLOSUM62 score.
    Blosum62,
    /// Same physicochemical class: RHK, DE, NQST, FWY, AVLIMC, GP.
    Group,
    /// Same IMGT side-chain volume class: GAS, CDPNT, QEHV, MILKR, FWY.
    Size,
}

impl SubstitutionScheme {
    /// Parse a CLI value; `None` for an unknown name.
    pub fn from_str(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "blosum62" | "blosum" => Some(Self::Blosum62),
            "group" => Some(Self::Group),
            "size" => Some(Self::Size),
            _ => None,
        }
    }

    /// Alternatives for `aa`, excluding `aa` itself. Empty for non-standard codes.
    pub fn alternatives(&self, aa: u8) -> Vec<u8> {
        let aa = aa as usize;
        if aa >= 20 {
            return Vec::new();
        }
        (0..20u8).filter(|&other| other as usize != aa && match self {
            Self::Blosum62 => BLOSUM62[aa][other as usize] > 0,
            Self::Group => GROUP[aa] == GROUP[other as usize],
            Self::Size => SIZE_CLASS[aa] == SIZE_CLASS[other as usize],
        }).collect()
    }
}

impl fmt::Display for SubstitutionScheme {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(match self {
            Self::Blosum62 => "blosum62",
            Self::Group => "group",
            Self::Size => "size",
        })
    }
}

#[rustfmt::skip]
const BLOSUM62: [[i8; 20]; 20] = [
    //A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0], // A
    [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3], // R
    [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3], // N
    [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3], // D
    [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1], // C
    [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2], // Q
    [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2], // E
    [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3], // G
    [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3], // H
    [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3], // I
    [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1], // L
    [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2], // K
    [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1], // M
    [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1], // F
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2], // P
    [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2], // S
    [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0], // T
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3], // W
    [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1], // Y
    [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4], // V
];

// 0 positive, 1 negative, 2 polar, 3 aromatic, 4 aliphatic, 5 special
#[rustfmt::skip]
const GROUP: [u8; 20] = [
//  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
    4, 0, 2, 1, 4, 2, 1, 5, 0, 4, 4, 0, 4, 3, 5, 2, 2, 3, 3, 4,
];

// IMGT volume classes: 0 very small, 1 small, 2 medium, 3 large, 4 very large
#[rustfmt::skip]
const SIZE_CLASS: [u8; 20] = [
//  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
    0, 3, 1, 1, 1, 2, 2, 0, 2, 3, 3, 3, 3, 4, 1, 0, 1, 4, 4, 2,
];

/// Amino acid pairs to query for one residue pair besides the observed one.
///
/// Each side's alternatives are crossed with the other side and with the observed residue,
/// so `A164:H,A200:ND` reaches (H,D), (H,N), (H,obs), (obs,D) and (obs,N).
/// Codes outside 0-19 (unknown letters, an unresolved `SCHEME_MARKER`) are skipped.
pub fn substitution_variants(
    observed: (f32, f32), sub_i: Option<&[u8]>, sub_j: Option<&[u8]>,
) -> Vec<(f32, f32)> {
    if sub_i.is_none() && sub_j.is_none() {
        return Vec::new();
    }
    let alternatives = |observed_aa: f32, subs: Option<&[u8]>| -> Vec<f32> {
        let mut out = vec![observed_aa];
        for &aa in subs.unwrap_or(&[]) {
            if aa >= 20 {
                continue;
            }
            let aa = aa as f32;
            if !out.contains(&aa) {
                out.push(aa);
            }
        }
        out
    };
    let alt_i = alternatives(observed.0, sub_i);
    let alt_j = alternatives(observed.1, sub_j);

    let mut variants = Vec::with_capacity(alt_i.len() * alt_j.len());
    for &aa_i in &alt_i {
        for &aa_j in &alt_j {
            if aa_i != observed.0 || aa_j != observed.1 {
                variants.push((aa_i, aa_j));
            }
        }
    }
    variants
}

/// Replace `SCHEME_MARKER` in a residue's substitution list with `scheme`'s
/// alternatives for the observed amino acid. With `apply_to_all`, residues without an
/// explicit list get the scheme too.
pub fn resolve_substitution(
    substitution: &Option<Vec<u8>>, observed_aa: u8, scheme: SubstitutionScheme, apply_to_all: bool,
) -> Option<Vec<u8>> {
    match substitution {
        Some(subs) if subs.contains(&SCHEME_MARKER) => {
            let mut resolved: Vec<u8> = subs.iter().copied().filter(|&aa| aa != SCHEME_MARKER).collect();
            for aa in scheme.alternatives(observed_aa) {
                if !resolved.contains(&aa) {
                    resolved.push(aa);
                }
            }
            Some(resolved)
        }
        Some(subs) => Some(subs.clone()),
        None if apply_to_all => Some(scheme.alternatives(observed_aa)),
        None => None,
    }
}

/// Everything an index expanded with `scheme` can pair with this query residue: the
/// scheme's alternatives of the observed residue and of every `substitution` code.
pub fn index_matching_substitution(
    substitution: &Option<Vec<u8>>, observed_aa: u8, scheme: SubstitutionScheme,
) -> Option<Vec<u8>> {
    let mut out: Vec<u8> = substitution.clone().unwrap_or_default();
    let sources: Vec<u8> = std::iter::once(observed_aa).chain(out.iter().copied()).collect();
    for aa in sources.into_iter().flat_map(|aa| scheme.alternatives(aa)) {
        if !out.contains(&aa) {
            out.push(aa);
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const ORDER: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

    fn code(letter: u8) -> u8 {
        ORDER.iter().position(|&c| c == letter).unwrap() as u8
    }

    fn letters(codes: &[u8]) -> String {
        let mut out: Vec<char> = codes.iter().map(|&c| ORDER[c as usize] as char).collect();
        out.sort_unstable();
        out.into_iter().collect()
    }

    #[test]
    fn order_matches_map_aa_to_u8() {
        use crate::utils::convert::map_one_letter_to_u8_vec;
        for (i, &letter) in ORDER.iter().enumerate() {
            assert_eq!(map_one_letter_to_u8_vec(letter as char), vec![i as u8]);
        }
    }

    #[test]
    fn blosum62_is_symmetric_with_positive_diagonal() {
        for i in 0..20 {
            assert!(BLOSUM62[i][i] > 0);
            for j in 0..20 {
                assert_eq!(BLOSUM62[i][j], BLOSUM62[j][i], "{} {}", ORDER[i] as char, ORDER[j] as char);
            }
        }
    }

    #[test]
    fn blosum62_positive_pairs_match_the_published_matrix() {
        let mut pairs = Vec::new();
        for i in 0..20 {
            for j in (i + 1)..20 {
                if BLOSUM62[i][j] > 0 {
                    pairs.push(format!("{}{}", ORDER[i] as char, ORDER[j] as char));
                }
            }
        }
        assert_eq!(pairs, [
            "AS", "RQ", "RK", "ND", "NH", "NS", "DE", "QE", "QK", "EK", "HY",
            "IL", "IM", "IV", "LM", "LV", "MV", "FW", "FY", "ST", "WY",
        ]);
        assert_eq!(BLOSUM62[code(b'W') as usize][code(b'W') as usize], 11);
        assert_eq!(BLOSUM62[code(b'L') as usize][code(b'D') as usize], -4);
    }

    #[test]
    fn scheme_alternatives() {
        let blosum = SubstitutionScheme::Blosum62;
        assert_eq!(letters(&blosum.alternatives(code(b'I'))), "LMV");
        assert_eq!(letters(&blosum.alternatives(code(b'D'))), "EN");
        assert_eq!(letters(&blosum.alternatives(code(b'H'))), "NY");
        assert!(blosum.alternatives(code(b'C')).is_empty());
        assert!(blosum.alternatives(code(b'G')).is_empty());
        assert_eq!(letters(&SubstitutionScheme::Group.alternatives(code(b'K'))), "HR");
        assert_eq!(letters(&SubstitutionScheme::Group.alternatives(code(b'D'))), "E");
        assert_eq!(letters(&SubstitutionScheme::Size.alternatives(code(b'W'))), "FY");
        assert_eq!(letters(&SubstitutionScheme::Size.alternatives(code(b'G'))), "AS");
        assert!(blosum.alternatives(255).is_empty());
    }

    #[test]
    fn groups_and_size_classes_partition_all_residues() {
        for table in [&GROUP, &SIZE_CLASS] {
            for i in 0..20 {
                let class: Vec<u8> = (0..20u8).filter(|&j| table[j as usize] == table[i]).collect();
                assert!(class.contains(&(i as u8)));
            }
        }
        assert_eq!((0..20).filter(|&i| GROUP[i] == 5).count(), 2); // G, P
        assert_eq!((0..20).filter(|&i| SIZE_CLASS[i] == 4).count(), 3); // F, W, Y
    }

    #[test]
    fn parse_scheme_names() {
        assert_eq!(SubstitutionScheme::from_str("BLOSUM62"), Some(SubstitutionScheme::Blosum62));
        assert_eq!(SubstitutionScheme::from_str("size"), Some(SubstitutionScheme::Size));
        assert_eq!(SubstitutionScheme::from_str("pam250"), None);
    }

    #[test]
    fn resolve_marker_and_apply_to_all() {
        let scheme = SubstitutionScheme::Blosum62;
        let d = code(b'D');
        // `:*` expands to the scheme, keeping explicit letters given with it
        assert_eq!(
            resolve_substitution(&Some(vec![SCHEME_MARKER]), d, scheme, false).map(|v| letters(&v)),
            Some("EN".to_string())
        );
        assert_eq!(
            resolve_substitution(&Some(vec![code(b'Q'), SCHEME_MARKER]), d, scheme, false).map(|v| letters(&v)),
            Some("ENQ".to_string())
        );
        // Explicit lists win over the global scheme
        assert_eq!(resolve_substitution(&Some(vec![code(b'H')]), d, scheme, true), Some(vec![code(b'H')]));
        assert_eq!(resolve_substitution(&None, d, scheme, false), None);
        assert_eq!(resolve_substitution(&None, d, scheme, true).map(|v| letters(&v)), Some("EN".to_string()));
    }

    #[test]
    fn index_matching_includes_alternatives_of_explicit_codes() {
        // `D:K` on a BLOSUM62 index reaches Arg through Lys
        let scheme = SubstitutionScheme::Blosum62;
        let resolved = index_matching_substitution(&Some(vec![code(b'K')]), code(b'D'), scheme).unwrap();
        assert_eq!(letters(&resolved), "EKNQR");
        assert_eq!(letters(&index_matching_substitution(&None, code(b'D'), scheme).unwrap()), "EN");
    }

    #[test]
    fn variants_cross_both_sides() {
        let (h, d, n) = (code(b'H') as f32, code(b'D') as f32, code(b'N') as f32);
        let observed = (0.0, 1.0);
        let variants = substitution_variants(observed, Some(&[code(b'H')]), Some(&[code(b'D'), code(b'N')]));
        assert_eq!(variants.len(), 5);
        for pair in [(h, d), (h, n), (h, 1.0), (0.0, d), (0.0, n)] {
            assert!(variants.contains(&pair));
        }
        assert!(substitution_variants(observed, None, None).is_empty());
        assert!(substitution_variants(observed, Some(&[255, SCHEME_MARKER]), None).is_empty());
    }
}
