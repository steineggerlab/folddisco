//! Chain identifiers wide enough for multi-character mmCIF asym IDs (`AA`, `10`).
//!
//! [`ChainId`] is inline and `Copy` to stay allocation-free in the retrieval hot path.
//! [`crate::structure::atom::Atom`] keeps a single-byte `chain` because its layout
//! must match Foldcomp's C `atom_t`.

use std::fmt;

/// Longest chain ID stored in full. PDB asym IDs are at most 4 characters; 8 costs
/// nothing extra since `Option<(ChainId, u64)>` is 24 bytes either way.
pub const CHAIN_ID_MAX_LEN: usize = 8;

/// Substitute for non-printable bytes, so [`ChainId::as_str`] never fails.
const REPLACEMENT: u8 = b'?';

/// A chain identifier of up to [`CHAIN_ID_MAX_LEN`] printable ASCII bytes,
/// null-padded on the right.
#[derive(Copy, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub struct ChainId {
    bytes: [u8; CHAIN_ID_MAX_LEN],
}

impl ChainId {
    /// Empty chain ID; never parsed from a file, so usable as a sentinel.
    pub const fn empty() -> Self {
        ChainId { bytes: [0; CHAIN_ID_MAX_LEN] }
    }

    /// Wrap a single byte, as read from column 22 of a PDB `ATOM` record.
    pub fn from_byte(byte: u8) -> Self {
        let mut bytes = [0u8; CHAIN_ID_MAX_LEN];
        bytes[0] = sanitize(byte);
        ChainId { bytes }
    }

    /// Take up to [`CHAIN_ID_MAX_LEN`] bytes; longer input is truncated.
    pub fn from_bytes(bytes: &[u8]) -> Self {
        let mut out = [0u8; CHAIN_ID_MAX_LEN];
        for (slot, &byte) in out.iter_mut().zip(bytes.iter()) {
            *slot = sanitize(byte);
        }
        ChainId { bytes: out }
    }

    /// Take up to [`CHAIN_ID_MAX_LEN`] bytes of a string; longer input is truncated.
    pub fn from_str(text: &str) -> Self {
        Self::from_bytes(text.as_bytes())
    }

    /// Number of bytes actually used, 0 for [`ChainId::empty`].
    pub fn len(&self) -> usize {
        self.bytes.iter().position(|&b| b == 0).unwrap_or(CHAIN_ID_MAX_LEN)
    }

    pub fn is_empty(&self) -> bool {
        self.bytes[0] == 0
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes[..self.len()]
    }

    pub fn as_str(&self) -> &str {
        // Every byte is forced into printable ASCII by `sanitize`, so this holds.
        std::str::from_utf8(self.as_bytes()).unwrap_or("?")
    }

    /// First byte, for the single-byte `Atom::chain`. Lossy for multi-character IDs.
    pub fn first_byte(&self) -> u8 {
        self.bytes[0]
    }

    /// Whether `chain` + `residue` is ambiguous without `_`: `A21` parses back,
    /// `1021` (chain `10`) and `AA250` do not.
    pub fn needs_separator(&self) -> bool {
        self.len() != 1 || !self.bytes[0].is_ascii_alphabetic()
    }
}

/// Map non-printable bytes onto [`REPLACEMENT`] to keep `ChainId` valid UTF-8.
fn sanitize(byte: u8) -> u8 {
    if byte.is_ascii_graphic() || byte == b' ' { byte } else { REPLACEMENT }
}

impl fmt::Display for ChainId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl fmt::Debug for ChainId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ChainId({:?})", self.as_str())
    }
}

impl From<u8> for ChainId {
    fn from(byte: u8) -> Self {
        ChainId::from_byte(byte)
    }
}

impl From<&str> for ChainId {
    fn from(text: &str) -> Self {
        ChainId::from_str(text)
    }
}

pub fn chain(text: &str) -> ChainId {
    ChainId::from_str(text)
}

/// Whether a residue list needs `_` separators. Decided per list, so one field
/// never mixes both spellings.
pub fn residue_list_needs_separator<'a, I>(chains: I) -> bool
where
    I: IntoIterator<Item = &'a ChainId>,
{
    chains.into_iter().any(|chain| chain.needs_separator())
}

/// `A21` or `A_21`, depending on `separator`.
pub fn format_chain_residue(chain: &ChainId, residue: u64, separator: bool) -> String {
    if separator {
        format!("{}_{}", chain.as_str(), residue)
    } else {
        format!("{}{}", chain.as_str(), residue)
    }
}

/// Split one query/output token into its chain and residue parts.
///
/// Accepts, in order of precedence:
///   `AA_250` / `A_250` / `10_250` — explicit separator, any chain ID
///   `A250`                        — legacy single leading letter
///   `250`                         — no chain, caller's default applies
///
/// Returns `(None, rest)` when the token carries no chain of its own.
pub fn split_chain_and_rest(token: &str) -> (Option<ChainId>, &str) {
    if let Some(sep) = token.find('_') {
        let (chain, rest) = (&token[..sep], &token[sep + 1..]);
        // A leading `_` is the "residue not matched" placeholder in output, and
        // carries no chain; treat it like a bare residue index.
        if chain.is_empty() {
            return (None, rest);
        }
        return (Some(ChainId::from_str(chain)), rest);
    }
    match token.as_bytes().first() {
        Some(first) if first.is_ascii_alphabetic() => {
            (Some(ChainId::from_byte(*first)), &token[1..])
        }
        _ => (None, token),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_str_roundtrips() {
        for text in &["A", "AA", "10", "1", "ABCD", "AAAAAAAA", "a"] {
            assert_eq!(ChainId::from_str(text).as_str(), *text);
            assert_eq!(ChainId::from_str(text).len(), text.len());
        }
    }

    #[test]
    fn from_byte_matches_from_str() {
        assert_eq!(ChainId::from_byte(b'A'), ChainId::from_str("A"));
        assert_eq!(ChainId::from_byte(b'A').first_byte(), b'A');
    }

    #[test]
    fn empty_is_distinct_from_every_real_chain() {
        assert!(ChainId::empty().is_empty());
        assert_eq!(ChainId::empty().len(), 0);
        assert_eq!(ChainId::empty().as_str(), "");
        assert_ne!(ChainId::empty(), ChainId::from_str("A"));
        assert_ne!(ChainId::empty(), ChainId::from_byte(b' '));
    }

    #[test]
    fn over_long_ids_truncate_rather_than_panic() {
        let long = ChainId::from_str("ABCDEFGHIJ");
        assert_eq!(long.len(), CHAIN_ID_MAX_LEN);
        assert_eq!(long.as_str(), "ABCDEFGH");
    }

    #[test]
    fn non_ascii_bytes_stay_printable() {
        assert_eq!(ChainId::from_bytes(&[0xff]).as_str(), "?");
        assert_eq!(ChainId::from_bytes(&[b'A', 0x01]).as_str(), "A?");
    }

    #[test]
    fn separator_is_needed_exactly_when_concatenation_is_ambiguous() {
        assert!(!ChainId::from_str("A").needs_separator());
        assert!(!ChainId::from_str("z").needs_separator());
        assert!(ChainId::from_str("AA").needs_separator());
        assert!(ChainId::from_str("1").needs_separator());
        assert!(ChainId::from_str("10").needs_separator());
        assert!(ChainId::empty().needs_separator());
    }

    #[test]
    fn residue_list_decision_is_all_or_nothing() {
        let legacy = vec![ChainId::from_str("A"), ChainId::from_str("B")];
        assert!(!residue_list_needs_separator(legacy.iter()));
        let mixed = vec![ChainId::from_str("A"), ChainId::from_str("10")];
        assert!(residue_list_needs_separator(mixed.iter()));
    }

    #[test]
    fn formatting_follows_the_separator_flag() {
        let chain = ChainId::from_str("A");
        assert_eq!(format_chain_residue(&chain, 21, false), "A21");
        assert_eq!(format_chain_residue(&chain, 21, true), "A_21");
        assert_eq!(format_chain_residue(&ChainId::from_str("10"), 21, true), "10_21");
    }

    #[test]
    fn split_accepts_both_grammars() {
        assert_eq!(split_chain_and_rest("A250"), (Some(ChainId::from_str("A")), "250"));
        assert_eq!(split_chain_and_rest("A_250"), (Some(ChainId::from_str("A")), "250"));
        assert_eq!(split_chain_and_rest("AA_250"), (Some(ChainId::from_str("AA")), "250"));
        assert_eq!(split_chain_and_rest("10_250"), (Some(ChainId::from_str("10")), "250"));
        assert_eq!(split_chain_and_rest("250"), (None, "250"));
        assert_eq!(split_chain_and_rest("_250"), (None, "250"));
        assert_eq!(split_chain_and_rest("A250-252"), (Some(ChainId::from_str("A")), "250-252"));
        assert_eq!(split_chain_and_rest("AA_250-252"), (Some(ChainId::from_str("AA")), "250-252"));
    }
}
