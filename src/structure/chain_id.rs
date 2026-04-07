/// A fixed-size chain identifier stored as up to 4 ASCII bytes, null-padded on the right.
///
/// Single-character chain IDs (the common case) are stored as `[b'A', 0, 0, 0]`.
/// Multi-character chain IDs like `"AA"` are stored as `[b'A', b'A', 0, 0]`.
pub type ChainId = [u8; 4];

/// Convert a string slice into a `ChainId`.
///
/// Up to 4 bytes are copied; any remaining positions are filled with `0`.
/// Characters beyond the 4th are silently truncated.
pub fn chain_id_from_str(s: &str) -> ChainId {
    let bytes = s.as_bytes();
    let mut id = [0u8; 4];
    for (i, &b) in bytes.iter().take(4).enumerate() {
        id[i] = b;
    }
    id
}

/// Wrap a single ASCII byte into a `ChainId`.
pub fn chain_id_from_byte(b: u8) -> ChainId {
    [b, 0, 0, 0]
}

/// Return the string representation of a `ChainId`, trimming trailing null bytes.
///
/// # Panics
/// Panics if the bytes are not valid UTF-8 (which should never happen when
/// chain IDs are constructed via `chain_id_from_str` or `chain_id_from_byte`).
pub fn chain_id_to_str(id: &ChainId) -> &str {
    let end = id.iter().position(|&b| b == 0).unwrap_or(4);
    std::str::from_utf8(&id[..end]).expect("ChainId contains non-UTF-8 bytes")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chain_id_from_str_single() {
        assert_eq!(chain_id_from_str("A"), [b'A', 0, 0, 0]);
    }

    #[test]
    fn test_chain_id_from_str_double() {
        assert_eq!(chain_id_from_str("AA"), [b'A', b'A', 0, 0]);
    }

    #[test]
    fn test_chain_id_from_str_four() {
        assert_eq!(chain_id_from_str("ABCD"), [b'A', b'B', b'C', b'D']);
    }

    #[test]
    fn test_chain_id_from_byte() {
        assert_eq!(chain_id_from_byte(b'B'), [b'B', 0, 0, 0]);
    }

    #[test]
    fn test_chain_id_to_str_single() {
        assert_eq!(chain_id_to_str(&[b'A', 0, 0, 0]), "A");
    }

    #[test]
    fn test_chain_id_to_str_double() {
        assert_eq!(chain_id_to_str(&[b'A', b'A', 0, 0]), "AA");
    }

    #[test]
    fn test_chain_id_roundtrip() {
        for s in &["A", "AA", "AAA", "AAAA", "10", "BB"] {
            assert_eq!(chain_id_to_str(&chain_id_from_str(s)), *s);
        }
    }
}
