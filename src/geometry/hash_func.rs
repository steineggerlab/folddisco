
// Invertible hash function for 64-bit unsigned integers
// Originally from https://aebou.rbind.io/post/a-rust-glimpse-at-thomas-wang-integer-hash-function/
#[inline]
pub fn hash_64(key: u64) -> u64 {
    let mut h_key = key;
    // key = (key << 21) - key - 1
    h_key = (!h_key).wrapping_add(h_key << 21);
    h_key = h_key ^ h_key >> 24;
    // key * 265
    h_key = h_key.wrapping_add(h_key << 3).wrapping_add(h_key << 8);
    h_key = h_key ^ h_key >> 14;
    // key * 21
    h_key = h_key.wrapping_add(h_key << 2).wrapping_add(h_key << 4);
    h_key = h_key ^ h_key >> 28;
    h_key = h_key.wrapping_add(h_key << 31);
    h_key
}

#[inline]
pub fn hash_64_inv(hashed_key: u64) -> u64 {
    let mut key = hashed_key;
    // Invert h_key = h_key.wrapping_add(h_key << 31)
    let mut tmp: u64 = key.wrapping_sub(key << 31);
    key = key.wrapping_sub(tmp << 31);
    // Invert h_key = h_key ^ h_key >> 28;
    tmp = key ^ key >> 28;
    key = key ^ tmp >> 28;
    // Invert h_key = h_key.wrapping_add(h_key << 2).wrapping_add(h_key << 4)
    key = key.wrapping_mul(14933078535860113213u64);
    // Invert h_key = h_key ^ h_key >> 14;
    tmp = key ^ key >> 14;
    tmp = key ^ tmp >> 14;
    tmp = key ^ tmp >> 14;
    key = key ^ tmp >> 14;
    // Invert h_key = h_key.wrapping_add(h_key << 3).wrapping_add(h_key << 8)
    key = key.wrapping_mul(15244667743933553977u64);
    // Invert h_key = h_key ^ h_key >> 24
    tmp = key ^ key >> 24;
    key = key ^ tmp >> 24;
    // Invert h_key = (!h_key).wrapping_add(h_key << 21)
    tmp = !key;
    tmp = !(key.wrapping_sub(tmp << 21));
    tmp = !(key.wrapping_sub(tmp << 21));
    key = !(key.wrapping_sub(tmp << 21));
    key
}

#[cfg(test)]
mod hash_func_test {

    #[test]
    fn test_hash_64_comparison_inline() {
        let start = std::time::Instant::now();
        for i in 0..10000000 {
            let key = i;
            let h_key = super::hash_64(key);
            let inv_key = super::hash_64_inv(h_key);
            assert_eq!(key, inv_key);
        }
        let end = std::time::Instant::now();
        println!("hash_64: {:?}", end - start);
    }

}

// pub struct U64Hasher;

// impl Hasher for U64Hasher {
//     fn finish(&self) -> u64 {
//         0
//     }

//     fn write(&mut self, _bytes: &[u8]) {
//         unimplemented!()
//     }

//     fn write_u64(&mut self, i: u64) {
//         unimplemented!()
//     }
// }