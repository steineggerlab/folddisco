//! Memory-mapped caches whose on-disk layout is the in-memory layout: loading is
//! an `mmap` plus a header check. Records are cast directly to `&[R]`; names live
//! in one trailing blob and are resolved on demand.
//!
//! ```text
//! magic [u8; 8] | version u32 | pad u32 | count u64 |
//! source length u64 | source mtime nanos u64 | names length u64 |   (48 bytes)
//! count * size_of::<R>() bytes of records                          (8-aligned)
//! names length bytes of UTF-8                                      (name blob)
//! ```
//!
//! A load checks magic, version, source length + mtime, and that the file size
//! matches the header counts; any failure means "rebuild from source".
//! No checksum (it would cost the full pass the cache avoids) and no
//! temp-file-and-rename: the bytes are deterministic, so concurrent writers emit
//! identical files and a partial write fails the size check. Keep the body
//! deterministic (no timestamps/thread ids) or both assumptions break.
//! Names are UTF-8-checked only when resolved; invalid names read as empty.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem::{align_of, size_of};

use memmap2::Mmap;

/// Header size; a multiple of 8 so the record table is aligned.
pub const HEADER_SIZE: usize = 48;

/// Marks a type as safe to reinterpret from arbitrary cache bytes.
///
/// # Safety
/// Implementors must be `#[repr(C)]`, padding-free, aligned to at most 8, sized
/// in multiples of 8, and valid for every bit pattern (integer/float fields only).
pub unsafe trait CacheRecord: Copy {
    /// Distinguishes cache kinds, so one kind's file is never read as another's.
    const MAGIC: &'static [u8; 8];
    /// Layout version; bump when the record layout changes.
    const VERSION: u32;
}

/// A cache's bytes: a file mapping, or an in-memory image when writing failed.
///
/// The image is stored as `u64` words so the record cast is 8-aligned; `len` is
/// the logical byte length.
enum Backing {
    Mapped(Mmap),
    Owned { words: Vec<u64>, len: usize },
}

impl Backing {
    #[inline]
    fn bytes(&self) -> &[u8] {
        match self {
            Backing::Mapped(mmap) => mmap,
            // SAFETY: `words` holds at least `len` bytes, and any bit pattern is
            // a valid `u8`.
            Backing::Owned { words, len } => unsafe {
                std::slice::from_raw_parts(words.as_ptr() as *const u8, *len)
            },
        }
    }

    /// Copy a serialised image into an 8-aligned word buffer.
    fn own(image: Vec<u8>) -> Self {
        let len = image.len();
        let mut words = vec![0u64; len.div_ceil(8)];
        // SAFETY: the destination holds `words.len() * 8 >= len` bytes.
        unsafe {
            std::ptr::copy_nonoverlapping(
                image.as_ptr(), words.as_mut_ptr() as *mut u8, len,
            );
        }
        Backing::Owned { words, len }
    }
}

/// `(length, mtime nanos)` of `path`, as stored in a cache header.
/// `None` if the file is missing or its mtime is unrepresentable.
pub fn source_identity(path: &str) -> Option<(u64, u64)> {
    let meta = std::fs::metadata(path).ok()?;
    let mtime = meta.modified().ok()?
        .duration_since(std::time::UNIX_EPOCH).ok()?
        .as_nanos();
    Some((meta.len(), u64::try_from(mtime).ok()?))
}

/// A loaded cache: a record table and a name blob over one backing buffer.
pub struct PodCache<R: CacheRecord> {
    backing: Backing,
    count: usize,
    names_offset: usize,
    names_len: usize,
    _record: std::marker::PhantomData<R>,
}

impl<R: CacheRecord> PodCache<R> {
    /// Map `cache_path` if it was built from the current `source_path`.
    /// `None` on any mismatch or corruption, so the caller rebuilds from source.
    pub fn map(cache_path: &str, source_path: &str) -> Option<Self> {
        let file = File::open(cache_path).ok()?;
        let mmap = unsafe { Mmap::map(&file).ok()? };
        Self::from_backing(Backing::Mapped(mmap), source_path)
    }

    fn from_backing(backing: Backing, source_path: &str) -> Option<Self> {
        let bytes = backing.bytes();
        if bytes.len() < HEADER_SIZE || &bytes[..8] != R::MAGIC {
            return None;
        }
        if u32::from_le_bytes(bytes[8..12].try_into().unwrap()) != R::VERSION {
            return None;
        }
        let count = u64::from_le_bytes(bytes[16..24].try_into().unwrap()) as usize;
        let src_len = u64::from_le_bytes(bytes[24..32].try_into().unwrap());
        let src_mtime = u64::from_le_bytes(bytes[32..40].try_into().unwrap());
        let names_len = u64::from_le_bytes(bytes[40..48].try_into().unwrap()) as usize;
        if source_identity(source_path)? != (src_len, src_mtime) {
            return None;
        }

        let names_offset = count.checked_mul(size_of::<R>())?.checked_add(HEADER_SIZE)?;
        if bytes.len() != names_offset.checked_add(names_len)? {
            return None;
        }
        // Guaranteed by construction, but checked because a misaligned cast is UB.
        if bytes.as_ptr() as usize % align_of::<R>() != 0 || HEADER_SIZE % align_of::<R>() != 0 {
            return None;
        }

        Some(PodCache {
            backing, count, names_offset, names_len,
            _record: std::marker::PhantomData,
        })
    }

    /// Load a cache from an in-memory image (used when the directory is read-only).
    fn from_image(image: Vec<u8>, source_path: &str) -> Option<Self> {
        Self::from_backing(Backing::own(image), source_path)
    }

    /// A cache with no entries, for a reader with no database attached.
    pub fn empty() -> Self {
        PodCache {
            backing: Backing::own(vec![0u8; HEADER_SIZE]),
            count: 0,
            names_offset: HEADER_SIZE,
            names_len: 0,
            _record: std::marker::PhantomData,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.count
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// The whole record table, cast out of the mapping.
    #[inline]
    pub fn records(&self) -> &[R] {
        if self.count == 0 {
            // `from_raw_parts` needs an aligned non-null pointer even for length 0.
            return &[];
        }
        let bytes = self.backing.bytes();
        let start = HEADER_SIZE;
        // SAFETY: `from_backing` checked alignment and that `count * size_of::<R>()`
        // bytes follow the header; `CacheRecord` types accept any bit pattern; the
        // slice borrows `self`, which owns the backing.
        unsafe {
            std::slice::from_raw_parts(
                bytes[start..].as_ptr() as *const R,
                self.count,
            )
        }
    }

    /// The raw name blob.
    #[inline]
    pub fn names_blob(&self) -> &[u8] {
        &self.backing.bytes()[self.names_offset..self.names_offset + self.names_len]
    }

    /// Resolve one name; empty if out of bounds or not UTF-8.
    #[inline]
    pub fn name_at(&self, offset: u64, len: u32) -> &str {
        let start = offset as usize;
        let end = start + len as usize;
        let blob = self.names_blob();
        if end > blob.len() {
            return "";
        }
        std::str::from_utf8(&blob[start..end]).unwrap_or("")
    }
}

/// Serialise a cache: header, `records` in reader order, then the `names` blob.
pub fn write_cache<R: CacheRecord, W: Write>(
    sink: &mut W, source_len: u64, source_mtime: u64,
    records: &[R], names: &[u8],
) -> std::io::Result<()> {
    sink.write_all(R::MAGIC)?;
    sink.write_all(&R::VERSION.to_le_bytes())?;
    sink.write_all(&0u32.to_le_bytes())?; // pad, keeps the counts 8-aligned
    sink.write_all(&(records.len() as u64).to_le_bytes())?;
    sink.write_all(&source_len.to_le_bytes())?;
    sink.write_all(&source_mtime.to_le_bytes())?;
    sink.write_all(&(names.len() as u64).to_le_bytes())?;
    // SAFETY: `CacheRecord` types are `#[repr(C)]` and padding-free.
    let record_bytes = unsafe {
        std::slice::from_raw_parts(
            records.as_ptr() as *const u8,
            std::mem::size_of_val(records),
        )
    };
    sink.write_all(record_bytes)?;
    sink.write_all(names)?;
    sink.flush()
}

/// Write a cache and map it back; falls back to an in-memory image if the
/// directory is not writable (read-only index directories are normal).
pub fn store_and_map<R: CacheRecord>(
    cache_path: &str, source_path: &str,
    records: &[R], names: &[u8],
) -> Option<PodCache<R>> {
    let (source_len, source_mtime) = source_identity(source_path)?;
    // 4 MB buffer: caches reach several GB.
    let written = File::create(cache_path).and_then(|file| {
        let mut writer = BufWriter::with_capacity(4 << 20, file);
        write_cache(&mut writer, source_len, source_mtime, records, names)
    });
    match written {
        Ok(()) => {
            if let Some(cache) = PodCache::map(cache_path, source_path) {
                return Some(cache);
            }
            // Written but not mappable: fall back to the in-memory image.
            crate::utils::log::print_log_msg(
                crate::utils::log::WARN,
                &format!("Wrote but could not map the cache {}; using memory", cache_path),
            );
        }
        Err(err) => {
            crate::utils::log::print_log_msg(
                crate::utils::log::WARN,
                &format!("Unable to write the cache {}: {}", cache_path, err),
            );
        }
    }
    let mut image = Vec::new();
    write_cache(&mut image, source_len, source_mtime, records, names).ok()?;
    PodCache::from_image(image, source_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[repr(C)]
    #[derive(Copy, Clone, Debug, PartialEq)]
    struct TestRecord {
        key: u64,
        name_offset: u64,
        name_len: u32,
        _pad: u32,
    }
    unsafe impl CacheRecord for TestRecord {
        const MAGIC: &'static [u8; 8] = b"FDTESTC1";
        const VERSION: u32 = 1;
    }

    /// The alignment and padding promises `CacheRecord` makes.
    #[test]
    fn record_layout_is_castable() {
        assert_eq!(size_of::<TestRecord>() % 8, 0);
        assert!(align_of::<TestRecord>() <= 8);
        assert_eq!(HEADER_SIZE % 8, 0);
    }

    fn temp_paths(tag: &str) -> (String, String) {
        let unique = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
            .unwrap().as_nanos();
        let dir = std::env::temp_dir().to_string_lossy().to_string();
        (format!("{}/fd_pod_{}_{}.src", dir, tag, unique),
         format!("{}/fd_pod_{}_{}.cache", dir, tag, unique))
    }

    fn sample() -> (Vec<TestRecord>, Vec<u8>) {
        let names: Vec<u8> = b"alpha beta gamma".to_vec();
        let records = vec![
            TestRecord { key: 7, name_offset: 0, name_len: 5, _pad: 0 },
            TestRecord { key: 9, name_offset: 6, name_len: 4, _pad: 0 },
            TestRecord { key: 11, name_offset: 11, name_len: 5, _pad: 0 },
        ];
        (records, names)
    }

    #[test]
    fn roundtrip_through_a_file() {
        let (src, cache) = temp_paths("roundtrip");
        std::fs::write(&src, b"source bytes").unwrap();
        let (records, names) = sample();

        let mapped = store_and_map(&cache, &src, &records, &names)
            .expect("cache should store and map");
        assert_eq!(mapped.len(), 3);
        assert_eq!(mapped.records(), records.as_slice());
        assert_eq!(mapped.name_at(0, 5), "alpha");
        assert_eq!(mapped.name_at(6, 4), "beta");
        assert_eq!(mapped.name_at(11, 5), "gamma");

        // And again from a cold map, i.e. the path a second process takes.
        let reopened = PodCache::<TestRecord>::map(&cache, &src).expect("should map");
        assert_eq!(reopened.records(), records.as_slice());

        let _ = std::fs::remove_file(&src);
        let _ = std::fs::remove_file(&cache);
    }

    #[test]
    fn a_source_that_changed_is_rejected() {
        let (src, cache) = temp_paths("stale");
        std::fs::write(&src, b"source bytes").unwrap();
        let (records, names) = sample();
        store_and_map(&cache, &src, &records, &names).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_some());

        // Same length, different content: only the mtime moves.
        std::fs::write(&src, b"SOURCE BYTES").unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        let _ = std::fs::remove_file(&src);
        let _ = std::fs::remove_file(&cache);
    }

    #[test]
    fn truncated_wrong_magic_and_wrong_version_are_rejected() {
        let (src, cache) = temp_paths("corrupt");
        std::fs::write(&src, b"source bytes").unwrap();
        let (records, names) = sample();
        store_and_map(&cache, &src, &records, &names).unwrap();
        let good = std::fs::read(&cache).unwrap();

        // Truncated: the size check fails.
        std::fs::write(&cache, &good[..good.len() - 3]).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        // Shorter than a header.
        std::fs::write(&cache, &good[..12]).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        // Wrong magic, i.e. another cache kind's file under this name.
        let mut wrong_magic = good.clone();
        wrong_magic[..8].copy_from_slice(b"FDOTHER1");
        std::fs::write(&cache, &wrong_magic).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        // Wrong version.
        let mut wrong_version = good.clone();
        wrong_version[8..12].copy_from_slice(&99u32.to_le_bytes());
        std::fs::write(&cache, &wrong_version).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        // A zeroed count must not decode as an empty database.
        let mut zero_count = good.clone();
        zero_count[16..24].copy_from_slice(&0u64.to_le_bytes());
        std::fs::write(&cache, &zero_count).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        // Nor may a wildly large one.
        let mut huge_count = good.clone();
        huge_count[16..24].copy_from_slice(&u64::MAX.to_le_bytes());
        std::fs::write(&cache, &huge_count).unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());

        let _ = std::fs::remove_file(&src);
        let _ = std::fs::remove_file(&cache);
    }

    #[test]
    fn a_missing_cache_is_not_an_error() {
        let (src, cache) = temp_paths("missing");
        std::fs::write(&src, b"source bytes").unwrap();
        assert!(PodCache::<TestRecord>::map(&cache, &src).is_none());
        let _ = std::fs::remove_file(&src);
    }

    #[test]
    fn names_are_bounds_and_utf8_checked_at_resolution() {
        let (src, cache) = temp_paths("names");
        std::fs::write(&src, b"source bytes").unwrap();
        let records = vec![TestRecord { key: 1, name_offset: 0, name_len: 2, _pad: 0 }];
        let mapped = store_and_map(&cache, &src, &records, &[0xff, 0xfe]).unwrap();
        assert_eq!(mapped.name_at(0, 2), "");   // not UTF-8
        assert_eq!(mapped.name_at(0, 99), "");  // past the blob
        let _ = std::fs::remove_file(&src);
        let _ = std::fs::remove_file(&cache);
    }

    #[test]
    fn an_empty_cache_has_no_records() {
        let empty = PodCache::<TestRecord>::empty();
        assert!(empty.is_empty());
        assert_eq!(empty.len(), 0);
        assert_eq!(empty.records(), &[]);
        assert_eq!(empty.names_blob(), &[] as &[u8]);
        assert_eq!(empty.name_at(0, 3), "");
    }

    /// An unwritable path degrades to an 8-aligned in-memory image.
    #[test]
    fn an_unwritable_cache_path_falls_back_to_memory() {
        let (src, _) = temp_paths("readonly");
        std::fs::write(&src, b"source bytes").unwrap();
        let (records, names) = sample();
        let cache = "/proc/definitely-not-writable/fd.cache";
        let mapped = store_and_map(cache, &src, &records, &names)
            .expect("should fall back to an in-memory image");
        assert_eq!(mapped.records(), records.as_slice());
        assert_eq!(mapped.name_at(0, 5), "alpha");
        assert_eq!(mapped.records().as_ptr() as usize % align_of::<TestRecord>(), 0);
        let _ = std::fs::remove_file(&src);
    }
}
