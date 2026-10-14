# Folddisco Index Structure

`folddisco index -i <prefix>` creates a matched index set. Rebuild the index
after filtering or reordering structures; do not edit one file independently.

## Folddisco Files

| File | Purpose |
| --- | --- |
| `<prefix>` | Compressed posting lists for geometric hashes. |
| `<prefix>.offset` | Maps a hash to its byte range in `<prefix>`. |
| `<prefix>.lookup` | Target names, metadata, and external database keys. |
| `<prefix>.type` | TOML hashing and input configuration. |
| `<prefix>.lookup.cache` | Regenerable memory-mapped cache of `.lookup`. |

The posting-list file stores sorted target IDs. The first ID in each list is
absolute; later IDs are deltas, encoded as unsigned 7-bit little-endian
varints. `<prefix>.offset` has the binary layout:

```text
hash_count: usize
hashes:     [u32; hash_count]
offsets:    [usize; hash_count + 1]
```

`hashes` is sorted. `offsets[i]..offsets[i + 1]` is the posting-list byte range
for `hashes[i]`. The native-width `usize` fields make this an
architecture-dependent format.

## Lookup IDs

`.lookup` is tab-separated UTF-8 text:

```text
internal_id<TAB>name<TAB>residue_count<TAB>mean_plddt<TAB>external_db_key
```

- `internal_id` is the posting-list ID and must be the dense lookup row number:
   `0, 1, ..., n - 1`.
- `name`, `residue_count`, and `mean_plddt` describe the indexed structure.
- `external_db_key` is optional for legacy four-column files, where it defaults
   to `internal_id`.

The first column cannot safely be sparse. The fifth column may be sparse when
it refers to a Foldcomp/MMseqs2-style database.

`.lookup.cache` is rebuilt automatically when `.lookup` changes or the cache
is invalid. It is disposable and should not be treated as source data.

## Foldcomp Support

For a Foldcomp input database `<db>`, `.type` records `foldcomp_db = "<db>"`.
Folddisco writes the Foldcomp key to column 5 of its `.lookup` and uses it to
retrieve matching structures.

| File | Purpose |
| --- | --- |
| `<db>` | Foldcomp structure data. |
| `<db>.lookup` | `key<TAB>name<TAB>file` records. |
| `<db>.index` | `key<TAB>offset<TAB>length` records. |
| `<db>.lookup.fdcache` | Regenerable mapped, key-sorted cache of `.lookup`. |
| `<db>.index.fdcache` | Regenerable mapped, key-sorted cache of `.index`. |

Foldcomp keys may be sparse. Its lookup and index records are sorted and
resolved by key, so their keys must remain a matching pair. The `.fdcache`
files are separate from Folddisco's `.lookup.cache` and may be deleted to
force regeneration.

Changing hash type, binning, multiple bins, or index-time expansion in
`<prefix>.type` requires rebuilding the Folddisco index.