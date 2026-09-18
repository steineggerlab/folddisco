# folddisco — update report

State of `feature-integration`. Measurements live in `feature_evaluation.md`; user docs in `README.md`.

## 1. This round

### Default filtering and sorting (latest)
- `--confident`: coverage ratio 0.8 on both filter stages (all residues for a 3-4 residue motif,
  one may be missing from five up) and RMSD ≤ 1.0 Å for queries of at most 12 residues; longer
  queries keep coverage only. Explicit filters win, `--skip-match` gets coverage only.
  `confident_filters` in `query_pdb.rs`.
- Default sort is now `match_score:desc,rmsd:asc` per match (`idf` × coverage² × TM-score) and
  `structure_score:desc,min_rmsd:asc` per structure (matched² × √`idf` / (1 + RMSD)); both are
  sort keys and output columns, as is `coverage_idf`.
- Selected from 281 lexicographic orders and 69 composite scores on re-sorted raw output
  (M-CSA q250 + motif), validated on 492 held-out M-CSA queries; see `feature_evaluation.md` §15.

### Substitution scoring (latest)
- The query map stores a weight per hash: 1 for the query's own residues, `SUBSTITUTION_WEIGHT`
  (0.75) per substituted side. A substituted hash's IDF is capped at the observed pair's.
- `count_query`: substituted hits add `weight · IDF` once per query edge (best hash), and only
  to targets without an exact hit on that edge.
- Per-match IDF: substituted and neighbouring-bin edges count only between the matched residues
  (`QueryHash::observed` marks the observed geometry); substituted ones are also weighted.
  Candidate scoring keeps summing neighbour hits: capping them there cost zinc-finger recall.
- An exact hash is no longer shadowed by another pair's substitution.
- Matching: above 200 query hashes the residue prefilter uses amino acid codes from
  `aa_dist_map` instead of scanning every pair; output is byte-identical.
- The per-pair 4096-hash cap never binds at default tolerances (radius 2: 50 bins × at most 36
  residue pairs); the earlier loss came from scoring, not the cap. Measurements: `feature_evaluation.md` §14.
- Paths without substitution are unchanged (rank-identical on M-CSA q60).

### `--nonrigid` renamed `--sensitive`
- Same preset (`--expand-radius 2`); the name describes the effect, not a flexible alignment.

### Novelty mode
- `--novelty-mode` prints one verdict row per query:
  `query_id, verdict, candidates, hits, index_coverage, best_hit, best_coverage, best_rmsd, best_residues, query_residues`.
- `verdict` is `KNOWN`, `PARTIAL`, `NOVEL` or `NO_HASHES`, gated by `--novelty-coverage` [0.8]
  and `--novelty-rmsd` [2.0]; `candidates`/`index_coverage` are taken before filters, and
  `best_residues` gives the residues the best hit matched. `--header` prints the column names once.

### Amino acid substitution schemes
- `src/controller/substitution.rs`: `blosum62` (positive score), `group` (RHK, DE, NQST, FWY,
  AVLIMC, GP), `size` (IMGT volume classes GAS, CDPNT, QEHV, MILKR, FWY).
- `query --aa-subst <MODE>` applies the scheme to every residue without an explicit `:ALT`;
  `:*` applies it to one residue (default `blosum62`), and `:Q*` adds Q to the scheme's set.
  Resolved against the residue observed in the query structure.
- Hash types without residue identity in the hash (Hybrid, TertiaryInteraction) ignore it with a warning.
- **Bug fixed:** matching only accepted target residue pairs whose amino acids were in the
  query's observed `aa_dist_map`, so any substituted residue — including the existing `:H`
  syntax — passed the index lookup and was then dropped at matching. Substituted pairs are now
  registered. Test: `substituted_residue_is_matched_not_only_looked_up` (4CHA with D102N).

### Index-time expansion
- `index --expand-radius <INT> --expand-distance <F> --expand-angle <F> --aa-subst <MODE>`
  also indexes each target pair's neighbour bins and substituted pairs. Off by default; for
  small databases, since the index grows several-fold.
- Stored in `<index>.type` (`expand_radius`, `expand_distance`, `expand_angle`, `aa_subst`);
  older `.type` files without these keys read as a plain index.
- Query and index share `expand::for_each_expanded_feature`, so both sides reach the same hashes.
- On an expanded index, query lookup uses exact hashes unless `-d`, `-a`, `--expand-radius` or
  `--aa-subst` is given. Matching uses `IndexExpansion::matching_tolerance` (offsets and radii
  summed) and the index scheme applied to every residue, so candidates returned by the index
  can still be matched.

### Tests and fixes
- The lib tests did not compile after the chain-ID merge (23 errors: `u8` vs `ChainId`, the old
  `make_query_map` signature, a test of the removed `_add_shifted_hashes`). Fixed.
- The merge had also dropped `F204-F215` (chain repeated on a range end), which the README
  documents; restored, with a different chain rejected.
- Query help and README now state the real default sort, `idf:desc,rmsd:asc`.

## 2. Sensitive search

`--expand-radius` (default 1) and `--sensitive` (= 2), `src/controller/expand.rs`.

- `PDBTrRosetta` hashes `[aa1, aa2, ca_dist, cb_dist, ca_cb_angle, theta1, theta2]` into 30 bits:
  residue identities exact; two distances at 1.2 Å per bin (16 bins over 2–20 Å); three angles as
  sin and cos, 4 bins each. Lookup is exact on the packed `u32`.
- Distances and angles form one list of tolerant dimensions. Each gets sub-stepped offsets
  (`ceil(tol / bin_width)` steps, max 8); only the widest `-d`/`-a` is used. Levels `1..=radius`
  choose that many distinct dimensions: a Hamming ball, not the full product.
- Perturbed torsions wrap at ±π and `acos` angles reflect; distances are not clamped.
- Hashes per pair at defaults (5 dimensions, 2 offsets): radius 0 → 1, 1 → 10, 2 → 50, 3 → 130,
  times `1 + substitution variants`, capped at `MAX_HASHES_PER_PAIR = 4096` attempts.
- Candidates must reproduce a query hash exactly on re-hashing the target, gated by `--ca-distance`.
- All quality filters default to off, and expansion only wins with `--max-node` (F1 0.9421 →
  0.9641 with it, 0.9265 → 0.9226 without). Cost +11% (human proteome) to ~73% (M-CSA).
- Removed: `--enm-sample` and the rare-hash IDF filter. `is_primary` survives with no production
  consumer. Older reports that say `--nonrigid` predate both removals.

Open concerns:
1. `MAX_SUBSTEPS = 8` clamps beyond 8 bins of tolerance (`-d > 9.6 Å`), reintroducing bin gaps.
2. Under the 4096 cap, surviving joint dimensions depend on order (distance-first), not geometry.
3. First-writer-wins dedup in `insert_binned_hash` lets pair A's variant shadow pair B's observed
   hash in residue mapping; RMSD is the only backstop.
4. Out-of-window Cβ distances can alias into the residue-identity bits (pre-existing).

## 3. IDF scoring

Per-structure `idf`: `nres^(-lp) · Σ log2(N / df(h))` over matched query hashes (`count_query.rs`),
a plain sum; `--length-penalty` `lp` default 0.5. Per-match `idf` (`calculate_subgraph_idf`) sums
over **every edge of the connected component** the match came from, not the matched residues'
edges, while `node_count` counts matched residues. Components are kept at `len() >= node_count`
with no upper bound, so repetitive structures produce very large scores (the EF-hand
852,272-vs-110.688 case; mechanism read from code, value not reproduced here).

Also:
- Geometric variants inherit the observed hash's IDF in `query.rs`; `count_query` recomputes per
  hash. Substituted variants keep the lower of their own and the observed IDF, weighted (§1).
- Per-structure `idf_sum` adds once per matching hash while `edge_count` counts edges, so many
  neighbour bins of one edge can outrank one bin of every edge.
- `--score` compares the per-structure sum before matching and the subgraph IDF per match.
- `--freq-filter` is off by default, but help shows `[0.0]`, and passing 0 skips every hash.
- E-value (`result.rs`) is an unvalidated fit of the subgraph IDF, only via `--format-output e_value`.
- At PDB scale raw IDF mis-ranks; gate with `--max-node`/`--covered-node` or sort by `node_count`.

What to fix, in order:
1. Score the subgraph induced by the reported residues, bounding it by `C(node_count, 2) · log2(N)`.
   Done for substituted and neighbouring-bin edges (§1); exact edges still sum over the component.
2. Use one IDF convention (per-hash recompute or inherited), not both.
3. Rename one of the two `idf` columns.

## 4. Lookup caches

`src/utils/pod_cache.rs`. Index lookup and both Foldcomp tables are memory-mapped caches whose
on-disk layout is the in-memory layout.

```text
offset  0   magic         [u8; 8]   FDLOOKUP | FDFCLKUP | FDFCIDXT
offset  8   version       u32
offset 12   pad           u32
offset 16   count         u64
offset 24   source_len    u64
offset 32   source_mtime  u64       nanoseconds
offset 40   names_len     u64
offset 48   records       count * size_of::<R>()   (8-aligned, cast as &[R])
            names         names_len bytes           (sliced, never copied)
```

| record | size | fields |
|---|---|---|
| `LookupRecord` | 40 B | `id u64`, `nres u64`, `db_key u64`, `name_offset u64`, `plddt f32`, `name_len u32` |
| `FoldcompLookupRecord` | 24 B | `key u64`, `name_offset u64`, `name_len u32`, `_pad u32` |
| `FoldcompIndexRecord` | 24 B | `key u64`, `offset u64`, `length u64` |

- Soundness: `unsafe trait CacheRecord` requires `#[repr(C)]`, no padding, align ≤ 8, size a
  multiple of 8, valid for any bit pattern. Alignment is asserted at load. Read-only directories
  build the image in memory as `Vec<u64>` for alignment.
- Validation: magic, version, source length and mtime, exact file size, then one parallel pass
  (non-empty names in bounds, `id < len`, strictly ascending keys, non-zero entry lengths). Any
  failure reparses the text. No checksum and no temp-and-rename, both relying on deterministic bytes.
- Foldcomp tables are sorted by key at build; name → key is a parallel scan (≤2 per query).

| | before | after |
|---|---|---|
| `afdb50_v4` lookup (53.7 M entries) load | 2.98 s, 8.55 GB anon | **17.5–17.9 ms**, 2.10 GB page cache |
| `afdb_uniprot_v6` Foldcomp lookup (241 M) load | 6.14 s every query, 27.9 GB | **45.8–50.3 ms**, 5.66 GB page cache |
| Foldcomp key → name | — | ~4.4–5.1 µs |
| Foldcomp name → key | — | ~98–102 ms |

To decide:
1. Foldcomp caches are 11.7 GB for a 10.78 GB lookup; a `u32` key would bring it to ~8 GB.
2. The validation pass is most of the load time; dropping it gives ~20 µs at the cost of
   torn-write protection.

## 5. Multi-character chain IDs

- `ChainId`: up to 8 printable ASCII bytes inline, `Copy`. `Atom::chain` stays one byte because
  it mirrors Foldcomp's C `atom_t`.
- The CIF reader now keeps multi-character `auth_asym_id` instead of falling back to `label_asym_id`.
- Output keeps `A21,A23`; a field switches to `AA_250` / `10_250` only when a chain in it is
  multi-character or numeric. `--chain-sep` forces `_`. `-q` accepts both spellings.
- Default output and index files are byte-identical to before; indices need no rebuild.
- Fixed: a residue ending a chain was filed under the next chain (e.g. 4CHA `A11` reported as
  `B11`). Changes results on multi-chain structures, correctly. The same one-residue shift in
  `b_factors` is left as is.

## 6. Tests

Run `cargo test --release`.
