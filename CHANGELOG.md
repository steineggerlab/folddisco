# Changelog

## 3.0.0

Indices built with 2.x are read unchanged.

### New
- `query --sensitive` (= `--expand-radius 2`) and `--expand-radius`: residue pairs may fall
  in neighbouring distance/angle bins, several features at once. `-d`/`-a` are sub-stepped.
- `query --aa-subst blosum62|group|size` and per-residue `:*`: substitution schemes.
  Substituted residues score below exact ones, so exact matches rank first.
- `index --expand-radius/--expand-distance/--expand-angle/--aa-subst`: index-time expansion.
- `query --confident`: keep only confident, full matches (≥ 80% of query residues within 1 Å;
  coverage only past 12 residues); filters given explicitly are left alone.
- `query --novelty-mode` with `--novelty-coverage`/`--novelty-rmsd`: one KNOWN/PARTIAL/NOVEL
  row per query, with candidates, coverage, best hit, RMSD and the best hit's matched residues.
- Multi-character and numeric chain IDs (`AA_250`, `10_250`); `--chain-sep`.
- Superposition-free `drmsd` and `max_dist_deviation` (columns, sort keys, `--drmsd`).
- Memory-mapped caches for index and Foldcomp lookups (`*.lookup.cache`, `*.fdcache`),
  built on first use.

### Changed
- Default sort: `match_score:desc,rmsd:asc` per match (was `idf:desc,rmsd:asc`) and
  `structure_score:desc,min_rmsd:asc` per structure (was `idf:desc,min_rmsd:asc`).
  `match_score` (`idf` × coverage² × TM-score), `structure_score` (matched² × √`idf` / (1 + RMSD))
  and `coverage_idf` are new sort keys and output columns.
- Matching prefilters target residues by amino acid code for large queries (same output, faster).
- Per-match IDF counts edges from neighbouring bins and substitutions only between the matched
  residues, so large components of similar residues no longer rank first.
- A residue that ends a chain is no longer labelled with the next chain's ID.
- `-d`/`-a` with several values use only the widest one.
- `-q` ranges may repeat the chain (`F204-F215`); malformed queries are reported, not panicked on.
- Help states the actual default sort, `idf:desc,rmsd:asc`.

### Fixed
- Substituted residues (`:H`) were found by the index but dropped at matching.
