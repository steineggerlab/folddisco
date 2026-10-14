# Folddisco 3.0 — update report

Printable version: [folddisco_3.0_report.pdf](folddisco_3.0_report.pdf). Benchmark scripts, raw
tables and figures live in [folddisco-analysis](https://github.com/steineggerlab/folddisco-analysis)
under `fd-branchbench/` (page source `report/index.html`, rendered by `scripts/make_report_pdf.py`);
measurement detail is in [feature_evaluation.md](feature_evaluation.md), branch state in
[UPDATE_REPORT.md](UPDATE_REPORT.md).

## What 3.0 adds

| feature | flag | what it does |
| --- | --- | --- |
| Sensitive search | `query --sensitive` (`--expand-radius`) | a residue pair may fall in neighbouring distance/angle bins, several features at once |
| Amino acid substitution | `query --aa-subst blosum62\|group\|size`, `:*` per residue | matches chemically similar residues, scored below exact ones |
| Confident hit list | `query --confident` | keeps only full, low-RMSD matches |
| Index-time expansion | `index --expand-radius/--expand-distance/--expand-angle/--aa-subst` | stores the neighbourhood in the index instead of expanding the query |
| Novelty screening | `query --novelty-mode` | one KNOWN/PARTIAL/NOVEL row per query with its evidence and the best hit’s residues |
| Multi-character chain IDs | `-q AA_250`, `10_250`, `--chain-sep` | mmCIF chains that are not a single letter |
| Deformation metrics | `drmsd`, `max_dist_deviation` | superposition-free deviation, as column, sort key and filter |
| Mapped lookup caches | automatic | `*.lookup.cache`, `*.fdcache` built on first use |

Indices built with 2.x are read unchanged.

## Benchmarks

| set | index | queries | answers | metric |
| --- | --- | --- | --- | --- |
| Motif | human AFDB, 23,391 structures | 8 commands over 4 motifs (3-23 residues) | zinc finger 761, MEROPS S01 124 accessions | precision/recall/F1, Sens@1FP |
| M-CSA | PDB-derived, 62,122 entries | 250 catalytic sites (3-21 residues), 492 held out | M-CSA homologues per site | Sens@1FP, average precision |
| Mutant M-CSA | same | the same 250 sites with one residue renamed to its closest BLOSUM62 alternative | unchanged | Sens@1FP |
| Motif-only index | 24,762 M-CSA motifs | 250 sites, published and mutant | same-entry motifs | Sens@1FP, top-1 |

Machine: 20 cores, local NVMe, warm cache. M-CSA timings run 5 queries in parallel with 4
threads each; motif timings are serial medians of 5 repeats.

## Results

Benchmark rerun at commit b907df5 (master 2a756d9); tables in `fd-branchbench/result/`.

### M-CSA, 250 catalytic sites

| config | Sens@1FP | Δ master | win / loss | no TP | s / query |
| --- | --- | --- | --- | --- | --- |
| master (2.x) | 0.4475 | — | — | 10 | 5.35 |
| 3.0 default | 0.4961 | +0.0486 | 133 / 42 | 1 | 2.65 |
| `--sensitive` | 0.5015 | +0.0540 | 134 / 42 | 1 | 3.02 |
| `--aa-subst blosum62` | 0.4764 | +0.0289 | 124 / 56 | 1 | 5.69 |
| `--aa-subst group` | 0.4734 | +0.0259 | 120 / 60 | 1 | 6.54 |
| `--aa-subst size` | 0.4662 | +0.0188 | 105 / 67 | 3 | 6.79 |
| `--sensitive --aa-subst blosum62` | 0.4807 | +0.0333 | 120 / 61 | 2 | 5.55 |

### Motif protocol (human proteome)

Same structures returned as 2.x on all eight commands (4/8 byte-identical, the rest reordered),
so set F1 is unchanged. `--sensitive` F1: zinc 4-residue 0.9421 → 0.9641, zinc 3-residue
0.9418 → 0.9577, triad prefilter 0.8831 → 0.9160, 23-residue segments 0.9204 → 0.9117.

### Mutated queries (one residue renamed per site)

| config | Sens@1FP | Δ vs mutant exact [95% CI] | win / loss |
| --- | --- | --- | --- |
| original query | 0.4961 | +0.194 [+0.162, +0.227] | 178 / 27 |
| mutant, exact | 0.3020 | — | — |
| mutant, `:*` | 0.4794 | +0.177 [+0.148, +0.208] | 172 / 22 |
| mutant, `--aa-subst blosum62` | 0.4546 | +0.153 [+0.124, +0.182] | 158 / 36 |
| mutant, `:*` + `--sensitive` | 0.4830 | +0.181 [+0.151, +0.214] | 176 / 23 |

### `--confident`

| M-CSA q250 | precision | recall | F1 | median hits | queries with a hit |
| --- | --- | --- | --- | --- | --- |
| default output | 0.058 | 0.751 | 0.084 | 6,000 | 100% |
| `--confident` | 0.757 | 0.475 | 0.508 | 55 | 92% |
| `--confident --sensitive` | 0.750 | 0.484 | 0.510 | 63 | 92% |

Motif protocol, F1 of one flag vs filters tuned per query: 0.944 vs 0.942 (zinc 4-residue),
0.930 vs 0.942 (zinc 3-residue), 0.911 vs 0.911 (triad), 0.525 vs 0.920 (23-residue segments,
where `--confident` is precision-first: 0.996 precision at 0.356 recall).

### `--confident` on the published motif protocols

| motif | filter | precision | recall | F1 | hits |
| --- | --- | --- | --- | --- | --- |
| Zinc finger, 4 res (`1G2F F207,F212,F225,F229`) | none | 0.492 | 0.974 | 0.654 | 7,425 |
| | protocol, by hand | 0.966 | 0.920 | 0.942 | 746 |
| | `--confident` | 0.961 | 0.928 | **0.944** | 4,304 |
| Zinc finger, 3 res (`1G2F F207,F225,F229`) | protocol | 0.959 | 0.925 | 0.942 | 756 |
| | `--confident` | 0.903 | 0.958 | 0.930 | 5,744 |
| Serine triad (`4CHA B57,B102,C195`) | none | 0.245 | 0.919 | 0.387 | 556 |
| | protocol | 0.964 | 0.863 | 0.911 | 113 |
| | `--confident` | 0.964 | 0.863 | 0.911 | 116 |
| Zinc, 23-res segments (`1G2F F204-215,F222-232`) | protocol | 0.972 | 0.874 | 0.920 | 697 |
| | `--confident` | 0.996 | 0.356 | 0.525 | 407 |

### Multi-character chain IDs

2.x panicked on mmCIF assemblies whose author chain IDs are not single letters
(`cif.rs: Chain name should be provided`). 3.0 keeps them (up to 8 characters) and accepts them
in `-q` with a separator:

```
folddisco query -i complexdir_folddisco -p 4V8S-assembly1.cif.gz -q AR_119,AR_364,AR_392
4V8S-assembly1.cif.gz   3   21.2189   0.0000   AR_119,AR_364,AR_392
4AYB-assembly1.cif.gz   3   15.4677   0.1150   B119,B364,B392
```

Output spells each chain in its own convention (`AR_119` vs `B119`, `_` for an unmatched query
residue); `--chain-sep` forces `A_81` everywhere. Missing separators are diagnosed.

### Ranking

Default order is now `match_score` (idf × coverage² × TM-score) per match and `structure_score`
(matched² × √idf / (1 + RMSD)) per structure, selected from 281 lexicographic orders and 69
composite scores and validated on 492 held-out M-CSA queries: +0.028 to +0.054 Sens@1FP and
+0.034 to +0.050 average precision over 2.x, every 95% interval excluding zero (§15).

### Index-time expansion

On a motif-only index of 24,762 M-CSA sites, index-side expansion scores 0.025–0.037 below
query-side expansion, returns a different ranking on every query (motif overlap 43–58%), and
grows the index 3.5× (radius 1) to 32× (radius 1 + blosum62).

## Compatibility

- Index files are unchanged; 2.x indices need no rebuild.
- `--nonrigid` is now `--sensitive`.
- Default ordering changed (see below), so hit lists come back in a different order than 2.x.
  On the eight motif commands the branch returns exactly the same structures as
  `origin/master`; the four `--per-structure` commands return them in a new order.
- The residue that ends a chain is no longer labelled with the next chain's ID, which changes
  matched-residue labels on multi-chain entries.
- `-d`/`-a` with several values now use only the widest one.

## Limits

- Sens@1FP is a k=1 metric: single accessions move it. Paired win/loss counts and average
  precision are reported alongside.
- M-CSA answer sets are homologous entries, so a partial match can be a true answer; the zinc
  and serine sets are defined by exact residue identity.
- `--confident` past 12 residues rests on four queries only.
- Index-time expansion was measured on one motif-only index.
- Single machine, one hash type (`PDBTrRosetta`, default binning).

## Reproduce

```bash
cd fd-branchbench   # in the folddisco-analysis checkout
REPEATS=5 WARMUP=1 ./scripts/run_all.sh 250     # motif, M-CSA, mutant, --confident, figures
defaults/scripts/collect.sh <binary>            # raw output for sort/filter selection
defaults/scripts/analyze.py && defaults/scripts/plot.py
index_expansion/scripts/run_queries.sh          # motif-only index study
```
