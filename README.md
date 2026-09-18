# Folddisco

<p align="center">
<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/steineggerlab/folddisco/master/.github/img/folddisco_logo_with_light.png">
<img src="https://raw.githubusercontent.com/steineggerlab/folddisco/master/.github/img/folddisco_logo_without_light.png" max-height="300px" height="300" display="block" margin-left="auto" margin-right="auto" display="block"/>
</picture>

Folddisco is tool for searching discontinuous motifs in protein structures.
It is designed to handle large-scale protein databases with efficiency, enabling the detection of structural motifs across thousands of proteomes or millions of structures.

## Publications
[Kim H, Kim RS, Mirdita M, Yoon J, Steinegger M. Structural motif search across the protein-universe with Folddisco. Nature Biotechnology, (2026)](https://www.nature.com/articles/s41587-026-03162-9)

What 3.0 changed: [docs/v3_report.md](docs/v3_report.md) ([PDF](docs/folddisco_3.0_report.pdf)).

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/folddisco.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/folddisco) [![Github All Releases](https://img.shields.io/github/downloads/steineggerlab/folddisco/total.svg)](https://github.com/steineggerlab/folddisco/releases/latest) 

## Webserver 
Search protein structures motifs against the [AlphaFoldDB](https://alphafold.ebi.ac.uk/) and [PDB](https://www.rcsb.org/) in seconds using the Folddisco webserver ([code](https://github.com/soedinglab/mmseqs2-app)): [search.foldseek.com/folddisco](https://search.foldseek.com/folddisco) 🚀

## Installation
```bash
# Install from Bioconda
conda create -n folddisco -c conda-forge -c bioconda folddisco

# Install through docker
docker pull ghcr.io/steineggerlab/folddisco:master

# Precompiled binary for Linux x86-64
wget https://mmseqs.com/folddisco/folddisco-linux-x86_64.tar.gz; tar xvfz folddisco-linux-x86_64.tar.gz; export PATH=$(pwd)/folddisco/bin/:$PATH

# Precompiled binary for Linux ARM64
wget https://mmseqs.com/folddisco/folddisco-linux-arm64.tar.gz; tar xvfz folddisco-linux-arm64.tar.gz; export PATH=$(pwd)/folddisco/bin/:$PATH

# macOS (universal, works on Apple Silicon and Intel Macs)
wget https://mmseqs.com/folddisco/folddisco-macos-universal.tar.gz; tar xvfz folddisco-macos-universal.tar.gz; export PATH=$(pwd)/folddisco/bin/:$PATH
```

**Compile from source**

Compiling from source requires the Rust toolchain (Cargo). Installation instructions are available [here](https://www.rust-lang.org/tools/install).

```bash
git clone https://github.com/steineggerlab/folddisco.git
cd folddisco
cargo install --features foldcomp --path .
```

## Quick start
Folddisco queries a database of precomputed geometric hashes computed from structures. 

### Download pre-build database 
You can download the pre-built human proteome index and use it to search for a common motif, like a zinc finger.

This example is fully self-contained. You can copy and paste the entire block into your terminal.

```bash
# Download human proteome index. Use wget or aria2 to download the index.
cd index
aria2c https://opendata.mmseqs.org/folddisco/h_sapiens_folddisco.tar.lz4

# Extract the index
lz4 -dc h_sapiens_folddisco.tar.lz4 | tar -xvf -
cd ..
```

#### Pre-built Indices

Download pre-built index files:
- [Human proteome](https://opendata.mmseqs.org/folddisco/h_sapiens_folddisco.tar.lz4)
- [E. coli proteome](https://opendata.mmseqs.org/folddisco/e_coli_folddisco.tar.lz4)
- [AFDB proteome of 16 model organisms](https://opendata.mmseqs.org/folddisco/afdb_proteome_v4_folddisco.tar.lz4)
- [Swiss-Prot](https://opendata.mmseqs.org/folddisco/afdb_swissprot_v4_folddisco.tar.lz4)
- [AFDB50](https://opendata.mmseqs.org/folddisco/afdb50_v4_folddisco.tar.lz4) 
- [ESM30](https://opendata.mmseqs.org/folddisco/highquality_clust30_folddisco.tar.lz4)
- [PDB](https://opendata.mmseqs.org/folddisco/pdb_folddisco.tar.lz4)
- To get the old version of Folddisco indices, please **visit** https://opendata.mmseqs.org/folddisco/
  - `*.tar.gz` indices are legacy indices (version 1.0), which are not compatible with version 2.0 or later.
    Please use `*.tar.lz4` indices; version 3.0 reads them unchanged.
  - **AFDB50** (`afdb50_v4_folddisco*` + `afdb50_v4*`)
  - **ESM30** (`highquality_clust30_folddisco*` + `highquality_clust30*`)

### Build an custom index 
The command below will read all PDB or mmCIF from `serine_peptidases` folder and generate an index `serine_peptidases_folddisco`.
```bash
folddisco index -p data/serine_peptidases -i index/serine_peptidases_folddisco
```

### Querying a Single Motif
To search for a specific structural motif, you'll use three main flags:
-   **`-p`**: Provides the query protein's structure file (PDB/mmCIF).
-   **`-q`**: Specifies the comma-separated list of residues that form your motif.
-   **`-i`**: Points to the target database index you want to search against.

If you omit the **`-q`** flag, `folddisco` defaults to a "whole structure" search. It will find all possible motifs from your entire query protein and search for them in the index.

```bash
# Search for the catalytic triad from 4CHA.pdb against the indexed peptidases.
folddisco query -i index/serine_peptidases_folddisco -p query/4CHA.pdb -q B57,B102,C195
```
#### Residue & motif syntax
We allow to customize the query motif using some motif syntax.
* **Residues:** `B57` = chain `B`, residue number `57`. Ranges are inclusive and may
  repeat the chain on the end: `1-10`, `F204-215` and `F204-F215` all work. A range
  cannot span two chains.
* **Lists:** comma-separated: `B57,B102,C195`.
* **Multi-character and numeric chain IDs:** large mmCIF entries use chain IDs
  such as `AA` or `10`, which cannot be pasted straight onto a residue number
  (`AA250` and `10250` are unreadable). Separate them with `_`:
  `-q 10_250,AA_312,AA_318`. The separator is optional for ordinary
  single-letter chains, so `B57` and `B_57` mean the same thing.
* **Substitutions:** `:<ALT>` allows alternatives:
  * Single amino acid: `164:H`
  * Set: `247:ND` (Asp or Asn)
  * Wildcard/categories:
    * `X`: any amino acid
    * `p`: positively charged (Arg, His, Lys)
    * `n`: negatively charged (Asp, Glu)
    * `h`: polar (Asn, Gln, Ser, Thr, Tyr)
    * `b`: hydrophobic (Ala, Cys, Gly, Ile, Leu, Met, Phe, Pro, Val)
    * `a`: aromatic (His, Phe, Trp, Tyr)
  * Scheme: `102:*` uses the `--aa-subst` scheme (default `blosum62`) for the observed residue

### Searching Multiple Motifs (Batch Mode)
To search for many motifs at once, you can provide a single query file to the **`-q`** flag (and omit the `-p` flag).

This file must be a **tab-separated** text file with these columns:
1.  **Column 1:** Path to the query structure (PDB/mmCIF).
2.  **Column 2:** Comma-separated list of motif residues.
3.  **Column 3:** (Optional) path to the output file (default: `stdout`).

```bash
# Search a zinc finger motif against pre-downloaded human proteome (see Download pre-build database)
folddisco query -i index/h_sapiens_folddisco -q query/serine_peptidase.txt
```

## Commands

### Usage of Query Module
```bash
folddisco query -i <INDEX> -p <QUERY_PDB> [-q <QUERY_RESIDUES> -d <DISTANCE_THRESHOLD> -a <ANGLE_THRESHOLD> --skip-match -t <THREADS>]
```

**Important parameter:**
- `-d`: Distance tolerance in Å, increase sensitivity during the prefilter (default: 0.5)
- `-a`: Angle tolerance in degrees, increase sensitivity during the prefilter (default: 5)
- `--sensitive`: Wider, slower search for deformed motifs (see [Sensitive search](#sensitive-search))
- `--expand-radius`: How many geometric features may fall in a neighbouring bin at once (default: 1)
- `--confident`: Keep only confident, full matches (see [Confident hits](#confident-hits))
- `--aa-subst`: Substitute every query residue by scheme (see [Amino acid substitution](#amino-acid-substitution))
- `--novelty-mode`: One KNOWN/PARTIAL/NOVEL row per query instead of a hit list (see [Novelty screening](#novelty-screening))
- `--skip-match`: Skips residue matching and RMSD calculation (prefilter only, much faster with same ranking)
- `--top`: Only report top N hits from the prefilter (controls speed and size of result)
- `-t`: Threads used for search
- `-v`: Verbose output

#### Example Querying
```bash
# Search with default settings (sorted by IDF, then RMSD)
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6
folddisco query -p query/1G2F.pdb -q F207,F212,F225,F229 -i index/h_sapiens_folddisco -d 0.5 -a 5 -t 6
folddisco query -p query/1LAP.pdb -q 250,255,273,332,334 -i index/h_sapiens_folddisco --skip-match -t 6 # Skip residue matching

# Query file given as separate text file
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 -d 0.5 -a 5

# Querying a whole structure
folddisco query -i index/h_sapiens_folddisco -p query/1G2F.pdb -t 6 --skip-match
# For a long query, low `--sampling-ratio` can be used to speed up the search
folddisco query -i index/h_sapiens_folddisco -p query/1G2F.pdb -t 6  --skip-match --sampling-ratio 0.3

# Using a query file with distance and angle thresholds
folddisco query -i index/h_sapiens_folddisco -q query/knottin.txt -d 0.5 -a 5 --skip-match -t 6

# Query with amino-acid substitutions and range. 
# Alternative amino acids can be given after colon. 
# X: substitute to any amino acid, p: positive-charged, n: negative-charged, h: hydrophilic, b: hydrophobic, a: aromatic
# Here's enolase query with 3 substitutions; Allow His at 164, Asp & Asn at 247, and His at 297. (Install e_coli_folddisco index first)
folddisco query -p query/2MNR.pdb -q 164:H,195,221,247:ND,297:H -i index/e_coli_folddisco -d 0.5 -a 5 --top 10 --header --per-structure
# Range can be given with dash. This will query first 10 residues and 11th residue with subsitution to any amino acid.
folddisco query -p query/4CHA.pdb -q 1-10,11:X -i index/h_sapiens_folddisco -t 6 --serial-index

# Advanced query with filtering and sorting
## Based on connected node and rmsd
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --connected-node 0.75 --rmsd 1.0

## Coverage based filtering & top N filtering without residue matching
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --covered-node 3 --top 1000 --per-structure --skip-match

# Print top 100 structures with sorting by score
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --top 100 --per-structure --sort-by idf
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --covered-node 4 --top 100 --sort-by idf --per-structure --skip-match

# Comprehensive filtering with multiple criteria
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 -d 0.5 -a 10.0 --ca-distance 1.0 --covered-node-ratio 0.3 --max-node-ratio 0.35 --rmsd 5.0 --tm-score 0.2 --gdt-ts 0.25 --gdt-ha 0.15 --chamfer 5.5 --hausdorff 12.0 --sort-by node_count,gdt_ts,rmsd,idf --format-output tid,node_count,gdt_ts,rmsd,idf,matching_residues,query_residues
```

### Sensitive search

A residue pair whose distance or angle drifts across a bin boundary gets a different hash.
`--expand-radius` sets how many features of a pair may sit in a neighbouring bin at once
(default 1, 0 = observed bins only); `--sensitive` is `--expand-radius 2`.

```bash
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --sensitive --max-node 3
```

- Pair it with `--max-node <n_residues>`: without it the extra recall costs more precision than it gains.
- Not for long segment queries, which are already saturated.
- `-d`/`-a` move one feature further; `--sensitive` lets more features move together.
- Rank deformed motifs by `drmsd` (superposition-free); `drmsd` and `max_dist_deviation` work in
  `--format-output`, `--sort-by` and as filters.

F1 on the human proteome (details in [feature_evaluation.md](docs/feature_evaluation.md)):

| query | default | `--sensitive` |
| --- | --- | --- |
| 4-residue zinc finger, matched | 0.9421 | **0.9641** |
| 3-residue zinc finger, matched | 0.9418 | **0.9577** |
| Ser-His-Asp triad, prefilter (MEROPS S01) | 0.8831 | **0.9160** |
| 23-residue two-segment query, matched | **0.9204** | 0.9117 |

### Confident hits

By default a query returns every structure that matched any part of the motif, partial matches
included. `--confident` keeps only the confident, full ones: at least 80% of the query residues
matched — all of them for a 3-4 residue motif, one may be missing from five residues up — within
1 Å RMSD. Past 12 residues only the coverage is required, because such matches are assembled
from several parts and run to several Å while the coverage alone is already precise. Filters you
give explicitly are left alone, so `--confident --rmsd 0.5` tightens only the RMSD. With
`--skip-match` there is no superposition, so only the coverage applies.

```bash
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --confident
```

On the M-CSA benchmark (250 catalytic sites, 3-21 residues) it raises precision from 0.06 to
0.76 and cuts the median hit list from 1,804 to 21, keeping 0.48 of the answers (0.74 unfiltered);
92% of queries still return something. See [feature_evaluation.md](docs/feature_evaluation.md) §15.

### Amino acid substitution

`--aa-subst <MODE>` lets every query residue without an explicit `:ALT` match similar residues;
`:*` does the same for one residue. Substitutions compose with the geometric tolerance.
A substituted residue scores 0.75 of an exact one and never above the query's own residue pair,
so exact matches rank first. Prefer `:*` on the residues expected to vary.

| mode | alternatives for the observed residue |
| --- | --- |
| `blosum62` | positive BLOSUM62 score (default for `:*`) |
| `group` | same class: RHK, DE, NQST, FWY, AVLIMC, GP |
| `size` | same IMGT side-chain volume class: GAS, CDPNT, QEHV, MILKR, FWY |

```bash
folddisco query -p query/4CHA.pdb -q B57,B102:*,C195 -i index/h_sapiens_folddisco -t 6
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --aa-subst group
```

### Novelty screening

`--novelty-mode` prints one verdict row per query instead of a hit list, with the evidence
behind the verdict.

```bash
folddisco query -q designs.txt -i index/pdb_folddisco -t 6 --novelty-mode --header
```

| column | meaning |
| --- | --- |
| `verdict` | `KNOWN` (best hit covers ≥ `--novelty-coverage` of the query within `--novelty-rmsd`), `PARTIAL` (covered but under either threshold), `NOVEL` (nothing covered), `NO_HASHES` (residues too far apart to search) |
| `candidates` | structures the index returned, before filters |
| `hits` | structures left after filters and matching |
| `index_coverage` | best hash-level residue coverage among candidates |
| `best_hit`, `best_coverage`, `best_rmsd` | highest-coverage hit after filters; RMSD is `NA` with `--skip-match` |
| `best_residues` | residues that hit matched, `_` where the query residue went unmatched |

Defaults are `--novelty-coverage 0.8` and `--novelty-rmsd 2.0`. A 3-4 residue motif covers
1.0 against something almost anywhere, so tighten `--novelty-rmsd` when screening short motifs.

Rows append to `-o`, so batch queries can share one file; rerunning replaces it. Filters
such as a high `--max-node` drop partial matches, which are often the most useful evidence.
A residue listed twice counts twice in coverage denominators.

### Lookup cache

The first load of an index writes `<index>.lookup.cache`; later loads decode it instead of
parsing the text lookup. It is validated against the lookup's size and mtime, falls back to
parsing when stale, and is safe to delete.

### Usage of Index Module
```bash
folddisco index -p <PDB_DIR|FOLDCOMP_DB> -i <INDEX_PATH> -t <THREADS> [-d <DISTANCE_BINS> -a <ANGLE_BINS> -y <HASH_TYPE>]
```

**Important parameter:**
- `-d`: Number of distance bins (default for `default` type: 16)
- `-a`: Number of angle bins (default for `default` type: 4)
- `-y`, `--type`: Hash type: `default`, `pdb`, `trrosetta`, `ppf`, `3di`
- `-t`: Threads
- `-v`: Verbose output

#### Examples
```bash
# Default indexing; a directory or a Foldcomp database
folddisco index -p h_sapiens -i index/h_sapiens_folddisco -t 12

# Custom hash type and bins
folddisco index -p h_sapiens -i index/h_sapiens_pdbtype -t 12 -y pdb -d 8 -a 3
```

### Index-time expansion

For small databases, the query expansion can be stored in the index instead:
each target pair is also indexed under its neighbouring bins and substituted residue pairs.
The index grows several-fold, so this is off by default.

```bash
folddisco index -p data/serine_peptidases -i index/serine_expanded -t 12 --expand-radius 1 --aa-subst blosum62
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/serine_expanded
```

`--expand-radius`, `--expand-distance` (0.5 Å) and `--expand-angle` (5°) mirror the query
options. Settings are stored in `<index>.type`. On such an index, query lookup uses exact
hashes unless `-d`, `-a`, `--expand-radius` or `--aa-subst` are given, and residue matching
applies the combined expansion.

## Output
### Match Result
Default output which prints out one matching motif per line
```
tid	node_count	idf	rmsd	matching_residues	query_residues
data/serine_peptidases/4cha.pdb	3	8.7616	0.0000	B57,B102,C195	B57,B102,C195
data/serine_peptidases/4cha.pdb	3	8.7616	0.0874	F57,F102,G195	B57,B102,C195
data/serine_peptidases/1pq5.pdb	3	4.1178	0.2609	A56,A99,A195	B57,B102,C195
data/serine_peptidases/1ju3.pdb	2	1.4739	0.7792	_,A223,A234	B57,B102,C195
data/serine_peptidases/1l7a.pdb	2	1.4739	0.7883	_,A146,A127	B57,B102,C195
data/serine_peptidases/1l7a.pdb	2	1.4739	0.8078	_,B146,B127	B57,B102,C195
data/serine_peptidases/1azw.pdb	2	4.6439	0.9234	A179,_,B176	B57,B102,C195
```
- `tid`: Identifier of the target protein structure
- `node_count`: Number of nodes in the match
- `idf`: Inverse document frequency score of matched structure
- `rmsd`: Root mean square deviation
- `matching_residues`: Residue indices in the match (comma-separated, _ for no match)
- `query_residues`: Residue indices in the query (comma-separated)

Residues are written as chain + residue number (`B57`), which is what the `-q`
grammar accepts. When a chain ID would make that unreadable -- a multi-character
`AA` or a numeric `10` -- the whole field switches to `AA_250` / `10_250`
instead, which `-q` also accepts. Pass `--chain-sep` to get the separated form
for every structure, if you would rather parse one fixed format.

### Structure Result
Output with one structure per line (`--per-structure`)
```
tid	idf	total_match_count	node_count	edge_count	max_node_cov	min_rmsd	nres	plddt	matching_residues	db_key	query_residues
data/serine_peptidases/4cha.pdb	0.6138	8	3	6	3	0.0000	477	13.5404	B57,B102,C195:0.0000;F57,F102,G195:0.0874	4	B57,B102,C195
data/serine_peptidases/1pq5.pdb	0.4869	4	3	4	3	0.2609	224	5.1340	A56,A99,A195:0.2609	3	B57,B102,C195
data/serine_peptidases/1ju3.pdb	0.0617	2	2	2	2	0.7792	570	19.4881	_,A223,A234:0.7792	1	B57,B102,C195
data/serine_peptidases/1l7a.pdb	0.0584	2	2	2	2	0.7883	636	11.7037	_,A146,A127:0.7883;_,B146,B127:0.8078	2	B57,B102,C195
data/serine_peptidases/1azw.pdb	0.1856	2	2	2	2	0.9234	626	34.2399	A179,_,B176:0.9234	0	B57,B102,C195
```
- `tid`: Identifier of the target protein structure
- `idf`: Inverse document frequency score with length penalty; Higher score indicates more matches within smaller structures
- `total_match_count`: Total number of matches
- `node_count`: Number of nodes in the structure
- `edge_count`: Number of edges in the structure
- `max_node_cov`: Maximum node coverage
- `min_rmsd`: Minimum root mean square deviation
- `nres`: Number of residues
- `plddt`: Predicted local distance difference test score
- `matching_residues`: Residue indices in the match (comma-separated, _ for no match, semicolon-separated for multiple matches with RMSD)
- `key`: Numeric identifier of the structure
- `query_residues`: Residue indices in the query (comma-separated)

### Display Options
- `--per-structure`: Outputs results per structure.
- `--per-match`: Outputs results per match.
- `--sort-by`: Sorts results by given columns (comma-separated). Defaults:
  `match_score:desc,rmsd:asc` per match, where `match_score` is `idf` × (matched fraction of the
  query)² × TM-score, and `structure_score:desc,min_rmsd:asc` per structure, where
  `structure_score` is matched residues² × √`idf` / (1 + RMSD). Both are also output columns.
  With `--skip-match` there is no match to score, so results stay in IDF order. For a family-level
  search whose filters already fix the coverage, `--sort-by max_node_count,min_rmsd` ranks by
  geometry alone.
- `--format-output`: Custom output format using column names.
- `--top <N>`: Outputs top N results.
- `--header`: Outputs header for the result.

## Contributions

<a href="https://github.com/steineggerlab/folddisco/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=steineggerlab/folddisco" />
</a>
