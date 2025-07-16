# Folddisco

Folddisco is a bioinformatics tool for indexing and searching discontinuous motifs in protein structures. 
It is designed to handle large-scale protein databases with unmatched speed and efficiency, 
enabling the detection of structural motifs across thousands of proteomes or millions of structures.

## Publications
[Kim H, Kim RS, Mirdita M, Steinegger M. Structural motif search across the protein-universe with Folddisco. bioRxiv, doi: 10.1101/2025.07.06.663357  (2025)](https://www.biorxiv.org/content/10.1101/2025.07.06.663357v1)

## Features
- Reduced index size, which enables large databases like AlphaFold to fit on a single disk
- Side-chain orientation-capturing feature and frequency-based scoring for higher precision
- Multi-threaded processing for fast indexing and querying

## Installation

### Default Installation
```bash
cargo install --features foldcomp --path .
```

### Build from Source
```bash
cargo build --release --features foldcomp
# Binary is located at target/release/folddisco
```

## Commands

### Indexing

#### Examples
```bash
# Default indexing for a small dataset
# h_sapiens directory or foldcomp database is indexed with default parameters
folddisco index -p h_sapiens -i index/h_sapiens -t 12

# Indexing big protein dataset
folddisco index -p swissprot -i index/swissprot -t 64 -m big -v

# Indexing with custom hash type and parameters
folddisco index -p h_sapiens -i index/h_sapiens -t 12 --type default -d 16 -a 4 # Default
folddisco index -p h_sapiens -i index/h_sapiens -t 12 --type pdb -d 8 -a 3 # PDB
```

#### Default Usage
```bash
folddisco index -p <PDB_DIR|FOLDCOMP_DB> -i <INDEX_PATH> -t <THREADS>
```

#### For Large Databases
```bash
folddisco index -p <PDB_DIR|FOLDCOMP_DB> -i <INDEX_PATH> -t <THREADS> -m big
```
- **Mode `big`:** Generates an 8GB fixed-size offset file suitable for datasets with more than 65,536 structures.

#### Custom Binning and Features
```bash
folddisco index -p <PDB_DIR|FOLDCOMP_DB> -i <INDEX_PATH> -t <THREADS> -d <DISTANCE_BINS> -a <ANGLE_BINS> -y <FEATURE_TYPE>
```

#### Example: Indexing the Human Proteome
```bash
folddisco index -p h_sapiens -i index/h_sapiens_folddisco -t 12
```

#### Pre-built Indices
Download pre-built index files:
- [Human proteome](https://foldcomp.steineggerlab.workers.dev/h_sapiens_folddisco.tar.gz)
- [E. coli proteome](https://foldcomp.steineggerlab.workers.dev/e_coli_folddisco.tar.gz)

### Querying

#### Examples
```bash
# Search with default settings. This will print out matching motifs with sorting by RMSD.
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
# Here's enolase query with 3 substitutions; Allow His at 164, Asp & Asn at 247, and His at 297.
folddisco query -p query/2MNR.pdb -q 164:H,195,221,247:ND,297:H -i index/e_coli_folddisco -d 0.5 -a 5 --top 10 --header --per-structure
# Range can be given with dash. This will query first 10 residues and 11th residue with subsitution to any amino acid.
folddisco query -p query/4CHA.pdb -q 1-10,11:X -i index/h_sapiens_folddisco -t 6 --serial-index

# Advanced query with filtering and sorting
## Based on connected node and rmsd
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --connected-node 0.75 --rmsd 1.0

## Coverage based filtering & top N filtering without residue matching
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --covered-node 3 --top 1000 --per-structure --skip-match

# Print top 100 structures with sorting by score
folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i index/h_sapiens_folddisco -t 6 --top 100 --per-structure --sort-by-score
folddisco query -q query/zinc_finger.txt -i index/h_sapiens_folddisco -t 6 --covered-node 4 --top 100 --sort-by-score --per-structure --skip-match
```


#### Default Usage
```bash
folddisco query -i <INDEX> -p <QUERY_PDB> -q <QUERY_RESIDUES> --skip-match -t <THREADS>
```
- `--skip-match`: Skips residue matching and RMSD calculation.
- `-v`: Verbose output.

#### Whole Structure as Query
```bash
folddisco query -i <INDEX> -p <QUERY_PDB> --skip-match -t <THREADS>
```

#### Using a Query File
```bash
folddisco query -i <INDEX> -q <QUERY_FILE> --skip-match -t <THREADS>
```

#### Distance and Angle Thresholds
```bash
folddisco query -i <INDEX> -p <QUERY_PDB> -q <QUERY_RESIDUES> -d <DISTANCE_THRESHOLD> -a <ANGLE_THRESHOLD> --skip-match -t <THREADS>
```

## Output
### Match Result
Default output which prints out one matching motif per line
```
id	node_count	avg_idf	rmsd	matching_residues	query_residues
AF-P00957-F1-model_v4.pdb	3	48.7694	0.2861	_,A666,A564,A568	F207,F212,F225,F229
AF-P0A6K3-F1-model_v4.pdb	3	58.2650	0.4315	A91,_,A133,A137	F207,F212,F225,F229
AF-P26649-F1-model_v4.pdb	2	36.0934	0.2204	A53,_,A22,_	F207,F212,F225,F229
AF-P05020-F1-model_v4.pdb	2	50.9269	0.3112	_,_,A17,A19	F207,F212,F225,F229
AF-P55798-F1-model_v4.pdb	2	62.1218	0.3725	_,_,A132,A14	F207,F212,F225,F229
```
- `id`: Identifier of the protein structure
- `node_count`: Number of nodes in the match
- `idf_score`: Inverse document frequency score of matched structure
- `rmsd`: Root mean square deviation
- `matching_residues`: Residue indices in the match (comma-separated, _ for no match)
- `query_residues`: Residue indices in the query (comma-separated)

### Structure Result
Output with one structure per line (`--per-structure`)
```
id	idf_score	total_match_count	node_count	edge_count	max_node_cov	min_rmsd	nres	plddt	matching_residues	query_residues
AF-P55798-F1-model_v4.pdb	62.1218	4	3	4	2	0.3725	218	95.4576	,,A132,A14:0.3725;,A39,A18,:0.6083	F207,F212,F225,F229
AF-P0A6K3-F1-model_v4.pdb	58.2650	4	3	4	3	0.4315	169	97.1329	A91,_,A133,A137:0.4315	F207,F212,F225,F229
AF-P05020-F1-model_v4.pdb	50.9269	4	3	4	3	0.4391	348	97.0974	,,A17,A19:0.3112;_,A222,A178,A203:0.4391	F207,F212,F225,F229
AF-P00957-F1-model_v4.pdb	48.7694	4	3	4	3	0.2861	876	90.7232	_,A666,A564,A568:0.2861	F207,F212,F225,F229
AF-P26649-F1-model_v4.pdb	36.0934	2	2	2	2	0.2204	66	75.4210	A53,,A22,:0.2204	F207,F212,F225,F229
```
- `id`: Identifier of the protein structure
- `idf_score`: Inverse document frequency score with length penalty; Higher score indicates more matches within smaller structures
- `total_match_count`: Total number of matches
- `node_count`: Number of nodes in the structure
- `edge_count`: Number of edges in the structure
- `max_node_cov`: Maximum node coverage
- `min_rmsd`: Minimum root mean square deviation
- `nres`: Number of residues
- `plddt`: Predicted local distance difference test score
- `matching_residues`: Residue indices in the match (comma-separated, _ for no match, semicolon-separated for multiple matches with RMSD)
- `query_residues`: Residue indices in the query (comma-separated)

### Display Options
- `--per-structure`: Outputs results per structure.
- `--per-match`: Outputs results per match.
- `--sort-by-score`: Sorts by score.
- `--sort-by-rmsd`: Sorts by RMSD.
- `--top <N>`: Outputs top N results.
- `--header`: Outputs header for the result.

## Example Index List
- **Human proteome:** `index/h_sapiens_folddisco` (23K structures, [Download](https://foldcomp.steineggerlab.workers.dev/h_sapiens_folddisco.tar.gz))
- **E. coli proteome:** `index/e_coli_folddisco` (4K structures, [Download](https://foldcomp.steineggerlab.workers.dev/e_coli_folddisco.tar.gz))

## Example Query List
- `query/`
  - `1G2F.pdb`: Zinc finger protein
  - `4CHA.pdb`: Serine protease
  - `1LAP.pdb`: Aminopeptidase
  - `zinc_finger.txt`: 1G2F.pdb F207,F212,F225,F229
  - `serine_protease.txt`: 4CHA.pdb B57,B102,C195
  - `aminopeptidase.txt`: 1LAP.pdb 250,255,273,332,334
  - `knottin.txt`: 2N6N.pdb 3,10,15,16,21,23,28,30
  - `enolase.txt`: 2MNR.pdb 164:H,195,221,247:ND,297:H

## Contributions

<a href="https://github.com/steineggerlab/motifsearch/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=steineggerlab/folddisco" />
</a>
