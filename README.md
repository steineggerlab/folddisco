# Folddisco

Folddisco is a bioinformatics tool for indexing and searching
for discontinous motifs in protein structures.

TODO: Write novel features of folddisco

## Command
### Installation
```bash
# Default
cargo install --features foldcomp --path .
# Build from source
cargo build --release --features foldcomp # folddisco binary is located in target/release/
```

### Indexing
```bash
# Default for small databases
folddisco index -p <PDB_DIR|FOLDCOMP_DB> -i <INDEX_PATH> -t <THREADS> 
# For large databases (number of structures > 65536). This will generate a 8GB fixed-size offset file
folddisco index -p <PDB_DIR|FOLDCOMP_DB> -i <INDEX_PATH> -t <THREADS> -m big

# Setting the number of bins and feature
folddisco index -p <PDB_DIR|FOLDCOMP_DB> -i <INDEX_PATH> -t <THREADS> -d <DISTANCE_BINS> -a <ANGLE_BINS> -y <FEATURE_TYPE>
```

```bash
# Example
# Indexing human proteome with 12 threads
folddisco index -p h_sapiens -i index/h_sapiens_folddisco -t 12
```

You can also download pre-built index files 
- [Human proteome](https://foldcomp.steineggerlab.workers.dev/h_sapiens_folddisco.tar.gz)
- [E. coli proteome](https://foldcomp.steineggerlab.workers.dev/e_coli_folddisco.tar.gz)

### Querying
- [ ] TODO: --retrieve flag has been removed. Update the documentation --> --skip-match
- [ ] TODO: Add information on --per-structure; per-match. --sort-by-rmsd 
- [ ] TODO: Add column information for the output file with example
- [ ] TODO: Add example in example directory
```bash
# default
folddisco query -i <INDEX> -p <QUERY_PDB> -q <QUERY_RESIDUES> -r -t <THREADS>
# -r flag is for residue matching & rmsd calculation
# -v flag is for verbose output (measures step-by-step runtime)

# Using the whole structure as a query
folddisco query -i <INDEX> -p <QUERY_PDB> -r -t <THREADS>

# Using query.txt file
folddisco query -i <INDEX> -q <QUERY_FILE> -r -t <THREADS>

# Using distance and angle threshold
folddisco query -i <INDEX> -p <QUERY_PDB> -q <QUERY_RESIDUES> -d <DISTANCE_THRESHOLD> -a <ANGLE_THRESHOLD> -r -t <THREADS>
```

```bash
# Example
# Zinc finger query to human proteome
folddisco query -i index/h_sapiens_folddisco -p query/1G2F.pdb -q F207,F212,F225,F229 -r -d 0.5 -a 5 -t 12
folddisco query -i index/h_sapiens_folddisco -q query/zinc_finger.txt -r -d 0.5 -a 5 -t 12

# Serine protease query to human proteome
folddisco query -i index/h_sapiens_folddisco -p query/4CHA.pdb -q B57,B102,C195 -r -t 12
folddisco query -i index/h_sapiens_folddisco -q query/serine_protease.txt -r -t 12
```


## Index list
- `index/`
  - `h_sapiens_folddisco`: Human proteome, 23K structures
  - `e_coli_folddisco`: E. coli proteome, 4K structures

## Example query list
- `query/`
  - `1G2F.pdb`: Zinc finger protein
  - `4CHA.pdb`: Serine protease
  - `1LAP.pdb`: Aminopeptidase
  - `zinc_finger.txt`: 1G2F.pdb F207,F212,F225,F229
  - `serine_protease.txt`: 4CHA.pdb B57,B102,C195
  - `aminopeptidase.txt`: 1LAP.pdb 250,255,273,332,334
  - `knottin.txt`: 2N6N.pdb 3,10,15,16,21,23,28,30


## Contributor
<a href="https://github.com/steineggerlab/motifsearch/graphs/contributors">
  <img src="https://contributors-img.firebaseapp.com/image?repo=steineggerlab/folddisco" />
</a>
