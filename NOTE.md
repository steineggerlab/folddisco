# Development note

## TODOs 240404
QUERYING
- [x] DONE: Comprehensive filtering parameters: node coverage, edge coverage, exact match, total match, grid count, check all grid is nearby
- [ ] Allow different amino acid pairs
  - [ ] TODO: Policies: Any, Exact, Same property
- [ ] Collect test query info / commands in QUERY.md

DEV
- [ ] CLI: polish grid related parameters
- [ ] CLI::index: Delete unncessary parameters
- [x] DONE: CLI::query_pdb: filter how??
- [ ] CLI: polish logging
- [x] DONE: CLI::query_pdb: Extract functions and measure time
- [ ] CLI::query_pdb: Log the original query
- [ ] CLI::query_pdb: Output option
- [x] DONE: CLI::query_pdb: Print only the base name of the pdb file
- [x] DONE: CLI::index: Save only the base name of the pdb file (optional)
- [ ] CLI::query_pdb: Multiple queries by input file

INDEXING
- [x] DONE: Add an option to save indices with different schemes
  - [x] 1. ID only
  - [x] 2. ID + grid
  - [ ] 3. ID + position

GEOMETRY
- [ ] TODO: Add 3Di hash 


## TODOs 

IMPORTANT: BENCHMARK 
- [ ] Setup module & script
- [ ] Build an index of PDB database
  - [ ] TODO: Rebuild one with nbin_dist = 16, nbin_angle = 4
- [ ] Build an index of Swissprot
- [ ] Check if the query from other lab works or not
- [ ] Read MASTER, PDB realtime motif, pyscomotif on how they benchmarked
- [ ] TODO: check SCOP database
- [ ] Compare with pyscomotif
  - [ ] TODO: IMPORTANT: Download and rerun pyscomotif

DEV
- [ ] TODO: Split and extract ranking module
  - [ ] TODO: Fill in and integrate with the rest of the code
- [ ] TODO: Polish logging
- [ ] TODO: Print original query
- [ ] TODO: expose necessary functions with prelude
- [ ] TODO: Setup output format --> easy to parse (tsv or other)

QUERYING
- [ ] TODO: IDEA: interactive mode that saves all the offsets in RAM
- [ ] TODO: FEATURE: multiple queries


BIG THINGS
- [ ] IMPORTANT: confirm if sin & cos representation is working 

QUERYING
- [ ] TODO: allowing different amino acid pairs.
- [ ] TODO: FEATURE: multiple queries
- [ ] TODO: measure time for querying with retrieval of matched positions

INDEX
- [ ] IMPORTANT: Reduce memory usage with delta encoding (Make this as an option). DELTA ENCODING!!!
- [ ] IMPORTANT: Concat multiple index tables

BENCHMARK
- [ ] TODO: Benchmarking -- IMPORTANT: build CLI for this
- [ ] Benchmarking -- set benchmarking dataset based on PDB's approach
- [ ] TODO: Gather data

IDEAS
- [ ] IDEA: Sukhwan: How powerful is amino acid pair? revival of aa_pair
- [ ] Longer motifs?
- [ ] IDEA: MINOR: residue matching strategy 
- [ ] IDEA: Naming of the project
- [ ] IndexTable with allocation (Priority: Middle)

DEV
- [ ] TODO:Polish logging
- [ ] Print original query
- [ ] TODO: expose necessary functions with prelude
- [ ] working examples
- [ ] TODO: Write rustdoc
- NOTE: Push only working code to the repository

TEST
- [ ] Check all structs and methods are working within tests
  - [ ] structure
  - [ ] geometry
  - [ ] index
  - [ ] unit test (object level)
    - [ ] controller
      - [ ] querying
  - [ ] integration test (module level)
      - [ ] cli
        - [x] build_index
        - [ ] query
        - [ ] benchmark
  - [ ] project level test


NOTE: Our algorithm is making a record-level inverted index, not a word-level inverted index.
TODO: Justify why we are using record-level inverted index and show the use case.
One of applications would be prefiltering for the next step of the analysis.



1ed8
https://pubs.acs.org/doi/10.1021/bi050155p

<!-- https://stats.stackexchange.com/questions/218407/encoding-angle-data-for-neural-network -->

> 2024-02-20 18:30:23 Organized 4 hash types
> - PDBMotif (angle binned with degree)
> - PDBMotifSinCos (angle binned with sin and cos)
> - FoldDiscoDefault (trrosetta feature + aa pair)
> - TrRosetta (trrosetta feature)

> 2024-02-23 17:15:47

testing commands
```sh
# Indexing
./target/release/motifsearch index -p data/serine_peptidases_filtered -H default -i data/index/serine_peptidases_filtered -t 4 -v

# Querying
./target/release/motifsearch query -i data/index/serine_peptidases_filtered data/serine_peptidases_filtered/1aq2.pdb A250,A232,A269

./target/release/motifsearch query -i <index> -p <pdb> data/serine_peptidases_filtered/4cha.pdb -q B57,B102,C195

# H. sapiens, serine peptidases
./target/release/motifsearch index -p analysis/h_sapiens_pdb -i analysis/h_sapiens_db/d8a3/index -t 8 -v -d 8 -a 3 -H pdb -c 65535
./target/release/motifsearch query -i analysis/h_sapiens_db/d8a3/index  -p data/serine_peptidases_filtered/4cha.pdb -q B57,B102,C195
```
