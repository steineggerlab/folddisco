# Development note

## TODOs 240409

DEV
- [x] DONE:Print original query
- [ ] Print node residue, 
- [ ] CLI::query_pdb: Output option
- [ ] CLI::query_pdb: Multiple queries by input file
- [ ] TODO: IMPORTANT: Restore retrieving functionality

IMPORTANT: BENCHMARK 
- [ ] Setup module & script
    - [x] DONE: benchmark module
    - [ ] benchmark script
- [x] DONE: Use hashset not to count duplicates
- [ ] TODO: check SCOP database
- [ ] Compare with pyscomotif
  - [ ] TODO: IMPORTANT: Download and rerun pyscomotif
  - [ ] TODO: Pyscomotif in two options

QUERYING
- [ ] Allow different amino acid pairs
  - [ ] TODO: Policies: Any, Exact, Same property
- [ ] Collect test query info / commands in QUERY.md

INDEXING
- [x] DONE: Add an option to save indices with different schemes
  - [ ] 3. ID + position

- [ ] Apply graph based algorithms (MASTER, ...) after running folddisco


## TODOs 

BENCHMARK
- [ ] Build an index of Swissprot & PDB
- [ ] Check if the query from other lab works or not
- [ ] Read MASTER, PDB realtime motif, pyscomotif on how they benchmarked

DEV
- [ ] TODO: expose necessary functions with prelude
- [ ] TODO: Setup output format --> easy to parse (tsv or other)
- [ ] CLI: polish grid related parameters
- [ ] CLI::index: Delete unncessary parameters
- [ ] CLI: polish logging

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

\time -v ~/Projects/06_Motifsearch/motifsearch/target/release/motifsearch query -p query/4CHA.pdb -q B57,B102,C195 -i h_sapiens_db/3di/index -t 32 > ~/serine.3di.tsv

```


```sh
(base) hyunbin@hulk:/fast/hyunbin/motif$ ~/Projects/06_Motifsearch/motifsearch/target/release/motifsearch index -p model_benchmark/h_sapiens/pdb/ -i h_sapiens_db/default32/d8a4/index -t 128 -y default32 -d 8 -a 4 -v --id uniprot -m id

░█▀▀░█▀█░█░░░█▀▄░█▀▄░▀█▀░█▀▀░█▀▀░█▀█
░█▀▀░█░█░█░░░█░█░█░█░░█░░▀▀█░█░░░█░█
░▀░░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀

[INFO] Indexing model_benchmark/h_sapiens/pdb/ with 128 threads and 1 chunks
[INFO] Hash type: Default32bit
[INFO] Indexing all PDB files in one chunk
[INFO] Collecting ids of the structures
[INFO] fold_disco.collect_hash_pairs: 835.665160066s
[INFO] Total 921636864 hashes collected (Allocated 21099.59MB)
[INFO] fold_disco.sort_hash_pairs: 9.796442734s
[INFO] Hash sorted (Allocated 21099.703MB)
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 11.134125762s
[INFO] Offset & values acquired (Allocated 9358.994MB)
[INFO] save_offset_vec: 20.450865102s
[INFO] write_usize_vector_in_bits: 11.350247254s
[INFO] save_lookup_to_file: 39.04581ms
[DONE] Indexing done for chunk 0 - h_sapiens_db/default32/d8a4/index
[DONE] Done.

```

> 2024-04-09 21:24:11
```sh
analysis/h_sapiens/d16a4/index_id       data/serine_pyscomotif.tsv      data/serine_answer.tsv  20504   198     124     PDBTrRosetta    16      4       108     20290   90      16      0.5455  0.8710  0.9948  0.6708
analysis/h_sapiens/d16a4/index_id       data/serine_folddisco.tsv       data/serine_answer.tsv  20504   105     124     PDBTrRosetta    16      4       101     20376   4       23      0.9619  0.8145  0.9987  0.8821
analysis/h_sapiens/d16a4/index_id       data/zinc_pyscomotif.tsv        data/zinc_answer.tsv    20504   304     1817    PDBTrRosetta    16      4       277     18669   27      1540    0.9112  0.1524  0.9236  0.2612
analysis/h_sapiens/d16a4/index_id       data/zinc_folddisco.tsv data/zinc_answer.tsv    20504   766     1817    PDBTrRosetta    16      4       753     18683   13      1064    0.9830  0.4144  0.9475  0.5830
```