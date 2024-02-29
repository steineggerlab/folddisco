# Development note

## TODOs 240229
RANKING & SCORING
- [x] DONE: Lookup original PDB file and retrieve the original structure. (Quite slow)
STRUCTURE
- [x] DONE: Last residue is not included in the motif
GEOMETRY
- [x] DONE: Modify PPF to use local CA as a reference point
- [x] DONE: SOLVED: PDBSinCos - uniqueness is not guaranteed
- [x] DONE: default hash type uses 64 bit integer and the hash seems to be too large. 
  - [x] Reduced default to use 32bit

## TODOs 
BIG THINGS
- [ ] TODO: Merge branch to main
- [ ] Reduce memory usage
- [ ] TODO: REMOVE MEMORY MONITORING PART AFTER THE MEMORY USAGE IS REDUCED ENOUGH
- [ ] IMPORTANT: confirm if sin & cos representation is working 

QUERYING
- [ ] TODO: allowing different amino acid pairs.
- [ ] TODO: FEATURE: multiple queries

GEOMETRY
- [ ] TODO: Add 3Di hash
- [ ] TODO: IMPORTANT: make binning and querying parameters configurable

INDEX
- [ ] IMPORTANT: Reduce memory usage with delta encoding (Make this as an option). DELTA ENCODING!!!
- [ ] IMPORTANT: Concat multiple index tables

BENCHMARK
- [ ] TODO: Benchmarking -- IMPORTANT: build CLI for this
- [ ] Benchmarking -- set benchmarking dataset based on PDB's approach
- [ ] TODO: Gather data

RANKING & SCORING
- [ ] TODO: Score output (RMSD maybe?)

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

<!-- https://stats.stackexchange.com/questions/218407/encoding-angle-data-for-neural-network -->

> 2024-02-16 19:27:55 Integration of geometry is done
```sh
    Finished test [optimized + debuginfo] target(s) in 0.63s
     Running tests/controller_test.rs (target/debug/deps/controller_test-28fe139faa771367)

running 1 test
0 | "data/homeobox/1akha-.pdb" | 2256 hashes
1 | "data/homeobox/1b72a-.pdb" | 4422 hashes
2 | "data/homeobox/1b72b-.pdb" | 5112 hashes
3 | "data/homeobox/1ba5--.pdb" | 2652 hashes
test test_fold_disco ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.01s
```

> 2024-02-20 18:30:23 Organized 4 hash types
> - PDBMotif (angle binned with degree)
> - PDBMotifSinCos (angle binned with sin and cos)
> - FoldDiscoDefault (trrosetta feature + aa pair)
> - TrRosetta (trrosetta feature)

> 2024-02-23 17:15:47

testing commands
```sh
# Indexing
./target/release/motifsearch index -d data/serine_peptidases_filtered -H default -i data/index/serine_peptidases_filtered -t 4 -v
# Querying
./target/release/motifsearch query -i data/index/serine_peptidases_filtered data/serine_peptidases_filtered/1aq2.pdb A250,A232,A269
```

## Built inverted index with sorting
> 2024-02-26 23:29:13
```sh
(base) hyunbin@ramda:/fast/hyunbin/motif/swissprot_benchmark$ \time -v ~/Projects/06_Motifsearch/motifsearch/target/release/motifsearch index2 -d ./swissprot_v4_raw/ -i ./folddisco_swissprot_v4/swissprot_v4_pdb -t 64 -v -H pdb

░█▀▀░█▀█░█░░░█▀▄░█▀▄░▀█▀░█▀▀░█▀▀░█▀█
░█▀▀░█░█░█░░░█░█░█░█░░█░░▀▀█░█░░░█░█
░▀░░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀

This is build_index2
[INFO] Building index table...
Time elapsed in "fold_disco.collect_hash_pairs()": 1861.898393297s # Sorting is included here
Time elapsed in "fold_disco.fill_numeric_id_vec()": 12.445695ms
Time elapsed in "fold_disco.set_index_table()": 743.989634429s     # This should be skipped
Time elapsed in "fold_disco.index_builder.convert_sorted_pairs_to_offset_and_values(fold_disco.hash_id_pairs)": 2970.673529343s # Current: Single threaded. NEED TO BE PARALLELIZED
Time elapsed in "save_offset_map(&offset_path,\n        &offset_table).expect(&log_msg(FAIL, \"Failed to save offset table\"))": 226.487304ms
Time elapsed in "write_usize_vector(&value_path,\n        &value_vec).expect(&log_msg(FAIL, \"Failed to save values\"))": 470.925270564s
Time elapsed in "save_lookup_to_file(&lookup_path, &fold_disco.path_vec,\n    &fold_disco.numeric_id_vec, None)": 185.856063ms
PDBMotifSinCos
[DONE] Done.
	Command being timed: "/home/hyunbin/Projects/06_Motifsearch/motifsearch/target/release/motifsearch index2 -d ./swissprot_v4_raw/ -i ./folddisco_swissprot_v4/swissprot_v4_pdb -t 64 -v -H pdb"
	User time (seconds): 33875.19
	System time (seconds): 2044.77
	Percent of CPU this job got: 593%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:40:57
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 1682107388
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 12
	Minor (reclaiming a frame) page faults: 150143064
	Voluntary context switches: 22725972
	Involuntary context switches: 224439
	Swaps: 0
	File system inputs: 22735136
	File system outputs: 559279392
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

```sh
[INFO] 
INITIAL - 0.031378746MB
SETUP - 0.03885269MB
collected: 3346869
memory_usage: 76MB
[INFO] fold_disco.collect_hash_pairs: 2.150621215s
AFTER COLLECTION - 76.672455MB
[INFO] fold_disco.sort_hash_pairs: 79.670447ms
AFTER SORTING - 76.68631MB
[INFO] fold_disco.fill_numeric_id_vec: 12.194µs
FILL NUM ID - 76.67482MB
[INFO] fold_disco.index_builder.convert_sorted_pairs_to_offset_and_values: 754.514141ms
OFFSET CONVERSION - 164.07326MB
[INFO] save_offset_map: 93.874368ms
SAVE OFFSET - 32.071117MB
[INFO] write_usize_vector: 10.31571ms
SAVE VALUE- 0.07117081MB
[INFO] save_lookup_to_file: 2.254281ms
SAVE LOOKUP - 0.07122421MB
FoldDiscoDefault
FINAL - 0.07127762MB
[DONE] Done.
```

```sh
/usr/local/bin/gtime -v ./target/release/motifsearch index -d analysis/raw_ecoli -i analysis/raw_ecoli_index -t 6 -v -H default
[INFO] ECOLI in Mac
INITIAL - 0.0429821MB
SETUP - 0.42183304MB
collected: 500747241
memory_usage: 11461MB
[INFO] fold_disco.collect_hash_pairs: 332.514011472s
AFTER COLLECTION - 11461.655MB
[INFO] fold_disco.sort_hash_pairs: 65.661235242s
AFTER SORTING - 11461.673MB
[INFO] fold_disco.fill_numeric_id_vec: 363.145µs
FILL NUM ID - 11461.727MB
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 45.420951796s
OFFSET CONVERSION - 4486.937MB
[INFO] save_offset_vec: 827.810023ms
SAVE OFFSET - 3820.9304MB
[INFO] write_usize_vector: 7.110334036s
SAVE VALUE- 0.53232193MB
[INFO] save_lookup_to_file: 13.876848ms
SAVE LOOKUP - 0.5323677MB
FoldDiscoDefault
FINAL - 0.5324135MB
[DONE] Done.
```

```sh
hbk@Hyunbinui-MacBookPro motifsearch % /usr/local/bin/gtime -v ./target/release/motifsearch index -d analysis/raw_ecoli -i analysis/raw_ecoli_pdb -t 6 -v -H pdb

░█▀▀░█▀█░█░░░█▀▄░█▀▄░▀█▀░█▀▀░█▀▀░█▀█
░█▀▀░█░█░█░░░█░█░█░█░░█░░▀▀█░█░░░█░█
░▀░░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀

[INFO] 
INITIAL - 0.04698944MB
SETUP - 0.4218254MB
collected: 92182758
memory_usage: 2109MB
[INFO] fold_disco.collect_hash_pairs: 75.935482089s
AFTER COLLECTION - 2110.338MB
[INFO] fold_disco.sort_hash_pairs: 3.475549827s
AFTER SORTING - 2110.385MB
[INFO] fold_disco.fill_numeric_id_vec: 32.214µs
FILL NUM ID - 2110.4282MB
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 1.262228495s
OFFSET CONVERSION - 712.70465MB
[INFO] save_offset_vec: 11.254174ms
SAVE OFFSET - 703.8309MB
[INFO] write_usize_vector: 468.940919ms
SAVE VALUE- 0.5323067MB
[INFO] save_lookup_to_file: 5.146061ms
SAVE LOOKUP - 0.53234863MB
PDBMotifSinCos
FINAL - 0.5323906MB
[DONE] Done.
        Command being timed: "./target/release/motifsearch index -d analysis/raw_ecoli -i analysis/raw_ecoli_pdb -t 6 -v -H pdb"
        User time (seconds): 380.70
        System time (seconds): 15.39
        Percent of CPU this job got: 487%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:21.24
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 3004036
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 3629307
        Voluntary context switches: 9553
        Involuntary context switches: 854724
        Swaps: 0
        File system inputs: 0
        File system outputs: 0
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

analysis/AF-P00776-F1-model_v4.pdb -q 149,171,253
trrosetta with aa pair
"analysis/raw_ecoli/AF-P39099-F1-model_v4.pdb"
> Changed to report non-singleton query results
analysis/raw_ecoli/AF-P39099-F1-model_v4.pdb: 3
analysis/raw_ecoli/AF-P0C0V0-F1-model_v4.pdb: 2
analysis/raw_ecoli/AF-P76176-F1-model_v4.pdb: 2
> All three are serine proteases

analysis/AF-P00776-F1-model_v4.pdb -q 149,171,253
pdb with aa pair
"analysis/raw_ecoli/AF-P0C0V0-F1-model_v4.pdb"
"analysis/raw_ecoli/AF-P09377-F1-model_v4.pdb"

analysis/AF-P00776-F1-model_v4.pdb -q 149,171,253
pdb with aa pair (binning just angle)
"analysis/raw_ecoli/AF-P09377-F1-model_v4.pdb"

> 2024-02-28 16:18:03
Confirmed zinc finger motif works

1ed8
https://pubs.acs.org/doi/10.1021/bi050155p