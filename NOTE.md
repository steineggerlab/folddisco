# Development note

## TODOs 240327

IMPORTANT: BENCHMARK 
- [ ] TODO: scoring metric that takes into account the length of the structure & IDFs

DEV
- [ ] TODO: Split and extract ranking module

CLI::Query
- [x] DONE: New cutoffs introduced.
  - [ ] matched count cutoff (both absolute and relative: if integer, it's absolute. If float, it's relative)
  - [ ] score cutoff (float, default: 0)
  - [ ] num_residue cutoff <- MAKE THIS AS QUERY OPTION
  - [ ] TODO: plddt_cutoff <- DONE: THIS SHOULD BE SAVED IN LOOKUP TOO
  - [ ] Integrate and test these cutoffs
  - [ ] TODO: IMPORTANT: Add option to do wider search (dist_bin_threshold, angle_bin_threshold)
- [x] DONE: Set default chunk size to max (65535)

GEOMETRY
- [ ] TODO: Compare half match is needed or not?
- [ ] TODO: Add 3Di hash 
  - [ ] TODO: Fill in and integrate with the rest of the code

IMPORTANT: BENCHMARK 
- [ ] Setup module & script
- [ ] Build an index of PDB database (Running)
  - [ ] TODO: Rebuild one with nbin_dist = 16, nbin_angle = 3
- [ ] Check if the query from other lab works or not
- [ ] Read MASTER, PDB realtime motif, pyscomotif on how they benchmarked
- [ ] TODO: check SCOP database
- [ ] Compare with pyscomotif
  - [ ] TODO: IMPORTANT: Download and rerun pyscomotif

QUERYING
- [ ] TODO: measure time for querying with retrieval of matched positions
- [x] DONE: IDF seems to be working now.
  - [ ] TODO: Need to be tested
- [ ] TODO: IDEA: interactive mode that saves all the offsets in RAM
- [ ] TODO: FEATURE: multiple queries

INDEX
- [ ] IN_PROGRESS: Build Swissprot index (with max chunk size)
- [ ] NOTE: Jointly indexing with different hash types increases accuracy even without scoring
    - [ ] Current combination: PDBMotifSinCos, TrRosetta
    - [ ] As this doubles the size of the index, let's try with binning with smaller bins and modifying the querying schemes

DEV
- [ ] TODO: Polish logging
- [ ] TODO: Print original query
- [ ] TODO: expose necessary functions with prelude
- [ ] TODO: Setup output format --> easy to parse (tsv or other)

## TODOs 
BIG THINGS
- [ ] IMPORTANT: confirm if sin & cos representation is working 

QUERYING
- [ ] TODO: allowing different amino acid pairs.
- [ ] TODO: FEATURE: multiple queries
- [ ] TODO: measure time for querying with retrieval of matched positions

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


NOTE: Our algorithm is making a record-level inverted index, not a word-level inverted index.
TODO: Justify why we are using record-level inverted index and show the use case.
One of applications would be prefiltering for the next step of the analysis.


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
./target/release/motifsearch index -p data/serine_peptidases_filtered -H default -i data/index/serine_peptidases_filtered -t 4 -v

# Querying
./target/release/motifsearch query -i data/index/serine_peptidases_filtered data/serine_peptidases_filtered/1aq2.pdb A250,A232,A269

./target/release/motifsearch query -i <index> -p <pdb> data/serine_peptidases_filtered/4cha.pdb -q B57,B102,C195

# H. sapiens, serine peptidases
./target/release/motifsearch index -p analysis/h_sapiens_pdb -i analysis/h_sapiens_db/d8a3/index -t 8 -v -d 8 -a 3 -H pdb -c 65535
./target/release/motifsearch query -i analysis/h_sapiens_db/d8a3/index  -p data/serine_peptidases_filtered/4cha.pdb -q B57,B102,C195
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


```sh
# TODO: Pass cargo test


warning: unused variable: `mmap`
  --> src/cli/workflows/query_pdb.rs:61:18
   |
Error:           let (mmap, value_vec) = measure_time!(read_usize_vector(&value_path).expect("[ERROR] Failed to load value vector"));
   |                  ^^^^ help: if this is intentional, prefix it with an underscore: `_mmap`

warning: unused variable: `path_vec`
  --> src/cli/workflows/query_pdb.rs:65:18
running 33 tests
test cli::workflows::build_index::tests::test_build_index_pdb ... FAILED
test cli::workflows::build_index::tests::test_build_index_ppf ... FAILED
test cli::workflows::build_index::tests::test_build_index_trrosetta ... FAILED
test cli::workflows::query_pdb::tests::test_query_pdb_workflow ... FAILED
test controller::io::tests::test_offset_map_io ... ok
test controller::io::tests::test_usize_vector_io ... ok
test controller::query::tests::test_make_query ... ok
test controller::query::tests::test_parse_query_string ... ok
test controller::query::tests::test_parse_query_string_with_space ... ok
test controller::query::tests::test_parse_query_string_with_space_and_no_chain ... ok
test geometry::core::tests::test_hash_type ... ok
test geometry::default::tests::test_default_hash_works ... FAILED
test geometry::default_32bit::tests::test_default_hash_works ... ok
test geometry::pdb_motif::tests::test_hash_works ... ok
test geometry::pdb_motif_sincos::tests::test_geometrichash_works ... ok
test geometry::ppf::tests::test_hashvalue ... ok
test cli::workflows::build_index::tests::test_build_index ... FAILED
```


running 1 test
[INFO] Indexing data/serine_peptidases_filtered with 4 threads
[INFO] fold_disco.collect_hash_pairs: 2.106823529s
[INFO] Total 3228509 hashes collected (Allocated 73.949684MB)
[INFO] fold_disco.sort_hash_pairs: 74.154872ms
[INFO] Hash sorted (Allocated 73.957146MB)
[INFO] fold_disco.fill_numeric_id_vec: 1.128µs
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 133.454276ms
[INFO] Converted to offsets (Allocated 74.8718MB)
[INFO] save_offset_vec: 39.108513ms
[INFO] write_usize_vector: 8.053151ms
[INFO] save_lookup_to_file: 219.252µs
[DONE] Done.

running 1 test
[INFO] Indexing data/serine_peptidases_filtered with 4 threads
[INFO] fold_disco.collect_hash_pairs: 2.607902484s
[INFO] Total 1195107 hashes collected (Allocated 27.411839MB)
[INFO] fold_disco.sort_hash_pairs: 20.52749ms
[INFO] Hash sorted (Allocated 27.414333MB)
[INFO] fold_disco.fill_numeric_id_vec: 1.303µs
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 21.74806ms
[INFO] Converted to offsets (Allocated 14.493169MB)
[INFO] save_offset_vec: 4.258504ms
[INFO] write_usize_vector: 3.430338ms
[INFO] save_lookup_to_file: 219.809µs
[DONE] Done.

running 1 test
[INFO] Indexing data/serine_peptidases_filtered with 4 threads
[INFO] fold_disco.collect_hash_pairs: 597.892004ms
[INFO] Total 607023 hashes collected (Allocated 13.939382MB)
[INFO] fold_disco.sort_hash_pairs: 11.692878ms
[INFO] Hash sorted (Allocated 13.954213MB)
[INFO] fold_disco.fill_numeric_id_vec: 954ns
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 17.915253ms
[INFO] Converted to offsets (Allocated 10.009801MB)
[INFO] save_offset_vec: 7.152273ms
[INFO] write_usize_vector: 1.818618ms
[INFO] save_lookup_to_file: 286.759µs
[DONE] Done.
test cli::workflows::build_index::tests::test_build_index ... ok


```sh

░█▀▀░█▀█░█░░░█▀▄░█▀▄░▀█▀░█▀▀░█▀▀░█▀█
░█▀▀░█░█░█░░░█░█░█░█░░█░░▀▀█░█░░░█░█
░▀░░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀

[INFO] Indexing analysis/raw_ecoli with 6 threads and 1 chunks
[INFO] Indexing all PDB files in one chunk
[INFO] fold_disco.collect_hash_pairs: 75.025886562s
[INFO] Total 92442306 hashes collected (Allocated 2116.5344MB)
[INFO] fold_disco.sort_hash_pairs: 2.835071084s
[INFO] Hash sorted (Allocated 2116.5508MB)
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 1.654009301s
[INFO] Offset & values acquired (Allocated 714.9261MB)
[INFO] save_offset_vec: 17.190798ms
[INFO] write_usize_vector_in_bits: 350.628364ms
[INFO] save_lookup_to_file: 6.863039ms
[DONE] Indexing done for chunk 0 - analysis/raw_ecoli_pdb
[DONE] Done.
-rw-r--r--      1 hbk  staff   221K Mar  6 12:31 raw_ecoli_pdb.lookup
-rw-r--r--      1 hbk  staff   5.5M Mar  6 12:31 raw_ecoli_pdb.offset
-rw-r--r--      1 hbk  staff    14B Mar  6 12:31 raw_ecoli_pdb.type
-rw-r--r--      1 hbk  staff   176M Mar  6 12:31 raw_ecoli_pdb.value
hbk@Hyunbinui-MacBookPro motifsearch % ./target/release/motifsearch index -d analysis/h_sapiens -i analysis/h_sapiens_index -t 6 -v -H pdb -c 2400

░█▀▀░█▀█░█░░░█▀▄░█▀▄░▀█▀░█▀▀░█▀▀░█▀█
░█▀▀░█░█░█░░░█░█░█░█░░█░░▀▀█░█░░░█░█
░▀░░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀

[INFO] Indexing analysis/h_sapiens with 6 threads and 10 chunks
[INFO] Indexing chunk 0
[INFO] fold_disco.collect_hash_pairs: 174.94383012s
[INFO] Total 66828939 hashes collected (Allocated 1531.52MB)
[INFO] fold_disco.sort_hash_pairs: 1.703404673s
[INFO] Hash sorted (Allocated 1531.5333MB)
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 769.431081ms
[INFO] Offset & values acquired (Allocated 520.72797MB)
[INFO] save_offset_vec: 9.896456ms
[INFO] write_usize_vector_in_bits: 192.254408ms
[INFO] save_lookup_to_file: 4.285185ms
```


-rw-rw-r-- 1 hyunbin steineggerlab 1904485972 Mar  6 13:09 swissprot_default32_84.value
-rw-rw-r-- 1 hyunbin steineggerlab     358986 Mar  6 03:29 swissprot_default32_9.lookup
-rw-rw-r-- 1 hyunbin steineggerlab  696379300 Mar  6 03:29 swissprot_default32_9.offset
-rw-rw-r-- 1 hyunbin steineggerlab         12 Mar  6 03:29 swissprot_default32_9.type
-rw-rw-r-- 1 hyunbin steineggerlab 1837116642 Mar  6 03:29 swissprot_default32_9.value

```sh
Command being timed: "/home/hyunbin/Projects/06_Motifsearch/motifsearch/target/release/motifsearch index -d ./swissprot_benchmark/swissprot_v4_raw/ -i ./swissprot_benchmark/default/swissprot_default32 -v -t 100 -H default32 -c 5000"
User time (seconds): 4258543.22
System time (seconds): 13097.09
Percent of CPU this job got: 8487%
Elapsed (wall clock) time (h:mm:ss or m:ss): 13:58:50
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 46282292
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 24987
Minor (reclaiming a frame) page faults: 732009010
Voluntary context switches: 88365756
Involuntary context switches: 12901302
Swaps: 0
File system inputs: 73534872
File system outputs: 528036232
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
  ```
  
```sh
Latest log
hbk@arm64-apple-darwin20 motifsearch % ./target/release/motifsearch index -p analysis/h_sapiens_pdb -i analysis/h_sapiens_db/d16a3/index -t 8 -v -d 16 -a 3 -H pdb -c 65535         

░█▀▀░█▀█░█░░░█▀▄░█▀▄░▀█▀░█▀▀░█▀▀░█▀█
░█▀▀░█░█░█░░░█░█░█░█░░█░░▀▀█░█░░░█░█
░▀░░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀

[INFO] Indexing analysis/h_sapiens_pdb with 8 threads and 1 chunks
[INFO] Indexing all PDB files in one chunk
[INFO] fold_disco.collect_hash_pairs: 2037.426889583s
[INFO] Total 508214706 hashes collected (Allocated 11635.57MB)
[INFO] fold_disco.sort_hash_pairs: 16.490825333s
[INFO] Hash sorted (Allocated 11635.64MB)
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 9.149903667s
[INFO] Offset & values acquired (Allocated 3886.0242MB)
[INFO] save_offset_vec: 4.082541ms
[INFO] write_usize_vector_in_bits: 1.458560792s
[INFO] save_lookup_to_file: 5.792792ms
[DONE] Indexing done for chunk 0 - analysis/h_sapiens_db/d16a3/index
[DONE] Done.
```

### H. sapiens PPF / TrRosetta
```sh
(base) hyunbin@hulk:~/Projects/06_Motifsearch/motifsearch$ cd /fast/hyunbin/motif/
(base) hyunbin@hulk:/fast/hyunbin/motif$ ~/Projects/06_Motifsearch/motifsearch/target/release/motifsearch index -p model_benchmark/h_sapiens/data/ -i ./h_sapiens_db/h_sapiens_ppf -t 200 -v -H ppf -c 65535

░█▀▀░█▀█░█░░░█▀▄░█▀▄░▀█▀░█▀▀░█▀▀░█▀█
░█▀▀░█░█░█░░░█░█░█░█░░█░░▀▀█░█░░░█░█
░▀░░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀

[INFO] Indexing model_benchmark/h_sapiens/data/ with 200 threads and 1 chunks
[INFO] Indexing all PDB files in one chunk
[INFO] fold_disco.collect_hash_pairs: 3740.329732938s
[INFO] Total 717406066 hashes collected (Allocated 16425.266MB)
[INFO] fold_disco.sort_hash_pairs: 8.021176973s
[INFO] Hash sorted (Allocated 16425.34MB)
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 8.130028772s
[INFO] Offset & values acquired (Allocated 5481.367MB)
[INFO] save_offset_vec: 93.807158ms
[INFO] write_usize_vector_in_bits: 3.687907908s
[INFO] save_lookup_to_file: 45.950435ms
[DONE] Indexing done for chunk 0 - ./h_sapiens_db/h_sapiens_ppf
[DONE] Done.
(base) hyunbin@hulk:/fast/hyunbin/motif$ ~/Projects/06_Motifsearch/motifsearch/target/release/motifsearch index -p model_benchmark/h_sapiens/data/ -i ./h_sapiens_db
/h_sapiens_trrosetta -t 200 -v -H trrosetta -c 65535

░█▀▀░█▀█░█░░░█▀▄░█▀▄░▀█▀░█▀▀░█▀▀░█▀█
░█▀▀░█░█░█░░░█░█░█░█░░█░░▀▀█░█░░░█░█
░▀░░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀

[INFO] Indexing model_benchmark/h_sapiens/data/ with 200 threads and 1 chunks
[INFO] Indexing all PDB files in one chunk
[INFO] fold_disco.collect_hash_pairs: 2091.266188905s
[INFO] Total 548859639 hashes collected (Allocated 12567.43MB)
[INFO] fold_disco.sort_hash_pairs: 11.914579738s
[INFO] Hash sorted (Allocated 12567.897MB)
[INFO] convert_sorted_pairs_to_offset_and_values_vec: 7.529235546s
[INFO] Offset & values acquired (Allocated 4195.5864MB)
[INFO] save_offset_vec: 24.456832ms
[INFO] write_usize_vector_in_bits: 2.80209816s
[INFO] save_lookup_to_file: 30.912753ms
[DONE] Indexing done for chunk 0 - ./h_sapiens_db/h_sapiens_trrosetta
[DONE] Done.

```