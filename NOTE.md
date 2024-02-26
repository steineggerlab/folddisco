# Development note

## TODOs 
- [ ] IMPORTANT:  working examples
- [ ] IMPORTANT:  querying
- [ ] IMPORTANT:  confirm if sin & cos representation is working 

  - [ ] lib.rs
    - [x] expose necessary functions with prelude
  - [ ] tests
- [ ] Push only working code to the repository
  - [x] pass cargo check - 2024-02-20 18:30:00
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

- [ ] After finishing the above tasks, start working on the followings
    - [ ] Add feature for handling multiple hash types to the project

- NEEDED FEATURES
  - [ ] IndexTable with allocation (Priority: Middle)
  - [ ] Multiple hash types (Priority: High)
  - [ ] Working solution for this project (Priority: High)

- ADDITIONAL TODOS
  - [ ] Write rustdoc

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