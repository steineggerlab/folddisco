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