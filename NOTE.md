# Development note

## TODOs 
- [ ] IMPORTANT: Clean up unnecessary files in these directories
  - [x] data: Keep minimal dataset for testing. REMOVE INTERMEDIATE FILES
  - [x] example
  - [ ] src
    - [ ] cli
    - [ ] controller
      - [ ] Integration 
        - [x] Structure
        - [x] Geometry
        - [ ] Index
      - [x] rename controller object as folddisco
    - [x] geometry
      - [x] core: changed to use enum for hash variants
      - [x] pdb_motif
      - [x] default
      - [x] changed angle representation to use sin and cos
        - [ ] confirm if sin & cos representation is working
      - [x] first keep pdb_motif & trrosetta / remove others
      - [x] confirmed
      - [x] moved hash-related utils to geometry::util
    - [ ] index
      - [ ] organize 
    - [ ] structure
    - [ ] utils
  - [ ] lib.rs
    - [x] moved hashablesync trait to top level
  - [ ] tests
- [ ] Push only working code to the repository
- [ ] Add tests to the code
  - [ ] unit test (object level)
    - [x] geometry: pdb_motif, default
    - [ ] controller
      - [x] hashing with pdb_motif
      - [x] hashing with default
      - [ ] indexing
      - [ ] querying
  - [ ] integration test (module level)
      - [ ] cli
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