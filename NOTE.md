# Development note

## TODOs 
- [ ] IMPORTANT: Clean up unnecessary files in these directories
  - [x] data: Keep minimal dataset for testing. REMOVE INTERMEDIATE FILES
  - [x] example
  - [ ] src
    - [ ] cli
    - [ ] controller
      - [ ] rename controller object as folddisco
    - [x] geometry
      - [x] core: changed to use enum for hash variants
      - [x] pdb_motif
      - [x] default
      - [x] changed angle representation to use sin and cos
      - [ ] trrosetta/trrosetta_subfamily
      - [x] first keep pdb_motif & trrosetta / remove others
      - [x] confirmed
      - [x] moved hash-related utils to geometry::util
    - [ ] index
    - [ ] structure
    - [ ] utils
  - [ ] lib.rs
  - [ ] tests
- [ ] Push only working code to the repository
- [ ] Add tests to the code
  - [ ] unit test (object level)
    - [x] geometry: pdb_motif, default
  - [ ] integration test (module level)
  - [ ] test & e (project level)

- [ ] After finishing the above tasks, start working on the followings
    - [ ] Add feature for handling multiple hash types to the project

- NEEDED FEATURES
  - [ ] IndexTable with allocation (Priority: Middle)
  - [ ] Multiple hash types (Priority: High)
  - [ ] Working solution for this project (Priority: High)

- ADDITIONAL TODOS
  - [ ] Write rustdoc



<!-- https://stats.stackexchange.com/questions/218407/encoding-angle-data-for-neural-network -->