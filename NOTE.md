# Development note

## TODOs 240522
QUERYING
- [x] DONE: IMPORTANT: RMSD calculation for matched positions
  - [x] DONE: NEED superposition for the matched positions
  - [x] Implemented QCP superposition
  - [x] DONE: Matching query and reference positions
  - [x] DONE: rayon based parallelization
- [x] DONE: For partial match, fill in with "_" in the output where the query does not match
- [x] DONE: SimpleHashMap that can be accessed with memmap2
- [x] DONE: Faster loading of offset using SimpleHashMap; 500ms - 2s --> 100-500 microseconds
- [ ] IMPORTANT: TODO: Rename features --> default goes to PDBTrRosetta
- [ ] IMPORTANT: TODO: Delete 64bit hashes; Only keep 32bit
  
DEV
- [ ] Write rustdoc

BENCHMARK
- [ ] TODO: Index should not be a mandatory parameter --> Fix this
- Benchmark plan
  - 1. Known motifs
    - Compare with pyscomotif, PDB
    - Serine peptidase, Zinc finger, Aminopeptidase, 
  - 2. Random motifs (SCOP)
    - SCOP or SCOP 40

SCORING
- [ ] MASTER

IMPORTANT: BENCHMARK 
  - [ ] TODO: Pyscomotif in two options

QUERYING
- [ ] Allow different amino acid pairs
  - [ ] TODO: Policies: Any, Exact, Same property
- [ ] Collect test query info / commands in QUERY.md
- [ ] Apply graph based algorithms (MASTER, ...) after running folddisco

INDEXING
- [x] DONE: Add an option to save indices with different schemes
  - [ ] 3. ID + position

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
analysis/h_sapiens/d16a4/index_id       temp.ser.tsv                    data/serine_answer.tsv  20504   643     124     PDBTrRosetta    16      4       114     19851   529     10      0.1773  0.9194  0.9737  0.2973
analysis/h_sapiens/d16a4/index_id       data/serine_folddisco.tsv       data/serine_answer.tsv  20504   105     124     PDBTrRosetta    16      4       101     20376   4       23      0.9619  0.8145  0.9987  0.8821
analysis/h_sapiens/d16a4/index_id       data/zinc_pyscomotif.tsv        data/zinc_answer.tsv    20504   304     1817    PDBTrRosetta    16      4       277     18669   27      1540    0.9112  0.1524  0.9236  0.2612
analysis/h_sapiens/d16a4/index_id       temp.zinc.tsv                   data/zinc_answer.tsv    20504   1409    1817    PDBTrRosetta    16      4       895     18182   514     922     0.6352  0.4926  0.9300  0.5549
analysis/h_sapiens/d16a4/index_id       data/zinc_folddisco.tsv         data/zinc_answer.tsv    20504   766     1817    PDBTrRosetta    16      4       753     18683   13      1064    0.9830  0.4144  0.9475  0.5830
```

# Things to keep in mind when developing
* No dependencies between huge modules
  * geometry ←→ structure ←→ index
* Make errors visible
  * Avoid unwrap, use Option/Result/expect/match
* Write tests
  * Unit tests at source files
  * Integration tests at `tests/`


analysis/e_coli/pdb/AF-P76176-F1-model_v4.pdb   76.3243 6       3       6       2       0       1       273     89.81777        A145,A223,A84   B57,B102,C195   analysis/e_coli/pdbtr/d16a4/index

ATOM    391  CA  HIS B  57       6.994   8.354  42.405  1.00  7.59           C  
ATOM    733  CA  ASP B 102       9.429   7.479  48.266  1.00  8.81           C  
ATOM   1392  CA  SER C 195       5.547   0.158  42.050  1.00  7.92           C  

ATOM   1084  CA  ASP A 145     -12.833   3.134  -7.780  1.00 97.66           C  
ATOM   1669  CA  SER A 223      -5.720  -2.218  -3.368  1.00 98.13           C  
ATOM    608  CA  HIS A  84     -13.958  -1.741  -4.223  1.00 98.13           C  

analysis/h_sapiens/pdb/AF-Q9BZJ3-F1-model_v4.pdb        68.08451        6       3       6       4       0       1       242     85.22195        A128,A231,A81   B57,B102,C195   analysis/h_sapiens/index

ATOM    612  CA  HIS A  81      -4.924   5.813  -9.485  1.00 96.40           C  
ATOM   1002  CA  ASP A 128      -0.499  10.073  -8.059  1.00 97.81           C  
ATOM   1796  CA  SER A 231      -0.792   0.658  -4.430  1.00 96.12           C  


---

> 2024-05-30 20:31:03
```sh
# Indexing E. coli
(base) hyunbin@super003:/fast/hyunbin/motif$ ~/Projects/06_Motifsearch/motifsearch/target/release/folddisco index -t 64 -p model_benchmark/e_coli/pdb/ -i e_coli_db/pdbtr/d16a4/index -d 16 -a 4 -v -m id -y pdbtr

░█▀▀░█▀█░█░░░█▀▄░█▀▄░▀█▀░█▀▀░█▀▀░█▀█
░█▀▀░█░█░█░░░█░█░█░█░░█░░▀▀█░█░░░█░█
░▀░░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀

[INFO] Indexing model_benchmark/e_coli/pdb/ with 64 threads and 1 chunks
[INFO] Hash type: PDBTrRosetta
[INFO] Indexing all PDB files in one chunk
[INFO] Collecting ids of the structures
[INFO] fold_disco.collect_hash_pairs: 46.725618018s
[INFO] Total 121589715 hashes collected (Allocated 2783.9915MB)
[INFO] fold_disco.sort_hash_pairs: 1.452377548s
[INFO] Hash sorted (Allocated 2784.2754MB)
[INFO] if let Some: 1.202668851s
[INFO] match index_mode {
    IndexMode::Id => {
        convert_sorted_hash_pairs_to_simplemap: 1.479611419s
[INFO] Offset & values acquired (Allocated 1788.4944MB)
[INFO] offset_map.dump_to_disk: 609.964463ms
[INFO] write_usize_vector_in_bits: 476.190861ms
[INFO] save_lookup_to_file: 4.392093ms
[DONE] Indexing done for chunk 0 - e_coli_db/pdbtr/d16a4/index

# Querying E. coli
(base) hyunbin@super003:/fast/hyunbin/motif$ ~/Projects/06_Motifsearch/motifsearch/target/release/folddisco query -p query/4CHA.pdb -q B57,B102,C195 -i e_coli_db/pdbtr/d16a4/index -d 0.5 -a 5 -v -r

░█▀▀░█▀█░█░░░█▀▄░█▀▄░▀█▀░█▀▀░█▀▀░█▀█
░█▀▀░█░█░█░░░█░█░█░█░░█░░▀▀█░█░░░█░█
░▀░░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀▀░▀▀▀░▀▀▀░▀▀▀

[INFO] Querying query/4CHA.pdb:B57,B102,C195 to e_coli_db/pdbtr/d16a4/index
[INFO] Found 1 index file(s). Querying with 1 threads
[INFO] SimpleHashMap::load_from_disk: 457.924µs
[INFO] load_lookup_from_file: 5.538099ms
[INFO] make_query_map: 2.055914ms
[INFO] read_u16_vector: 440.195µs
[INFO] count_query_idmode: 269.826µs
[INFO] queried_from_indices.par_sort_by: 9.82µs
id	idf_score	total_match_count	node_count	edge_count	exact_match_count	overflow_count	grid_count	nres	plddt	matching_residues	query_residues	query_file	index_path
model_benchmark/e_coli/pdb/AF-P76176-F1-model_v4.pdb	76.3243	6	3	6	2	0	1	273	89.8178	A84,A145,A223:0.1256	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P39099-F1-model_v4.pdb	49.3797	4	3	4	2	0	1	455	86.1556	A109,A139,A214:0.1160	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0C0V0-F1-model_v4.pdb	49.1437	4	3	4	2	0	1	474	85.4206	A131,A161,A236:0.1310	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P21517-F1-model_v4.pdb	29.7054	2	2	2	0	0	1	604	95.6333	A584,_,A598:0.3684	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0AEE3-F1-model_v4.pdb	28.0689	2	2	2	2	0	1	355	87.3617	A96,A126,_:0.1319	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P76222-F1-model_v4.pdb	27.5519	2	2	2	2	0	1	182	94.4896	_,A95,A106:0.6922	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0AE63-F1-model_v4.pdb	27.1905	2	2	2	0	0	1	76	87.1307	_,A70,A12:0.2013	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P52139-F1-model_v4.pdb	26.5164	2	2	2	0	0	1	152	90.6160	A29,A26,_:0.6851	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P11553-F1-model_v4.pdb	26.4169	2	2	2	2	0	1	472	95.6012	A389,_,A271:0.4888	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P23890-F1-model_v4.pdb	26.1822	2	2	2	2	0	1	512	86.9527	A344,_,A218:0.4947	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P46859-F1-model_v4.pdb	26.1099	2	2	2	0	0	1	175	93.3436	A43,A40,_:0.4780	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P30178-F1-model_v4.pdb	26.0206	2	2	2	0	0	1	361	97.4951	_,A46,A334:0.3938	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P17444-F1-model_v4.pdb	25.9443	2	2	2	2	0	1	556	96.2182	A417,_,A391:0.6174	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P39208-F1-model_v4.pdb	25.9185	2	2	2	0	0	1	187	92.1190	A38,A35,_:0.4767	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P52125-F1-model_v4.pdb	25.6114	2	2	2	0	0	1	208	88.2736	A91,A88,_:0.3894	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P32672-F1-model_v4.pdb	25.5918	2	2	2	2	0	1	359	91.9060	_,A51,A36:0.4375	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P78067-F1-model_v4.pdb	25.4825	2	2	2	0	0	1	435	94.7385	_,A263,A74:0.5848	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0AGG4-F1-model_v4.pdb	25.4485	2	2	2	0	0	1	139	95.4991	_,A31,A82:0.3645	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P60778-F1-model_v4.pdb	25.3307	2	2	2	2	0	1	393	92.5210	_,A138,A308:0.4686	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P77554-F1-model_v4.pdb	25.3213	2	2	2	0	0	1	460	90.9688	_,A125,A191:0.4789	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P11664-F1-model_v4.pdb	25.2348	2	2	2	0	0	1	237	95.6755	A85,A82,_:0.3820	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P55138-F1-model_v4.pdb	25.1273	2	2	2	0	0	1	492	94.2079	_,A209,A228:0.5283	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P23883-F1-model_v4.pdb	25.1097	2	2	2	0	0	1	495	97.6602	_,A120,A199:0.3781	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P76319-F1-model_v4.pdb	25.0606	2	2	2	0	0	1	159	88.7526	_,A142,A29:0.4372	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P39835-F1-model_v4.pdb	25.0179	2	2	2	2	0	1	438	87.6747	_,A261,A60:0.3139	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P37685-F1-model_v4.pdb	25.0123	2	2	2	0	0	1	512	96.8134	_,A123,A202:0.3721	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P21693-F1-model_v4.pdb	24.8954	2	2	2	2	0	1	457	90.9402	_,A42,A195:0.4543	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P11663-F1-model_v4.pdb	24.8846	2	2	2	0	0	1	169	93.5550	_,A66,A92:0.6719	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0AEX5-F1-model_v4.pdb	24.6624	2	2	2	0	0	1	289	93.2634	A45,A42,_:0.5656	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P16431-F1-model_v4.pdb	24.2629	2	2	2	2	0	1	569	92.1165	_,A481,A237:0.4493	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P09980-F1-model_v4.pdb	23.7786	2	2	2	2	0	1	673	88.4221	_,A109,A379:0.1914	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P36667-F1-model_v4.pdb	23.5976	2	2	2	0	0	1	264	91.6702	_,A195,A130:0.2351	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0A6C1-F1-model_v4.pdb	23.3767	2	2	2	0	0	1	285	97.1442	_,A70,A87:0.4835	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P77286-F1-model_v4.pdb	23.2839	2	2	2	0	0	1	466	91.8885	A177,A174,_:0.4359	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P05793-F1-model_v4.pdb	23.1331	2	2	2	0	0	1	491	95.7749	A111,A108,_:0.7473	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0AFE8-F1-model_v4.pdb	23.0292	2	2	2	0	0	1	509	93.7413	A490,A487,_:0.7770	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0A9T4-F1-model_v4.pdb	22.8171	2	2	2	0	0	1	346	96.6186	_,A127,A9:0.6353	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P77360-F1-model_v4.pdb	22.7593	2	2	2	0	0	1	353	96.3956	_,A140,A157:0.5120	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P29018-F1-model_v4.pdb	22.6129	2	2	2	0	0	1	588	89.6803	A134,A131,_:0.5743	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0A8M3-F1-model_v4.pdb	22.3594	2	2	2	0	0	1	642	94.2909	A195,A225,_:0.9356	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P76235-F1-model_v4.pdb	22.2102	2	2	2	0	0	1	427	79.4287	_,A267,A42:0.1757	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P00370-F1-model_v4.pdb	22.0781	2	2	2	0	0	1	447	95.8961	_,A161,A86:0.6623	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P77437-F1-model_v4.pdb	21.6085	2	2	2	0	0	1	526	92.8459	_,A309,A358:0.6906	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P63389-F1-model_v4.pdb	21.0561	2	2	2	0	0	1	637	84.6540	_,A422,A393:0.5600	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0A8F8-F1-model_v4.pdb	20.8974	2	2	2	0	0	1	673	88.5640	_,A410,A43:0.7293	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P37330-F1-model_v4.pdb	20.6907	2	2	2	0	0	1	723	96.6547	_,A315,A328:0.2507	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P77172-F1-model_v4.pdb	20.5964	2	2	2	0	0	1	747	87.5405	_,A492,A543:0.1283	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P75780-F1-model_v4.pdb	20.5466	2	2	2	0	0	1	760	90.9388	_,A61,A395:0.7400	B57,B102,C195	query/4CHA.pdb	e_coli_db/pdbtr/d16a4/index

[INFO] Querying query/1G2F.pdb:F207,F212,F225,F229 to e_coli_db/pdbtr/d16a4/index
[INFO] Found 1 index file(s). Querying with 1 threads
[INFO] SimpleHashMap::load_from_disk: 297.596µs
[INFO] load_lookup_from_file: 4.404964ms
[INFO] make_query_map: 1.754608ms
[INFO] read_u16_vector: 413.015µs
[INFO] count_query_idmode: 199.667µs
[INFO] queried_from_indices.par_sort_by: 10.35µs

id	idf_score	total_match_count	node_count	edge_count	exact_match_count	overflow_count	grid_count	nres	plddt	matching_residues	query_residues	query_file	index_path
model_benchmark/e_coli/pdb/AF-P05020-F1-model_v4.pdb	50.9269	4	3	4	0	0	1	348	97.0979	_,A222,A178,_:0.4606;_,_,A17,A19:0.3354	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P00957-F1-model_v4.pdb	48.7694	4	3	4	2	0	1	876	90.7229	_,A666,A564,A568:0.2691	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P64631-F1-model_v4.pdb	34.4415	2	2	2	2	0	1	117	96.7825	_,A73,A86,_:0.3155	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P76027-F1-model_v4.pdb	31.3890	2	2	2	0	0	1	337	90.4659	A185,_,A134,_:0.2929	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P06611-F1-model_v4.pdb	30.2622	2	2	2	2	0	1	249	96.7035	A180,_,_,A201:0.2231	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P39336-F1-model_v4.pdb	29.7054	2	2	2	0	0	1	604	64.3381	_,A386,A402,_:0.5127	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P23367-F1-model_v4.pdb	29.6533	2	2	2	0	0	1	615	81.9514	_,_,A297,A313:0.3344	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P76323-F1-model_v4.pdb	29.5204	2	2	2	2	0	1	92	74.5873	_,_,A84,A88:0.1495	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P77161-F1-model_v4.pdb	28.6326	2	2	2	0	0	1	292	97.7363	_,A259,A235,_:0.4761	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0AFD1-F1-model_v4.pdb	27.8174	2	2	2	2	0	1	166	93.2819	_,_,A24,A28:0.1244	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P45570-F1-model_v4.pdb	27.6982	2	2	2	2	0	1	173	95.9589	_,_,A57,A62:0.3499	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P77671-F1-model_v4.pdb	27.3656	2	2	2	0	0	1	453	97.9395	_,A287,A242,_:0.4745	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P33195-F1-model_v4.pdb	26.3775	2	2	2	2	0	1	957	96.7614	A624,_,_,A657:0.2361	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P16525-F1-model_v4.pdb	26.0246	2	2	2	2	0	1	309	95.3719	_,_,A27,A31:0.1273	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0A8Q3-F1-model_v4.pdb	25.8967	2	2	2	0	0	1	119	96.6612	_,_,A84,A88:0.2962	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0C8J8-F1-model_v4.pdb	25.1390	2	2	2	2	0	1	420	96.1708	_,_,A186,A189:0.5203	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0AFK0-F1-model_v4.pdb	24.9399	2	2	2	2	0	1	450	96.7774	_,_,A349,A296:0.6387	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0A6K3-F1-model_v4.pdb	24.8846	2	2	2	0	0	1	169	97.1330	_,_,A133,A137:0.3534	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P04805-F1-model_v4.pdb	24.8083	2	2	2	2	0	1	471	94.6752	_,_,A322,A326:0.1492	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P75914-F1-model_v4.pdb	23.8131	2	2	2	0	0	1	245	97.6758	_,_,A7,A9:0.3130	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0A908-F1-model_v4.pdb	23.7780	2	2	2	0	0	1	248	88.7910	_,_,A97,A101:0.3019	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P77766-F1-model_v4.pdb	23.2968	2	2	2	0	0	1	293	94.4010	_,_,A13,A15:0.2955	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P33591-F1-model_v4.pdb	23.0971	2	2	2	0	0	1	314	90.3978	_,_,A309,A313:0.2845	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P77302-F1-model_v4.pdb	22.3274	2	2	2	0	0	1	410	90.9846	_,_,A191,A195:0.3709	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P76641-F1-model_v4.pdb	22.1302	2	2	2	0	0	1	439	96.7650	_,_,A82,A84:0.3604	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P15034-F1-model_v4.pdb	22.1171	2	2	2	0	0	1	441	98.6170	_,_,A72,A74:0.3752	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-Q46812-F1-model_v4.pdb	22.1105	2	2	2	0	0	1	442	97.4955	_,_,A62,A64:0.3049	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0AEH1-F1-model_v4.pdb	22.0588	2	2	2	0	0	1	450	92.9109	_,_,A22,A26:0.3413	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-Q46806-F1-model_v4.pdb	21.9891	2	2	2	0	0	1	461	97.5706	_,_,A59,A61:0.3441	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P33919-F1-model_v4.pdb	21.2968	2	2	2	0	0	1	586	89.6611	_,_,A483,A485:0.3252	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P0AFA2-F1-model_v4.pdb	21.2383	2	2	2	0	0	1	598	83.9818	_,_,A494,A498:0.2825	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P76562-F1-model_v4.pdb	20.9060	2	2	2	0	0	1	671	90.5473	_,_,A308,A310:0.4780	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-Q46814-F1-model_v4.pdb	19.8846	2	2	2	0	0	1	956	96.3027	_,_,A362,A364:0.3246	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P21513-F1-model_v4.pdb	19.5839	2	2	2	0	0	1	1061	65.7316	_,_,A483,A485:0.3404	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P38394-F1-model_v4.pdb	16.2837	1	2	1	1	0	1	56	87.1077	A12,A7,_,_:0.2372	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P45771-F1-model_v4.pdb	14.5993	1	2	1	1	0	1	180	82.3846	A173,A168,_,_:0.3801	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P06612-F1-model_v4.pdb	12.3345	1	2	1	1	0	1	865	90.2204	A689,A683,_,_:0.3782;A736,A731,_,_:0.3956	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
model_benchmark/e_coli/pdb/AF-P39393-F1-model_v4.pdb	12.2677	1	2	1	1	0	1	906	86.2010	A124,A119,_,_:0.2776	F207,F212,F225,F229	query/1G2F.pdb	e_coli_db/pdbtr/d16a4/index
```