# Serine protease
-p query/4CHA.pdb -q B57,B102,C195

# Zinc finger
 -p query/1G2F.pdb -q F207,F212,F225,F229

# Aminopeptidase
-p query/1LAP.pdb -q 250,255,273,332,334

# Alkaline phosphatase
-p query/1ED8.pdb -q 51,102,153,155,166,322,327,328,331,369,370,412

# Commands to build index
## Swissprot
```sh
# PDBTRrosetta, grid mode, distance 16, angle 4, 128 threads
\time -v ~/Projects/06_Motifsearch/motifsearch/target/release/folddisco index -p swissprot_benchmark/swissprot_v4_raw -i swissprot_benchmark/pdbtr/grid/d16a4 -t 128 -y pdbtr -d 16 -a 4 -v --id uniprot -m grid
```

## PDB
```sh
\time -v ~/Projects/06_Motifsearch/motifsearch/target/release/folddisco index -p <> -i <> -t 128 -y pdbtr -d 16 -a 4 -v --id pdb -m grid
```

# Commands to query
```sh
~/Projects/06_Motifsearch/motifsearch/target/release/folddisco query -p query/20240426_WHY.pdb -i scop_benchmark/pdbtr/d16a4/index -t 32 -d 0.5 -a 5 > ~/temp.tsv
```