# File: quick_start.sh
# Created: 2025-08-07 02:59:59
# Author: Hyunbin Kim (khb7840@gmail.com)
# Description:
#     Shell script to quickly build and run the folddisco tool with example data.
# Copyright © 2025 Hyunbin Kim, All rights reserved

# Build
cargo build --release --features foldcomp
echo "[INFO] Build complete. The folddisco binary is located at target/release/folddisco"
# Start with example data
echo "[INFO] Indexing data/serine_peptidases to index/serine_peptidases_folddisco"
target/release/folddisco index -p data/serine_peptidases -i index/serine_peptidases_folddisco
echo "[INFO] Indexing complete."
echo "[INFO] Search query/serine_peptidase.txt against the index built"

# Query a motif against the indexed serine peptidases
target/release/folddisco query -i index/serine_peptidases_folddisco -q query/serine_peptidase.txt
echo "[INFO] Querying complete. Results are printed to the console."
