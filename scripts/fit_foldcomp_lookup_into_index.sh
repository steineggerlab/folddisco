#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <lookup_file> <index_file> <subset_file>"
  exit 1
fi

# Assign arguments to variables
lookup_file="$1"
index_file="$2"
subset_file="$3"

# Use awk to extract matched rows from lookup based on index file
echo "Extracting $lookup_file based on $index_file..."
awk 'NR==FNR {index_ids[$1]; next} $1 in index_ids' "$index_file" "$lookup_file" > "$subset_file"
echo "Subset extracted"
# Sort the output based on the second column alphabetically
echo "Sorting subset..."
sort -t$'\t' -k2,2 "$subset_file" -o "$subset_file"
echo "Subset sorted"
# Print success message
echo "Subset saved to $subset_file"
