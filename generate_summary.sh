#!/bin/bash

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <outprefix_1> <outprefix_2> [<outprefix_3> ...]"
    exit 1
fi

# Extract header
cols="${1}_summary.txt"
header=$(head -n 1 "$cols")

# Create summary file
echo "$header" > all_summary.txt
for test in "$@"; do
    sed -n '2p' "${test}_summary.txt" >> all_summary.txt
done
