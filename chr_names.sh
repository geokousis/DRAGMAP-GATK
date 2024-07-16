#!/bin/bash

# Define the annotation file
annotation_file="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/genomic.gff"
# Define the output file
output_file="chromosome_list.txt"
# Extract chromosome names and put them in a list
awk '{print $1}' "$annotation_file" | sort | uniq | grep '^NC' > "$output_file"

