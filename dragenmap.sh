#!/bin/bash

start_time=$(date +"%Y-%m-%d %H:%M:%S")
cowsay -f dragon "Welcome To The DRAGMAP mapping script"
echo -e "\nThis section of the script performs mapping (+indexing, sort, flagstat) on a given list of files using dragen-os mapping aka DRAGMAP\n\n\n\n"

T_N=20
reference="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
ans="no"
samples="/home/kousis/work/meth/samples.txt"

total_lines=$(wc -l < "$samples")
mkdir -p hash_files
# Hash Table Dragen
if [ "$ans" == "yes" ]; then
    echo "Building Hash-table..."
    dragen-os --build-hash-table true --ht-reference $reference  --num-threads $T_N --output-directory hash_files 
  if [ $? -ne 0 ]; then
    echo "Hashing failed"
    exit 1
  fi
  echo "Hash Table created"
fi

# Initialize counter
counter=1
# DRAGMAP
while IFS= read -r line1 && IFS= read -r line2; do
    echo "Processing file $counter of $((total_lines / 2)): $line1 & $line2"

    # BWA MEM2
    filename=$(basename "$line1" | sed 's/_.*//')
    echo "Mapping with DRAGMAP on $filename"
    dragen-os -r hash_files --RGID $filename --RGSM $filename -1 "$line1" -2 "$line2" --num-threads $T_N > "${filename}.sam"
    if [ $? -ne 0 ]; then
      echo "DRAGMAP mapping failed for $filename"
      exit 1
    fi
    cowsay "DRAGMAP done"

    # SAM Files Processing
    samtools view -@ "$T_N" -b "${filename}.sam" > "${filename}.bam"
    samtools sort -@ "$T_N" "${filename}.bam" -o "sorted_${filename}.bam"
    samtools flagstat -@ "$T_N" "sorted_${filename}.bam" > "LOG_${filename}.txt"

    # Remove temporary data
    rm "${filename}.sam" "${filename}.bam"

    counter=$((counter + 1))
done < "$samples"

# File management
mkdir -p BWA_LOGS BWA_files
mv LOG_* BWA_LOGS
mv sorted_*.bam BWA_files

# Calculate and print the total execution time
end_time=$(date +"%Y-%m-%d %H:%M:%S")
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
runtime_seconds=$((end_seconds - start_seconds))
echo "Total runtime: $((runtime_seconds / 3600)) hours $(((runtime_seconds / 60) % 60)) minutes $((runtime_seconds % 60)) seconds"


