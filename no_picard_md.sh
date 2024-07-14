#!/bin/bash

# Variables
bwa_files="/home/kousis/work/meth/Dra-GATK/BWA_files/"
N_T=20

# Enable alias expansion
shopt -s expand_aliases
alias picard='java -jar /home/kousis/work/tools/picard/build/libs/picard.jar'

# Change to the directory with BAM files
pushd $bwa_files > /dev/null
ls -1 -d $(pwd)/*.bam > tmp.txt
popd > /dev/null

# Function to mark duplicates using Picard
MarkDuplicates() {
    local input_bam="$1"
    base_name=$(basename "$input_bam" | sed 's/^sorted_//')
    output_bam="md_${base_name}"
    metrics_file="marked_dup_metrics_${base_name%.bam}.txt"

    echo "Processing $input_bam -> $output_bam"

    picard MarkDuplicates \
        CREATE_INDEX=true \
        QUIET=true \
        INPUT="$input_bam" \
        OUTPUT="$bwa_files/$output_bam" \
        METRICS_FILE="$bwa_files/$metrics_file"

    if [ $? -ne 0 ]; then
        echo "Error processing $input_bam"
        exit 1
    fi
}

# Record start time
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "This section of the script performs Duplication marking with Picard's MarkDuplicates"

# Process each BAM file
count=0
while read -r bam_file; do
    # Control the number of parallel executions
    while [ "$(jobs -p | wc -l)" -ge "$N_T" ]; do
        sleep 1
    done

    MarkDuplicates "$bam_file" &
    count=$((count + 1))
    echo "Processed $count BAM files"
done < "${bwa_files}/tmp.txt"
wait

echo "All BAM files have been marked"

# Calculate and print the total execution time
end_time=$(date +"%Y-%m-%d %H:%M:%S")
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
runtime_seconds=$((end_seconds - start_seconds))

echo "Total runtime: $((runtime_seconds / 3600)) hours $(((runtime_seconds / 60) % 60)) minutes $((runtime_seconds % 60)) seconds"
