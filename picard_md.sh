#!/bin/bash

# Variables
bwa_files="/home/kousis/work/meth/Dra-GATK/BWA_files/"
N_T=20

pushd "$bwa_files" > /dev/null
ls -1 -d "$(pwd)"/*.bam > tmp.txt
popd > /dev/null
mkdir -p md_tmp
# Function to mark duplicates using GATK to use SPARK
MarkDuplicates() {
    local input_bam="$1"
    local base_name=$(basename "$input_bam" | sed 's/^sorted_//')
    local output_bam="md_${base_name}"
    local metrics_file="marked_dup_metrics_${base_name%.bam}.txt"

    echo "Processing $input_bam -> $output_bam"

    gatk --java-options "-XX:ConcGCThreads=1" MarkDuplicatesSpark \
        -I "$input_bam" \
        -O "$bwa_files/$output_bam" \
        -M "$bwa_files/$metrics_file" \
	--tmp-dir md_tmp \
        --create-output-bam-splitting-index true \
	--spark-master local[15] \
	--conf 'spark.executor.cores=5' \
	--conf 'spark.executor.instances=3'

    if [ $? -ne 0 ]; then
        echo "Error processing $input_bam"
        exit 1
    fi
}


start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "This section of the script performs Duplication marking with Picard's MarkDuplicates"

# Process each BAM file
count=0
while read -r bam_file; do
    MarkDuplicates "$bam_file"
    count=$((count + 1))
    echo "Processed $count BAM files"
done < "${bwa_files}/tmp.txt"

echo "All BAM files have been marked"

rm -r md_tmp
end_time=$(date +"%Y-%m-%d %H:%M:%S")
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
runtime_seconds=$((end_seconds - start_seconds))

echo "Total runtime: $((runtime_seconds / 3600)) hours $(((runtime_seconds % 3600) / 60)) minutes $((runtime_seconds % 60)) seconds"
