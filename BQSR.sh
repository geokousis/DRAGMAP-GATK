#!/bin/bash

# Variables
bwa_files="/home/kousis/work/meth/Dra-GATK/BWA_files/"
N_T=20
reference="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
db_SNP="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/SNP_db/ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.vcf"
pushd "$bwa_files" > /dev/null
ls -1 -d "$(pwd)"/md_*.bam > tmp.txt
popd > /dev/null
mkdir -p md_tmp
# Function to Recalibrate Bases using GATK-Spark
BQSR() {
    local input_bam="$1"
    local base_name
    base_name=$(basename "$input_bam" | sed 's/^sorted_//')
    local output_bam="md_${base_name}"
    local metrics_file="marked_dup_metrics_${base_name%.bam}.txt"
    local recal_table="${base_name}.table"
    local output_bqsr_bam="${base_name}_BQSR.bam"

    echo "Processing $input_bam -> $output_bam"

    gatk --java-options "-XX:ConcGCThreads=1" BaseRecalibratorSpark \
        -I "$input_bam" \
        -R "$reference" \
        --known-sites "$db_SNP" \
        -O "$recal_table" \
	--tmp-dir md_tmp \
        -- \
        --spark-master local[15] \
        --conf 'spark.executor.cores=5' \
        --conf 'spark.executor.instances=3'

    if [ $? -ne 0 ]; then
        echo "Error during BaseRecalibratorSpark for $input_bam"
        exit 1
    fi

    gatk --java-options "-XX:ConcGCThreads=1" ApplyBQSRSpark \
        -I "$input_bam" \
        -bqsr "$recal_table" \
        -O "$output_bqsr_bam" \
        -- \
        --spark-master local[15] \
        --conf 'spark.executor.cores=5' \
        --conf 'spark.executor.instances=3'

    if [ $? -ne 0 ]; then
        echo "Error during ApplyBQSRSpark for $input_bam"
        exit 1
    fi

    echo "Successfully processed $input_bam to $output_bqsr_bam"
}

start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "This section of the script performs Duplication marking with Picard's MarkDuplicates"

# Process each BAM file
count=0
while read -r bam_file; do
    BQSR "$bam_file"
    count=$((count + 1))
    echo "Processed $count BAM files"
done < "${bwa_files}/tmp.txt"

echo "All BAM files have been BQSR"

rm -r md_tmp
end_time=$(date +"%Y-%m-%d %H:%M:%S")
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
runtime_seconds=$((end_seconds - start_seconds))

echo "Total runtime: $((runtime_seconds / 3600)) hours $(((runtime_seconds % 3600) / 60)) minutes $((runtime_seconds % 60)) seconds"
