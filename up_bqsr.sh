#!/bin/bash

# Variables
bwa_files="/home/kousis/work/meth/Dra-GATK/BWA_files/"
NUM_CORES=20
reference="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
db_SNP="/home/kousis/work/meth/annotation/known/ncbi_1000G_phase1.snps.high_confidence.hg38.vcf"
GFF_FILE="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/genomic.gff"
db_indels="/home/kousis/work/meth/annotation/known/ncbi_Mills_and_1000G_gold_standard.indels.hg38.vcf"
pushd "$bwa_files" > /dev/null
ls -1 -d "$(pwd)"/md_*.bam > tmp.txt
popd > /dev/null
mkdir -p md_tmp

# Step 1: Extract Exon Regions from GFF and Convert to BED
awk '$3 == "exon"' $GFF_FILE | awk '{OFS="\t"; print $1, $4-1, $5}' | sort -k1,1 -k2,2n -u > tmp.exons.bed
echo "BED created"

# Step 2: Calculate Optimal Split Size
TOTAL_LINES=$(wc -l < tmp.exons.bed)
SPLIT_SIZE=$(( (TOTAL_LINES + NUM_CORES - 1) / NUM_CORES ))

# Split BED File
python3 fix_bed.py
split -l $SPLIT_SIZE -d --additional-suffix=.bed merged_exons.bed exons_split_

rm tmp.exons.bed
rm merged_exons.bed

ls exons_split_*.bed > tmp.exons.list

# Function to create recalibration tables
BRC() {
    local input_bam="$1"
    local input_interval="$2"
    local base_name=$(basename "$input_bam" | sed 's/^sorted_//')
    local recal_table="tmp_${base_name}_${input_interval}.table"

    echo "Processing $input_bam -> $output_bam"

    gatk BaseRecalibrator \
        -I "$input_bam" \
        -R "$reference" \
        -L "$input_interval" \
        --known-sites "$db_SNP" \
        --known-sites "$db_indels" \
        -O "$recal_table" \
        --tmp-dir md_tmp

    if [ $? -ne 0 ]; then
        echo "Error during BaseRecalibrator for $input_bam"
        exit 1
    fi
}

# Function to Merge Recalibration Tables
MergeRT() {
    local input_bam="$1"
    local base_name=$(basename "$input_bam" | sed 's/^sorted_//')
    local output_table="${base_name}.table"

    # Find all matching table files and join them with -I flag
    local input_tables=$(ls tmp*${base_name}*.table | xargs -I {} echo -I {})
    
    gatk GatherBQSRReports $input_tables -O "$output_table"

    if [ $? -ne 0 ]; then
        echo "Error during Merging Tables for $output_table"
        exit 1
    fi
}


# Function to Apply BaseQualityScoreRecalibration
ApplyBQSR() {
    local input_bam="$1"
    local base_name=$(basename "$input_bam" | sed 's/^sorted_//')
    local input_table="${base_name}.table"
    local output_bam="Recalibrated_${base_name}"
    
    gatk --java-options "-XX:ConcGCThreads=1" ApplyBQSRSpark \
        -I "$input_bam" \
        -bqsr "$input_table" \
        -O "$output_bam" \
        -- \
        --spark-master local[15] \
        --conf 'spark.executor.cores=5' \
        --conf 'spark.executor.instances=3'
}

# Main script logic
start_time=$(date +"%Y-%m-%d %H:%M:%S")

echo "This section of the script performs BaseQualityScoreRecalibration with a slight variation of GATK best practices"

# Process each BAM file
count=0
while read -r bam_file; do
    while read -r exons_sites; do
        running_jobs=$(jobs -p | wc -l)
        while [ "$running_jobs" -ge "$NUM_CORES" ]; do
            sleep 1
            running_jobs=$(jobs -p | wc -l)
        done
        BRC "$bam_file" "$exons_sites" &
        count=$((count + 1))
    done < tmp.exons.list
    wait
    MergeRT "$bam_file"
    ApplyBQSR "$bam_file"
    echo "Processed $count BAM files"
done < "${bwa_files}/tmp.txt"

echo "All BAM files have been BQSR"

rm -r md_tmp
rm tmp.exons.list
rm tmp*.table
rm exons_split*
end_time=$(date +"%Y-%m-%d %H:%M:%S")
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
runtime_seconds=$((end_seconds - start_seconds))

echo "Total runtime: $((runtime_seconds / 3600)) hours $(((runtime_seconds % 3600) / 60)) minutes $((runtime_seconds % 60)) seconds"
