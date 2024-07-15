#!/bin/bash

bwa_files="/home/kousis/work/meth/Dra-GATK/BWA_files/"
NUM_CORES=25
GFF_FILE="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/genomic.gff"
REFERENCE="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"

# Create chromosome list
output_file="chromosome_list.list"
grep "^>" "$REFERENCE" | awk '{print substr($1, 2)}' | sort | uniq | grep '^NC' > "$output_file"

pushd "$bwa_files" > /dev/null
ls -1 -d "$(pwd)"/md_*.bam > tmp.txt
popd > /dev/null

awk '$3 == "exon"' "$GFF_FILE" | awk '{OFS="\t"; print $1, $4-1, $5}' | sort -k1,1 -k2,2n -u > tmp.exons.bed
echo "BED created"

TOTAL_LINES=$(wc -l < tmp.exons.bed)
SPLIT_SIZE=$(( (TOTAL_LINES + NUM_CORES - 1) / NUM_CORES ))

python3 fix_bed.py
rm tmp.exons.bed
mkdir -p tmp_v
mkdir -p vcf

Data_Pre_Processing() {
    local reference="$1"
    gatk ComposeSTRTableFile \
         -R "$reference" \
         -O vcf/str_table.tsv
    if [ $? -ne 0 ]; then
        echo "Error during Data pre processing"
        exit 1
    fi
}

Haplotype_caller() {
    local reference="$1"
    local input_bam="$2"
    local base_name=$(basename "$input_bam" | sed 's/^md_//')
    gatk CalibrateDragstrModel \
         -R "$reference" \
         -I "$input_bam" \
         -str vcf/str_table.tsv \
         -O vcf/"${base_name}_dragstr_model.txt"
    if [ $? -ne 0 ]; then
        echo "Error during Calibration"
        exit 1
    fi

    gatk HaplotypeCaller \
         -R "$reference" \
         -I "$input_bam" \
         -L merged_exons.bed \
         --interval-padding 100 \
         -O vcf/"${base_name}.g.vcf" \
         -ERC GVCF \
         --tmp-dir tmp_v \
         --create-output-variant-index \
         --dragen-mode true \
         --dragstr-params-path vcf/"${base_name}_dragstr_model.txt"
    if [ $? -ne 0 ]; then
        echo "Error during call"
        exit 1
    fi

    # Uncomment if per-sample filtration is needed
    # gatk VariantFiltration \
    #      -V vcf/"${base_name}.g.vcf" \
    #      --filter-expression "QUAL < 50" \
    #      --filter-name "DRAGENHardQUAL" \
    #      --filter-expression "DP < 60" \
    #      --filter-name "LowDP" \
    #      -O vcf/filtered_"${base_name}.g.vcf" \
    #      --create-output-variant-index \
    #      --tmp-dir tmp_v
}

Merge() {
    local reference="$1"
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
         -L chromosome_list.list \
         --genomicsdb-workspace-path vcf/my_database \
         --tmp-dir tmp_v \
         --sample-name-map vcf/cohort.sample_map \
         --create-output-variant-index
    if [ $? -ne 0 ]; then
        echo "Error during GDB"
        exit 1
    fi
    gatk --java-options "-Xmx4g" GenotypeGVCFs \
         -R "$reference" \
         -V gendb://vcf/my_database \
         -O vcf/Merged.vcf \
         --create-output-variant-index
    if [ $? -ne 0 ]; then
        echo "Error during call"
        exit 1
    fi

    # Select SNPs from the merged VCF
    gatk SelectVariants \
         -R "$reference" \
         -V vcf/Merged.vcf \
         --select-type-to-include SNP \
         -O SNPs.vcf
    if [ $? -ne 0 ]; then
        echo "Error during Variant Selection SNPs"
        exit 1
    fi

    # Select INDELs from the merged VCF
    gatk SelectVariants \
         -R "$reference" \
         -V vcf/Merged.vcf \
         --select-type-to-include INDEL \
         -O INDELs.vcf
    if [ $? -ne 0 ]; then
        echo "Error during Variant Selection Indels"
        exit 1
    fi
}

# Create cohort.sample_map file
create_sample_map() {
    while read -r bam_file; do
        base_name=$(basename "$bam_file" | sed 's/^md_//')
        echo -e "${base_name}\tvcf/${base_name}.g.vcf"
    done < "${bwa_files}/tmp.txt" > vcf/cohort.sample_map
}

# Data pre-processing
Data_Pre_Processing "$REFERENCE"

# Process each BAM file
count=0
while read -r bam_file; do
    running_jobs=$(jobs -p | wc -l)
    while [ "$running_jobs" -ge "$NUM_CORES" ]; do
        sleep 1
        running_jobs=$(jobs -p | wc -l)
    done
    Haplotype_caller "$REFERENCE" "$bam_file" &
    count=$((count + 1))
done < "${bwa_files}/tmp.txt"

wait

# Create sample map after all HaplotypeCaller jobs are done
create_sample_map

# Merge the GVCF files
Merge "$REFERENCE"

echo "All BAM files have been processed and merged."
