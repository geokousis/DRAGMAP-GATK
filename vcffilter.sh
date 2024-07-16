#!/bin/bash
REFERENCE="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
vcf="/media/storage/kousis/meth/Dra-GATK/SNPs.vcf"

filter_vcf(){
    local input_vcf="$1"
    local base_name=$(basename "$input_vcf")
    mkdir -p vcf tmp_v  # Ensure the directories exist
    
    # Apply the initial filter
    gatk VariantFiltration \
         -V "$input_vcf" \
         --filter-expression "DP < 600" \
         --filter-name "LowDP" \
         -O vcf/filtered_"${base_name}" \
         --create-output-variant-index \
         --tmp-dir tmp_v

    # Select only the variants that passed the filter
    gatk SelectVariants \
         -V vcf/filtered_"${base_name}" \
         --exclude-filtered \
         -O vcf/filtered_passed_"${base_name}" \
         --create-output-variant-index \
         --tmp-dir tmp_v

    # Split multi-allelic sites
    bcftools norm -m - vcf/filtered_passed_"${base_name}" > vcf/filtered_passed_split_"${base_name}"
}

# Call the function with the vcf variable
filter_vcf "$vcf"
