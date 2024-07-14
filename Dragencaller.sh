#!/bin/bash

bwa_files="/home/kousis/work/meth/Dra-GATK/BWA_files/"
NUM_CORES=20
GFF_FILE="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/genomic.gff"
REFERENCE="/home/kousis/work/meth/annotation/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"

pushd "$bwa_files" > /dev/null
ls -1 -d "$(pwd)"/md_*.bam > tmp.txt
popd > /dev/null

awk '$3 == "exon"' $GFF_FILE | awk '{OFS="\t"; print $1, $4-1, $5}' | sort -k1,1 -k2,2n -u > tmp.exons.bed
echo "BED created"

TOTAL_LINES=$(wc -l < tmp.exons.bed)
SPLIT_SIZE=$(( (TOTAL_LINES + NUM_CORES - 1) / NUM_CORES ))

python3 fix_bed.py

mkdir -p tmp_v
mkdir -p vcf

Data_Pre_Processing(){
    local reference="$1"
    local input_reads=$(ls ${bwa_files}/md_*.bam | xargs -I {} echo -I {})
    gatk ComposeSTRTableFile \
         -R $reference \
         -O vcf/str_table.tsv

    gatk CalibrateDragstrModel \
         -R $reference \
         -I $input_reads \
         -str vcf/str_table.tsv \
         -O vcf/dragstr_model.txt
    if [ $? -ne 0 ]; then
        echo "Error during Data pre processing"
        exit 1
    fi
}

Haplotype_caller(){
    local reference="$1"
    local input_bam="$2"
    local base_name=$(basename "$input_bam" | sed 's/^md_//')
    gatk HaplotypeCaller \
         -R "$reference" \
         -I "$input_bam" \
         -L tmp.exons.bed \
         -O vcf/${base_name}.g.vcf \
         -ERC GVCF \
         --tmp-dir tmp_v \
         --create-output-variant-index \
         --interval-padding 100 \
         --dragen-mode true \
         --dragstr-params-path vcf/dragstr_model.txt

    gatk VariantFiltration \
         -V vcf/${base_name}.g.vcf \
         --filter-expression "QUAL < 10" \
         --filter-name "DRAGENHardQUAL" \
         -O vcf/filtered_${base_name}.g.vcf \
         --create-output-variant-index \
         --tmp-dir tmp_v
}

Merge(){
    local reference="$1"
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
         --genomicsdb-workspace-path vcf/my_database \
         --tmp-dir tmp_v \
         --sample-name-map vcf/cohort.sample_map \
         --create-output-variant-index

    gatk --java-options "-Xmx4g" GenotypeGVCFs \
         -R $reference \
         -V gendb://vcf/my_database \
         -O vcf/Merged.vcf \
         --create-output-variant-index
}

# Create cohort.sample_map file
create_sample_map() {
    while read -r bam_file; do
        base_name=$(basename "$bam_file" | sed 's/^md_//')
        echo -e "${base_name}\t${bwa_files}/vcf/filtered_${base_name}.g.vcf"
    done < "${bwa_files}/tmp.txt" > vcf/cohort.sample_map
}

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