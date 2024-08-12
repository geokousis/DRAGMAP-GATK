#!/bin/bash
start_time=$(date +"%Y-%m-%d %H:%M:%S")

echo -e "\033[1mStarting Base Quality Score Recalibration (BQSR) process.\033[0m"
echo -e "\033[1mSplitting Intervals ---> BQSR (Parallel by Intervals) ---> CombineBQSRReports ---> ApplyBQSR (Parallel by BAM files)\033[0m"

# Variables
bwa_files="/home/kousis/work/meth/Dra-GATK/mutec/BWA_files/"
N_T=28
reference="/home/kousis/work/meth/Dra-GATK/mutec/reference/Homo_sapiens_assembly38.fasta"
db_SNP="/home/kousis/work/meth/Dra-GATK/mutec/res/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
thousend_SNP="/home/kousis/work/meth/Dra-GATK/mutec/res/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
exome_db="/home/kousis/work/meth/Dra-GATK/mutec/res/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
mills_db="/home/kousis/work/meth/Dra-GATK/mutec/res/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
knownindels_db="/home/kousis/work/meth/Dra-GATK/mutec/res/Homo_sapiens_assembly38.known_indels.vcf.gz"
intervals="/home/kousis/work/meth/Dra-GATK/exome_calling_regions.v1.unpadded.interval_list"

mkdir -p intervals logs_bqsr bqsr_tmp

# Function to split intervals
SplitIntervals() {
    echo -e "\033[1mSplitting Intervals...\033[0m"
    gatk SplitIntervals \
        -R "$reference" \
        -L "$intervals" \
	--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
        --scatter-count "$N_T" \
        -O intervals

    if [ $? -ne 0 ]; then
        echo "Error during interval splitting"
        exit 1
    fi
    echo -e "\033[1mInterval Splitting Completed\033[0m"
}

# Function to perform Base Quality Score Recalibration
BQSR() {
    local input_bam="$1"
    local interval_split="$2"
    local base_name=$(basename "$input_bam" .bam)
    local interval_code=$(basename "$interval_split" | sed 's/\.interval_list$//')
    local recal_table="${base_name}_${interval_code}.table"
    local log_file="logs_bqsr/bqsr_${base_name}_${interval_code}.log"

    echo "Running BQSR on $(basename "$input_bam") for interval $interval_code" | tee "$log_file"

    gatk BaseRecalibrator \
        -I "$input_bam" \
        -R "$reference" \
        -L "$interval_split" \
        --interval-padding 50 \
        --known-sites "$db_SNP" \
        --known-sites "$thousend_SNP" \
        --known-sites "$exome_db" \
        --known-sites "$mills_db" \
        --known-sites "$knownindels_db" \
        -O "$recal_table" \
        --tmp-dir bqsr_tmp \
        &>> "$log_file"

    if [ $? -ne 0 ]; then
        echo "Error during BQSR for $input_bam on $interval_split" | tee -a "$log_file"
        exit 1
    fi
}

# Function to combine the BQSR reports
CombineTables() {
    local base_name="$1"
    shift
    local recal_tables=("$@")
    local output_combined_table="${bwa_files}${base_name}_combined.table"
    local log_file="logs_bqsr/combine_${base_name}.log"
    local arguments_file="logs_bqsr/${base_name}_recal_tables.args"

    echo "Combining BQSR tables for $base_name" | tee "$log_file"

    # Write the list of recalibration tables to an arguments file
    for table in "${recal_tables[@]}"; do
        echo "-I $table" >> "$arguments_file"
    done

    gatk GatherBQSRReports \
        --arguments_file "$arguments_file" \
        -O "$output_combined_table" \
        &>> "$log_file"

    if [ $? -ne 0 ]; then
        echo "Error during CombineTables for $base_name" | tee -a "$log_file"
        exit 1
    fi

    echo "Combined BQSR tables into $output_combined_table" | tee -a "$log_file"

    # Clean up the arguments file
    rm "$arguments_file"
}

# Function to apply BQSR
ApplyBQSR() {
    local input_bam="$1"
    local base_name=$(basename "$input_bam" .bam)  # Correctly get the base name without the .bam extension
    local combined_table="${bwa_files}${base_name}_combined.table"
    local output_bqsr_bam="${bwa_files}${base_name}_recalibrated.bam"
    local log_file="logs_bqsr/applybqsr_${base_name}.log"

    echo "Applying BQSR to $(basename "$input_bam") using recalibration table $combined_table" | tee "$log_file"

    # Check if the recalibration table exists
    if [ ! -f "$combined_table" ]; then
        echo "Error: Recalibration table $combined_table not found!" | tee -a "$log_file"
        exit 1
    fi

    gatk ApplyBQSR \
        -R "$reference" \
        -I "$input_bam" \
        --bqsr-recal-file "$combined_table" \
        -O "$output_bqsr_bam" \
        &>> "$log_file"

    if [ $? -ne 0 ]; then
        echo "Error during ApplyBQSR for $input_bam" | tee -a "$log_file"
        exit 1
    fi

    echo "Applied BQSR to $input_bam, output to $output_bqsr_bam" | tee -a "$log_file"
}

# Variables
export bwa_files
export N_T
export reference
export db_SNP
export thousend_SNP
export exome_db
export mills_db
export knownindels_db
export intervals


# Step 1: Split Intervals (sequential, no parallelization)
SplitIntervals

# Step 2: Run BQSR in parallel across intervals for each BAM file
export -f BQSR
for bam_file in "${bwa_files}"/f_*.bam; do
    echo "Running BQSR in parallel for $(basename "$bam_file")"
    parallel --jobs "$N_T" BQSR ::: "$bam_file" ::: intervals/*.interval_list
done

# Step 3: Combine the BQSR tables for each BAM file
for bam_file in "${bwa_files}"/f_*.bam; do
    tables=()
    base_name=$(basename "$bam_file" .bam)
    for interval_split in intervals/*.interval_list; do
        interval_code=$(basename "$interval_split" | sed 's/\.interval_list$//')
        tables+=("${base_name}_${interval_code}.table")
    done
    CombineTables "$base_name" "${tables[@]}"
done

# Step 4: Apply BQSR in parallel across all BAM files with the prefix f_
export -f ApplyBQSR

echo "Applying BQSR in parallel for all BAM files with prefix f_"
parallel --jobs "$N_T" ApplyBQSR ::: "${bwa_files}"/f_*.bam

# Clean up temporary directories
rm -r bqsr_tmp intervals
rm f_*scattered.table


end_time=$(date +"%Y-%m-%d %H:%M:%S")
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
runtime_seconds=$((end_seconds - start_seconds))

echo "Total runtime: $((runtime_seconds / 3600)) hours $(((runtime_seconds % 3600) / 60)) minutes $((runtime_seconds % 60)) seconds"
