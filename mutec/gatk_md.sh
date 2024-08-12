#!/bin/bash
start_time=$(date +"%Y-%m-%d %H:%M:%S")

echo -e "\033[1mStarting Duplication marking process with GATK Spark implementation of Picard's MarkDuplicates.\033[0m"
echo -e "\033[1mBAM list ---> QuerySort BAMs ---> MarkDuplicates ---> SetNmMdAndUqTags\033[0m"

# Variables
bwa_files="/home/kousis/work/meth/Dra-GATK/mutec/BWA_files/"
N_T=28
reference="/home/kousis/work/meth/Dra-GATK/mutec/reference/Homo_sapiens_assembly38.fasta"
mkdir -p md_tmp qs_tmp fix_tmp logs_mark_dup

# Function to sort BAM files by queryname
QuerySort() {
    local input_bam="$1"
    local base_name=$(basename "$input_bam" | sed 's/^sorted_//')
    local output_bam="${bwa_files}qs_${base_name}"
    local log_file="logs_mark_dup/qs_${base_name}.log"

    echo "Query Sorted $(basename "$input_bam") -> $(basename "$output_bam")" | tee "$log_file"

    # Running the Picard SortSam command
    java -jar /home/kousis/work/tools/picard/build/libs/picard.jar SortSam \
        I="$input_bam" \
        O="$output_bam" \
        SORT_ORDER=queryname \
        TMP_DIR=qs_tmp \
        &>> "$log_file"

    if [ $? -ne 0 ]; then
        echo "Error processing $input_bam" | tee -a "$log_file"
        exit 1
    fi
}

# Function to mark duplicates using GATK (runs sequentially)
MarkDuplicates() {
    local input_bam="$1"
    local base_name=$(basename "$input_bam" | sed 's/^qs_//')
    local output_bam="${bwa_files}md_${base_name}"
    local metrics_file="${bwa_files}marked_dup_metrics_${base_name%.bam}.txt"
    local log_file="logs_mark_dup/md_${base_name}.log"

    echo "Marking duplicates in $(basename "$input_bam") -> $(basename "$output_bam")" | tee "$log_file"

    gatk MarkDuplicatesSpark \
        -I "$input_bam" \
        -O "$output_bam" \
        -M "$metrics_file" \
        --tmp-dir md_tmp \
        --create-output-bam-splitting-index true \
        --spark-runner LOCAL \
        --conf "spark.executor.cores=${N_T}" \
        &>> "$log_file"

    if [ $? -ne 0 ]; then
        echo "Error processing $input_bam" | tee -a "$log_file"
        exit 1
    fi
}

# Function to set NM, MD, and UQ tags
SetNmMdAndUqTags() {
    local input_bam="$1"
    local base_name=$(basename "$input_bam" | sed 's/^md_//')
    local output_bam="${bwa_files}f_${base_name}"
    local log_file="logs_mark_dup/f_${base_name}.log"

    echo "Set NM, MD, and UQ tags in $(basename "$input_bam") -> $(basename "$output_bam")" | tee "$log_file"
    java -jar /home/kousis/work/tools/picard/build/libs/picard.jar SetNmMdAndUqTags \
        R="$reference" \
        I="$input_bam" \
        O="$output_bam" \
        TMP_DIR=fix_tmp \
        CREATE_INDEX=true \
        &>> "$log_file"

    if [ $? -ne 0 ]; then
        echo "Error processing $input_bam" | tee -a "$log_file"
        exit 1
    fi
}

# Export functions and variables for parallel execution
export bwa_files
export N_T
export reference
export -f QuerySort
export -f SetNmMdAndUqTags

# Sorting
echo -e "\033[1mQueryname Sorting...\033[0m"
# Creating list for execution
pushd "$bwa_files" > /dev/null
ls -1 -d "$(pwd)"/*.bam > bam_list.txt
popd > /dev/null

# Query sort all BAM files in parallel
#parallel --jobs "$N_T" QuerySort :::: "${bwa_files}/bam_list.txt"
echo -e "\033[1mQueryname Sorting Completed\033[0m"

# Mark dup
echo -e "\033[1mMarking Duplicates...\033[0m"
# Creating list for execution
pushd "$bwa_files" > /dev/null
ls -1 -d "$(pwd)"/qs*.bam > qs_list.txt
popd > /dev/null

# Mark duplicates sequentially (one by one)
# while read -r bam_file; do
#     MarkDuplicates "$bam_file"
#     echo "Marked duplicates for $(basename "$bam_file")"
# done < "${bwa_files}/qs_list.txt"
echo -e "\033[1mMarking Duplicates Completed\033[0m"

# Fix Bam according to GATK recommendation
echo -e "\033[1mFixing BAMs...\033[0m"
# Creating list for execution
pushd "$bwa_files" > /dev/null
ls -1 -d "$(pwd)"/md*.bam > md_list.txt
popd > /dev/null

# Set NM, MD, and UQ tags in parallel
parallel --jobs "$N_T" SetNmMdAndUqTags :::: "${bwa_files}/md_list.txt"
echo -e "\033[1mFixing BAMs Completed\033[0m"

# Move initial list log to logs directory
mv "${bwa_files}/bam_list.txt" "${bwa_files}/md_list.txt" "${bwa_files}/qs_list.txt" logs_mark_dup/
echo "Logs can be found in the ${bwa_files} directory, processed files start with f_"

# Clean up temporary directories and environment
rm -r md_tmp qs_tmp fix_tmp
unset bwa_files
unset N_T
unset reference

# Script Duration
end_time=$(date +"%Y-%m-%d %H:%M:%S")
start_seconds=$(date -d "$start_time" +%s)
end_seconds=$(date -d "$end_time" +%s)
duration=$(( end_seconds - start_seconds ))
echo "Duplication marking completed in $duration seconds."
echo -e "\033[1mScript ended\033[0m"
