#!/bin/bash

# RNAseq_preprocessing.sh - A script to preprocess RNAseq data. Last edited by Allison on 2025/03/21
# Usage: RNAseq_preprocessing.sh list_of_fastq_paths.txt [/path/to/parent_dir] [nproc_per_sample] [strandness]
# list_of_fastq_paths.txt: A file with full path to all R1 fastq.gz or fq.gz samples containing _R1. or _1.
# /path/to/parent_dir: Directory for analysis, default=current directory
# nproc_per_sample: Number of cores per sample, default=1
# strandness: Strandness, default=unset (0=unstranded, 1=forward, 2=reverse)

# Print help message
usage="Usage: $(basename "$0") list_of_fastq_paths.txt [/path/to/parent_dir] [nproc_per_sample] [strandness]

list_of_fastq_paths.txt: A file with full path to all R1 fastq.gz or fq.gz samples containing _R1. or _1.
/path/to/parent_dir: Directory for analysis, default=current directory
nproc_per_sample: Number of cores per sample, default=1
strandness: Strandness, default=unset (0=unstranded, 1=forward, 2=reverse)
-h: Display this help message"

# Check if the help option is invoked
if [ "$1" == "-h" ]; then
    echo "
    $usage"
    exit 0
fi

# Input arguments
fastq_list=$1 # FIRST ARGUMENT, file with list of full path to all samples ending in _R1/fastq or .fastq.gz
parent_dir=${2:-$(pwd)} # SECOND ARGUMENT, directory for analysis, default=current directory
nproc_per_sample=${3:-"1"} # THIRD ARGUMENT, number of cores per sample, default=1
strandness=${4:-"100"} # FOURTH ARGUMENT, strandness, default=unset (0=unstranded, 1=forward, 2=reverse)

# Check if the list of fastq paths file is provided
if [ "$#" -lt 1 ]; then
    echo "Error: No fastq paths file provided.
    
    $usage"
    exit 1
fi

# Check if the list of fastq paths file exists
if [ ! -f "$fastq_list" ]; then
    echo "Error: The fastq paths file does not exist.
    
    $usage"
    exit 1
fi

# Count the number of lines in the list_of_fastq_paths.txt
total_files=$(wc -l < "$fastq_list")
batch_size=20

if [ "$total_files" -lt "$batch_size" ]; then
    /home/shared/pipelines/RNAseq_preprocessing_PE_mouse_testing.sh "$fastq_list" "$parent_dir" "$nproc_per_sample"
    exit 0
fi

batch_count=$(( (total_files + batch_size - 1) / batch_size ))  # Round up division

# Iterate over each batch
for ((batch_num=1; batch_num<=batch_count; batch_num++)); do
    # Generate the file for this batch
    batch_file=$parent_dir/"$(basename "$fastq_list" .txt)_batch${batch_num}.txt"

    # Calculate the line range for this batch
    start_line=$(( (batch_num - 1) * batch_size + 1 ))
    end_line=$(( batch_num * batch_size ))
    if [ "$end_line" -gt "$total_files" ]; then
        end_line=$total_files
    fi

    # Extract the corresponding lines for this batch and save them to a new file
    sed -n "${start_line},${end_line}p" "$fastq_list" > "$batch_file"

    # Run RNAseq_preprocessing.sh for this batch
    echo "Processing batch ${batch_num} with ${batch_size} files (or less for last batch)..."
    /home/shared/pipelines/RNAseq_preprocessing_PE_mouse_testing.sh "$batch_file" "$parent_dir" "$nproc_per_sample" "$batch_num" "$strandness"
done

echo "All batches processed."
