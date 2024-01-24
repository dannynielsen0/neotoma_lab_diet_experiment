#!/bin/bash

# Check if the input directory argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

# Set the input directory and output manifest file
input_directory="$1"
output_manifest="manifest.txt"

# Change to the input directory
cd "$input_directory" || exit

# Remove existing manifest file if it exists
rm -f "$output_manifest"

# Create an array to track unique sample IDs
declare -a unique_samples

# Get the absolute path of the input directory
absolute_path=$(realpath "$input_directory")

# Add column headers to the manifest file
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$output_manifest"


# Loop through the forward files and create the manifest
for forward_file in *_R1.fastq; do
    sample_id=$(echo "$forward_file" | cut -d'_' -f1)
    reverse_file=$(echo "$forward_file" | sed 's/_R1/_R2/')

    # Get the absolute paths of forward and reverse files
    forward_path="$absolute_path/$forward_file"
    reverse_path="$absolute_path/$reverse_file"

    # Check if sample ID is already in the manifest
    if [[ ! " ${unique_samples[@]} " =~ " $sample_id " ]]; then
        echo -e "$sample_id\t$forward_path\t$reverse_path" >> "$output_manifest"
        unique_samples+=("$sample_id")
    fi
done

echo "Manifest file created: $output_manifest"

