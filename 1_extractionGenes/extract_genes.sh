#!/bin/bash

# Usage:
# bash extract_genes.sh config.txt

configFile_path=${1:-../data/config.txt}

echo "Using config file: $configFile_path"

# Check if config file exists
if [ ! -f "$configFile_path" ]; then
    echo "ERROR: Config file not found: $configFile_path"
    exit 1
fi

# Clean result directory
rm -rf results
mkdir -p results

# Loop over GTF files listed in config file (column 4)
tail -n +2 "$configFile_path" | cut -f4 | while read GTF
do
    echo ""
    echo "Processing: $GTF"

    # Check if GTF exists
    if [ ! -f "$GTF" ]; then
        echo "WARNING: File not found -> $GTF"
        echo "Skipping..."
        continue
    fi

    # Extract filename
    base=$(basename "$GTF" .gtf)
    output_name="${base}_genesOnly.gtf"

    # Extract gene features
    awk -F"\t" '
    $0 !~ /^#/ && $3=="gene" {print}
    ' "$GTF" > "results/$output_name"

    echo "Output written to results/$output_name"

done

echo ""
echo "Gene extraction completed."

echo "Running R script: extract_geneInfo.R"
Rscript extract_geneInfo.R