#!/usr/bin/env bash

###############################################################################
# Aim: Run Ensembl Compara queries for all species listed in config.txt
#
# Input:
#   ../data/config.txt
#   input_formatted/lncRNA_formatted_<species>.tsv
#
# Output:
#   results_raw/<ensemblName>_compara.tsv
#
# Usage:
#   bash run_compara_queries.sh [config_file]
###############################################################################

set -e
set -o pipefail

CONFIG_FILE=${1:-"../data/config.txt"}
PERL_SCRIPT="2_MP_customGene.pl"
INPUT_DIR="input_formatted"
OUTPUT_DIR="results_raw"

echo "========================================="
echo " COMPARA QUERIES - ALL SPECIES"
echo "========================================="

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: config file not found: $CONFIG_FILE"
    exit 1
fi

if [[ ! -f "$PERL_SCRIPT" ]]; then
    echo "ERROR: Perl script not found: $PERL_SCRIPT"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Using config file: $CONFIG_FILE"
echo "Using input directory: $INPUT_DIR"
echo "Using output directory: $OUTPUT_DIR"
echo "-----------------------------------------"

# Read config file line by line
# Expected columns:
# shortName   completeName   ensemblName   gtfPath   
while IFS=$'\t' read -r shortName completeName ensemblName gtfPath 
do
    # Skip header
    if [[ "$shortName" == "shortName" ]]; then
        continue
    fi

    input_file="${INPUT_DIR}/lncRNA_formatted_${ensemblName}.tsv"
    output_file="${OUTPUT_DIR}/${ensemblName}_compara.tsv"

    echo ""
    echo "Processing species: $completeName ($ensemblName)"

    if [[ ! -f "$input_file" ]]; then
        echo "WARNING: input file not found: $input_file"
        echo "Skipping..."
        continue
    fi

    perl "$PERL_SCRIPT" \
        --name "$ensemblName" \
        --input "$input_file" \
        --output "$output_file"

    echo "Output written to: $output_file"

done < "$CONFIG_FILE"

echo ""
echo "-----------------------------------------"
echo "All Compara queries completed"
echo "-----------------------------------------"