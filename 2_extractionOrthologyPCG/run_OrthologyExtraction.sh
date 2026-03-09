#!/bin/bash

set -euo pipefail

# Usage:
# bash run_OrthologyExtraction.sh [ensembl_version]

CONFIG_FILE="../data/config.txt"

ENSEMBL_VERSION=${1:-}

echo "========================================="
echo " PCG ORTHOLOGY EXTRACTION FROM BIOMART"
echo "========================================="

# Check config file

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: config file not found: $CONFIG_FILE"
    exit 1
fi


# Clean results folder

rm -rf results
mkdir -p results


# Run R script

if [ -z "$ENSEMBL_VERSION" ]; then
    
    echo "Running with latest Ensembl version"
    Rscript orthologyFromBiomart.R
    
else
    
    echo "Running with Ensembl version: $ENSEMBL_VERSION"
    Rscript orthologyFromBiomart.R "$ENSEMBL_VERSION"
    
fi


echo ""
echo "Orthology extraction completed"