#!/bin/bash

set -euo pipefail

echo "========================================="
echo "METHOD 2 - FEELnc CLASSIFICATION"
echo "========================================="


CONFIG_FILE=${1:-../data/config.txt}

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: config file not found: $CONFIG_FILE"
    exit 1
fi


# Output directories
rm -rf results
mkdir -p results
mkdir -p results/logs
mkdir -p results_geneLevel

# FEELnc executable (must be in PATH)
FEELNC_CLASSIFIER=${FEELNC_CLASSIFIER:-FEELnc_classifier.pl}


# Regex definitions
regEx_lncRNA='gene_biotype\s"(lncRNA|lincRNA|sense_intronic|sense_exonic|sense_overlapping|antisense)"'
regEx_mRNA='gene_biotype\s"protein_coding"'


echo "Extracting lncRNA and PCG annotations..."


while IFS=$'\t' read -r shortName completeName ensemblName gtfPath
do

    if [[ "$shortName" == "shortName" ]]; then
        continue
    fi

    if [ ! -f "$gtfPath" ]; then
        echo "WARNING: GTF not found: $gtfPath"
        continue
    fi

    species=$(basename "$gtfPath" .gtf)

    echo "Processing $species"


    grep -v "^#" "$gtfPath" | grep -P "$regEx_lncRNA" > results/${species}_lncRNA.tmp.gtf
    grep -v "^#" "$gtfPath" | grep -P "$regEx_mRNA" > results/${species}_mRNA.tmp.gtf


    echo "Running FEELnc classifier for $species"

    $FEELNC_CLASSIFIER \
        -i results/${species}_lncRNA.tmp.gtf \
        -a results/${species}_mRNA.tmp.gtf \
        -l results/logs/${species}_feelncclassifier.log \
        > results/${species}_classes_feelncclassifier.txt


    echo "Converting transcript level to gene level"

    Rscript 0_FEELnc_tpLevel2gnLevelClassification.R 


done < "$CONFIG_FILE"


echo "Cleaning temporary files..."
rm results/*tmp.gtf


echo "-----------------------------------------"
echo "FEELnc classification completed"
echo "Results in: results/"
echo "-----------------------------------------"