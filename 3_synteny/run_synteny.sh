#!/bin/bash

set -euo pipefail

echo "======================================="
echo "   METHOD 1 - SYNTENY ORTHOLOGY"
echo "======================================="

CONFIG=${1:-"../data/config.txt"}

# Check config
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: config file not found: $CONFIG"
    exit 1
fi


# Check required directories
if [ ! -d "../1_extractionGenes/results" ]; then
    echo "ERROR: missing ../1_extractionGenes/results"
    echo "Run step 1 first."
    exit 1
fi

if [ ! -d "../2_extractionOrthologyPCG/results" ]; then
    echo "ERROR: missing ../2_extractionOrthologyPCG/results"
    echo "Run step 2 first."
    exit 1
fi


# Create output directories
mkdir -p results_table
mkdir -p results_synteny
mkdir -p results_syntenyMerged


echo ""
echo "STEP 1/3 : Creating lncRNA-PCG neighborhood tables"
echo "---------------------------------------------------"

Rscript 1_creationTableLncRNAbetweenPCG.R


echo ""
echo "STEP 2/3 : Inferring pairwise lncRNA orthology by synteny"
echo "---------------------------------------------------------"

Rscript 2_syntenyBySpecies.R


echo ""
echo "STEP 3/3 : Merging orthology evidence across species"
echo "----------------------------------------------------"

Rscript 3_syntenyMerge.R


echo ""
echo "======================================="
echo "SYNTENY METHOD COMPLETED"
echo "Results:"
echo " - results_table/"
echo " - results_synteny/"
echo " - results_syntenyMerged/"
echo "======================================="
