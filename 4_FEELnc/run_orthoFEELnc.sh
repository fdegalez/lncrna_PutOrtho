#!/usr/bin/env bash

###############################################################################
# Aim: Run the complete FEELnc-based orthology inference workflow
#
# Steps:
#   1. Extract lncRNA/mRNA and run FEELnc classifier & Convert transcript-level FEELnc annotation to gene-level
#   2. Infer orthology relationships based on FEELnc configuration
#   3. Merge orthology results by reference species
#
# Author: Fabien Degalez
###############################################################################

set -e
set -o pipefail

############################
# Arguments
############################

CONFIG_FILE=${1:-"../data/config.txt"}

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: config file not found: $CONFIG_FILE"
    exit 1
fi

echo "Using config file: $CONFIG_FILE"
echo "------------------------------------"


############################
# STEP 1 — FEELnc extraction
############################

echo ""
echo "STEP 1 — Running FEELnc classification"
echo "------------------------------------"

bash 1_extractionFEELnc.sh "$CONFIG_FILE"

echo "STEP 1 completed"
echo ""


############################
# STEP 2 — Orthology inference
############################

echo ""
echo "STEP 2 — Inferring orthology from FEELnc configurations"
echo "------------------------------------"

Rscript 2_orthoFeelncBySpecies.R

echo "STEP 2 completed"
echo ""


############################
# STEP 3 — Merge results by species
############################

echo ""
echo "STEP 3 — Merging orthology results by species"
echo "------------------------------------"

Rscript 3_FEELncMerge.R

echo "STEP 3 completed"
echo ""


echo "------------------------------------"
echo "FEELnc orthology pipeline completed successfully"
echo "------------------------------------"