#!/usr/bin/env bash

###############################################################################
# Aim: Run the complete Compara / Mercator-Pecan orthology workflow
#
# Steps:
#   1. Format lncRNA coordinates for Compara queries
#   2. Query Ensembl Compara for all species
#   3. Summarize Compara matches across species
#
# Usage:
#   bash run_compara_pipeline.sh [config_file]
#
# Input:
#   ../data/config.txt
#
# Output:
#   input_formatted/
#   results_raw/
#   output_isMatching/
#   + downstream outputs from later parsing/filtering scripts
###############################################################################

set -e
set -o pipefail

CONFIG_FILE=${1:-"../data/config.txt"}

INPUT_FORMAT_SCRIPT="1_inputFormatting.R"
QUERY_SCRIPT="2_run_compara_queries.sh"
SUMMARY_SCRIPT="3_isMatching.R"   

echo "========================================="
echo " METHOD 3 - ENSEMBL COMPARA"
echo "========================================="

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: config file not found: $CONFIG_FILE"
    exit 1
fi

if [[ ! -f "$INPUT_FORMAT_SCRIPT" ]]; then
    echo "ERROR: missing script: $INPUT_FORMAT_SCRIPT"
    exit 1
fi

if [[ ! -f "$QUERY_SCRIPT" ]]; then
    echo "ERROR: missing script: $QUERY_SCRIPT"
    exit 1
fi

if [[ ! -f "$SUMMARY_SCRIPT" ]]; then
    echo "ERROR: missing script: $SUMMARY_SCRIPT"
    exit 1
fi

echo "Using config file: $CONFIG_FILE"
echo "-----------------------------------------"

###############################################################################
# STEP 1 - Format input
###############################################################################

echo ""
echo "STEP 1/3 - Formatting lncRNA coordinates for Compara"
echo "----------------------------------------------------"

Rscript "$INPUT_FORMAT_SCRIPT"

echo "STEP 1 completed"

###############################################################################
# STEP 2 - Run Compara queries
###############################################################################

echo ""
echo "STEP 2/3 - Running Compara queries for all species"
echo "--------------------------------------------------"

bash "$QUERY_SCRIPT" "$CONFIG_FILE"

echo "STEP 2 completed"

###############################################################################
# STEP 3 - Summarize results
###############################################################################

echo ""
echo "STEP 3/3 - Summarizing Compara matches"
echo "--------------------------------------"

Rscript "$SUMMARY_SCRIPT"

echo "STEP 3 completed"


echo ""
echo "========================================="
echo "Compara pipeline completed successfully"
echo "Generated directories:"
echo " - input_formatted/"
echo " - results_raw/"
echo " - results_isMatching/"
echo "========================================="