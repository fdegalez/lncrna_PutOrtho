#!/usr/bin/env bash

###############################################################################
# Aim: Run the complete lncRNA orthology workflow
#
# Pipeline steps:
#   0. Gene extraction and annotation formatting
#   1. PCG orthology extraction (BioMart)
#   2. Method 1 — Synteny
#   3. Method 2 — FEELnc
#   4. Method 3 — Genome alignment (Compara PECAN)
#
# Usage:
#   bash run_all_methods.sh [options] [config_file]
#
# Options:
#   --all        Run the full pipeline (default)
#   --synteny    Run method 1 only
#   --feelnc     Run method 2 only
#   --compara    Run method 3 only
###############################################################################

set -e
set -o pipefail

###############################################################################
# Default parameters
###############################################################################

ENSEMBL_VERSION=""  # If empty, will use latest version in orthology extraction step

RUN_EXTRACTION=true
RUN_ORTHOLOGY=true
RUN_SYNTENY=false
RUN_FEELNC=false
RUN_COMPARA=false

CONFIG_FILE="data/config.txt"

###############################################################################
# Parse arguments
###############################################################################

while [[ $# -gt 0 ]]; do
    case $1 in
        --synteny)
            RUN_SYNTENY=true
            RUN_EXTRACTION=true
            RUN_ORTHOLOGY=true
            shift
            ;;
        --feelnc)
            RUN_FEELNC=true
            RUN_EXTRACTION=true
            RUN_ORTHOLOGY=true
            shift
            ;;
        --compara)
            RUN_COMPARA=true
            RUN_EXTRACTION=true
            RUN_ORTHOLOGY=true
            shift
            ;;
        --all)
            RUN_SYNTENY=true
            RUN_FEELNC=true
            RUN_COMPARA=true
            shift
            ;;
        *)
            CONFIG_FILE="$1"
            shift
            ;;
    esac
done

###############################################################################
# If no specific method selected → run all
###############################################################################

if ! $RUN_SYNTENY && ! $RUN_FEELNC && ! $RUN_COMPARA; then
    RUN_SYNTENY=true
    RUN_FEELNC=true
    RUN_COMPARA=true
fi

###############################################################################
# Check config file
###############################################################################

echo "=========================================="
echo " LncRNA Orthology Detection Pipeline"
echo "=========================================="

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: config file not found: $CONFIG_FILE"
    exit 1
fi

echo "Using config file: $CONFIG_FILE"
echo "------------------------------------------"

###############################################################################
# STEP 0 — Gene extraction
###############################################################################

if $RUN_EXTRACTION; then

    echo ""
    echo "=========================================="
    echo " STEP 0 — Gene extraction"
    echo "=========================================="

    cd 1_extractionGenes
    bash extract_genes.sh "../$CONFIG_FILE"
    cd ..

    echo "Gene extraction completed"
fi

###############################################################################
# STEP 1 — PCG orthology (BioMart)
###############################################################################

if $RUN_ORTHOLOGY; then

    echo ""
    echo "=========================================="
    echo " STEP 1 — PCG orthology extraction"
    echo "=========================================="

    cd 2_extractionOrthologyPCG
    bash run_OrthologyExtraction.sh "$ENSEMBL_VERSION"
    cd ..

    echo "PCG orthology extraction completed"
fi

###############################################################################
# METHOD 1 — SYNTENY
###############################################################################

if $RUN_SYNTENY; then

    echo ""
    echo "=========================================="
    echo " METHOD 1 — SYNTENY"
    echo "=========================================="

    cd 3_synteny
    bash run_synteny.sh "../$CONFIG_FILE"
    cd ..

    echo "Method 1 completed"
fi

###############################################################################
# METHOD 2 — FEELnc
###############################################################################

if $RUN_FEELNC; then

    echo ""
    echo "=========================================="
    echo " METHOD 2 — FEELnc"
    echo "=========================================="

    cd 4_FEELnc
    bash run_orthoFEELnc.sh "../$CONFIG_FILE"
    cd ..

    echo "Method 2 completed"
fi

###############################################################################
# METHOD 3 — COMPARA
###############################################################################

if $RUN_COMPARA; then

    echo ""
    echo "=========================================="
    echo " METHOD 3 — COMPARA (Mercator / PECAN)"
    echo "=========================================="

    cd 5_compara
    bash run_compara_pipeline.sh "../$CONFIG_FILE"
    cd ..

    echo "Method 3 completed"
fi

###############################################################################
# END
###############################################################################

echo ""
echo "=========================================="
echo " PIPELINE COMPLETED"
echo "=========================================="