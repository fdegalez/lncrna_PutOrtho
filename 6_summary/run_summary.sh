#!/usr/bin/env bash

###############################################################################
# Run summary analyses
# 1) Pairwise summary between species
# 2) Species-centered summary (one species vs all others)
###############################################################################

set -e

echo "----------------------------------------"
echo "Running summary module"
echo "----------------------------------------"

PAIRWISE_SCRIPT="1_summaryPairwise.R"
ALL_SCRIPT="2_summaryAll.R"

# Check scripts exist
if [ ! -f "$PAIRWISE_SCRIPT" ]; then
    echo "ERROR: $PAIRWISE_SCRIPT not found"
    exit 1
fi

if [ ! -f "$ALL_SCRIPT" ]; then
    echo "ERROR: $ALL_SCRIPT not found"
    exit 1
fi

echo ""
echo "Step 1: Pairwise summary"
echo "----------------------------------------"
Rscript "$PAIRWISE_SCRIPT"

echo ""
echo "Step 2: Species-level summary"
echo "----------------------------------------"
Rscript "$ALL_SCRIPT"

echo ""
echo "----------------------------------------"
echo "Summary module completed"
echo "----------------------------------------"