#!/usr/bin/env Rscript

# ==============================================================================
# Script: OrthologyFromBiomart.R
#
# Aim:
#   Extract pairwise PCG orthology relationships from Ensembl BioMart
#
# Usage:
#   Rscript OrthologyFromBiomart.R [ensembl_version]
#
# Example:
#   Rscript OrthologyFromBiomart.R 109
#
# If no version is provided, the latest Ensembl release is used.
#
# Input:
#   ../data/config.txt
#
# Output:
#   results/source-target_homology.tsv
#
# Author: Fabien Degalez
# ==============================================================================


# Libraries --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(biomaRt)
})


# Parameters -------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

ensembl_version <- NA
if(length(args) >= 1){
  ensembl_version <- as.numeric(args[1])
}

config_path <- "../data/config.txt"
output_dir  <- "results"

dir.create(output_dir, showWarnings = FALSE)


# Read config ------------------------------------------------------------------

config <- read.delim(config_path, stringsAsFactors = FALSE)
species <- config$ensemblName


# Generate pairwise combinations -----------------------------------------------

pairs <- expand.grid(
  source_species = species,
  target_species = species,
  stringsAsFactors = FALSE
)

pairs <- pairs[pairs$source_species != pairs$target_species, ]

cat(nrow(pairs), "pairwise comparisons to process\n\n")


# Connect to Ensembl -----------------------------------------------------------

if(is.na(ensembl_version)){
  
  cat("Connecting to latest Ensembl release\n")
  
  mart <- useEnsembl(
    biomart = "genes"
  )
  
} else {
  
  cat("Connecting to Ensembl version:", ensembl_version, "\n")
  
  mart <- useEnsembl(
    biomart = "genes",
    version = ensembl_version
  )
}


# Main loop --------------------------------------------------------------------

for (i in seq_len(nrow(pairs))) {

  source <- pairs$source_species[i]
  target <- pairs$target_species[i]

  cat(i, "/", nrow(pairs), " - ", source, " vs ", target, "\n", sep="")

  dataset <- paste0(source, "_gene_ensembl")

  # Connect to species dataset
  mart_dataset <- tryCatch(
    useDataset(dataset, mart = mart),
    error = function(e) {
      cat("ERROR connecting to dataset:", dataset, "\n")
      return(NULL)
    }
  )

  if (is.null(mart_dataset)) next


  attributes <- c(
    "ensembl_gene_id",
    paste0(target, "_homolog_ensembl_gene"),
    paste0(target, "_homolog_orthology_type")
  )


  # Query BioMart
  homology <- tryCatch(
    getBM(attributes = attributes, mart = mart_dataset),
    error = function(e) {
      cat("BioMart query failed\n")
      return(NULL)
    }
  )


  if (is.null(homology)) next


  # Save results
  out_file <- file.path(
    output_dir,
    paste0(source, "-", target, "_homology.tsv")
  )

  write.table(
    homology,
    out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

}