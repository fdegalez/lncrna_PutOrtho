#!/usr/bin/env Rscript

# ==============================================================================
# Script: extract_geneInfo.R
# Aim:
#   1. Extract gene information (gene_id, gene_name, gene_biotype) from GTF files
#   2. Produce one TSV per species
#   3. Produce a merged TSV across all species
#
# Input:
#   results/*genesOnly.gtf
#
# Output:
#   results_gnInfo/*_gnInfo.tsv
#   results_gnInfo/0_allMerged_gnInfo.tsv
#
# Author: Fabien Degalez
# ==============================================================================


# Libraries --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(stringr)
})


# Input directory --------------------------------------------------------------

input_dir  <- "results"
output_dir <- "results_gnInfo"

dir.create(output_dir, showWarnings = FALSE)


# List GTF files ---------------------------------------------------------------

gtf_files <- list.files(
  input_dir,
  pattern = "genesOnly.gtf$",
  full.names = TRUE
)

cat(length(gtf_files), "GTF files detected\n\n")


# Function to extract attribute -----------------------------------------------

extract_attr <- function(attr, field) {
  pattern <- paste0("(^|;\\s*)", field, "\\s+([^;]+)")
  m <- regexec(pattern, attr)
  res <- regmatches(attr, m)
  out <- sapply(res, function(x) {
    if (length(x) >= 3) trimws(x[3]) else NA_character_
  })
  return(out)
}


# Storage for merged results ---------------------------------------------------

merged_list <- list()


# Main loop --------------------------------------------------------------------

for (gtf in gtf_files) {

  species <- str_split(str_remove(basename(gtf), "_genesOnly.gtf"), "\\.", simplify = T)[,1]

  cat("Processing:", species, "\n")

  gtf_df <- read.delim(
    gtf,
    header = FALSE,
    stringsAsFactors = FALSE,
    sep = "\t"
  )

  colnames(gtf_df) <- c(
    "seqname","source","feature",
    "start","end","score",
    "strand","frame","attribute"
  )


  # Extract attributes ---------------------------------------------------------

  gene_id      <- extract_attr(gtf_df$attribute, "gene_id")
  gene_name    <- extract_attr(gtf_df$attribute, "gene_name")
  gene_biotype <- extract_attr(gtf_df$attribute, "gene_biotype")


  gene_df <- data.frame(
    gene_id      = gene_id,
    gene_name    = gene_name,
    gene_biotype = gene_biotype,
    seqname      = gtf_df$seqname,
    start        = gtf_df$start,
    end          = gtf_df$end,
    strand       = gtf_df$strand,
    species      = species,
    stringsAsFactors = FALSE
  )


  # Save per-species file ------------------------------------------------------

  out_file <- file.path(
    output_dir,
    paste0(species, "_gnInfo.tsv")
  )

  write.table(
    gene_df,
    out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )


  merged_list[[species]] <- gene_df

}


# Merge all species ------------------------------------------------------------

cat("\nMerging all species...\n")

merged_df <- do.call(rbind, merged_list)

write.table(
  merged_df,
  file.path(output_dir, "0_allMerged_gnInfo.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done.\n")
