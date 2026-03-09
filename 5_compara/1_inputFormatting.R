#!/usr/bin/env Rscript

###############################################################################
# Aim: Format lncRNA coordinates to generate input files for Ensembl Compara
#      genome alignment queries.
#
# Input:
#   ../1_extractionGenes/results_gnInfo/*_gnInfo.tsv
#
# Output:
#   input_formatted/lncRNA_formatted_<species>.tsv
#
# Format expected by the Perl Compara script:
#   seqname  start  end  strand  gene_id
#
# Author: Fabien Degalez
###############################################################################

suppressPackageStartupMessages({
  library(stringr)
})

############################
# Parameters
############################

input_dir  <- "../1_extractionGenes/results_gnInfo"
output_dir <- "input_formatted"

lncRNA.regex <- c(
  "lncRNA",
  "lincRNA",
  "antisense",
  "sense_overlapping",
  "sense_intronic"
)

############################
# Create output directory
############################

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

############################
# Input files
############################

file.list <- list.files(
  input_dir,
  full.names = TRUE,
  pattern = "_gnInfo.tsv"
)

# Remove merged files if present
file.list <- file.list[!grepl("allMerged", file.list)]

cat("Number of species detected:", length(file.list), "\n\n")

############################
# Processing
############################

for (file.path in file.list){

  sp.name <- basename(file.path) |>
    str_remove("_gnInfo.tsv")

  cat("Processing:", sp.name, "\n")

  file <- read.delim(
    file.path,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  required_cols <- c("seqname","start","end","strand","gene_id","gene_biotype")

  if(!all(required_cols %in% colnames(file))){
    stop(
      paste("Missing required columns in", file.path)
    )
  }

  # Keep only lncRNA genes
  file <- file[file$gene_biotype %in% lncRNA.regex, ]

  if(nrow(file) == 0){
    warning(paste("No lncRNA found for", sp.name))
    next
  }

  file <- file[, c("seqname","start","end","strand","gene_id")]

  # Convert strand to numeric format expected by Compara
  file$strand[file$strand == "+"] <- 1
  file$strand[file$strand == "-"] <- -1

  write.table(
    file,
    file.path(output_dir, paste0("lncRNA_formatted_", sp.name, ".tsv")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

}

cat("\nFormatting completed.\n")
