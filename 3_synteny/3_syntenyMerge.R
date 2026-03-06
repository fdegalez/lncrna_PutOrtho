#!/usr/bin/env Rscript

# ==============================================================================
# Script: syntenyMerge.R
# Aim:
#   Merge pairwise synteny results for each reference species
# ==============================================================================


# Libraries --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(stringr)
  library(pbapply)
})


# Functions --------------------------------------------------------------------

uncollapse_lnc <- function(x){

  n <- str_count(x[1], ";") + 1

  res <- data.frame(matrix(rep(x, n), nrow = n, byrow = TRUE),
                    stringsAsFactors = FALSE)

  res[,1] <- str_split(x[1], ";")[[1]]

  return(res)

}


totalIndication <- function(x){

  sum(!is.na(x))

}


profilIndication_fun <- function(x){

  if(all(is.na(x))) return(NA)

  x <- x[!is.na(x)]

  tab <- data.frame(table(x), stringsAsFactors = FALSE)

  tab <- tab[order(tab$Freq),]

  paste0(
    apply(tab,1,paste0,collapse=":"),
    collapse=";"
  )

}


# Input ------------------------------------------------------------------------

config <- read.delim("../data/config.txt", stringsAsFactors = FALSE)

synteny_files <- list.files(
  "results_synteny",
  pattern="synteny.tsv",
  full.names=TRUE
)

dir.create("results_syntenyMerged", showWarnings = FALSE)


# Main loop --------------------------------------------------------------------

for(i in seq_len(nrow(config))){

  shortName    <- config$shortName[i]
  ensemblName  <- config$ensemblName[i]
  completeName <- config$completeName[i]

  cat("Processing:",completeName,"\n")


  # Select relevant synteny files

  synteny_interest <- synteny_files[
    grep(paste0(ensemblName,"-"), synteny_files)
  ]


  res <- data.frame(
    matrix(NA,nrow=0,ncol=1),
    stringsAsFactors=FALSE
  )

  colnames(res) <- paste0("lncRNA.",completeName)


  # Merge pairwise results -----------------------------------------------------

  for(path in synteny_interest){

    synteny_file <- read.delim(path, stringsAsFactors=FALSE)

    synteny_file <- synteny_file[,c(1,10,20)]

    expanded <- pbapply(synteny_file,1,uncollapse_lnc)

    expanded <- do.call(rbind,expanded)

    res <- merge(
      res,
      expanded,
      by.x=paste0("lncRNA.",completeName),
      by.y=paste0("lncRNA.",completeName),
      all=TRUE
    )

  }


  # Clean categories -----------------------------------------------------------

  res[res=="one_to_zero"]  <- NA
  res[res=="many_to_zero"] <- NA


  tmp <- res[,grepl("orthology.type",colnames(res))]


  # Compute summary metrics ----------------------------------------------------

  res$nbTotalIndication <- pbapply(tmp,1,totalIndication)

  res$profilIndication  <- pbapply(tmp,1,profilIndication_fun)


  # Output ---------------------------------------------------------------------

  write.table(
    res,
    paste0(
      "results_syntenyMerged/",
      ensemblName,
      "_syntenyMerged.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

}