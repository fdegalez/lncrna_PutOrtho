# Header ---------------------------------------------------------------------------------------------
## Aim : Extract orthology relation for couple of species (pariwise) and from BioMart
## Date : 13/03/2023 (last update)
## Author(s): -Fabien DEGALEZ


# Working Directory ---------------------------------------------------------------------------
setwd("~/2_homology_2023/2_extractionOrthology")

# Libraries -----------------------------------------------------------------------------------
library(biomaRt)
library(tidyr)
library(stringi)
library(stringr)

# Variables -----------------------------------------------------------------------------------
config.path <- "../data/config.txt"


# Functions -----------------------------------------------------------------------------------


# Input --------------------------------------------------------------------------------------
config <- read.delim(config.path, header = T, stringsAsFactors = F)

# Script --------------------------------------------------------------------------------------

dir.create("results", showWarnings = FALSE)
## Generation of all combination 
## There is duplicated file but need to keep them (ex: sp1-sp2 vs. sp2-sp1) to deal with the many cases

nameEnsembl.list <- config$ensemblName #TODO: Upgrade the researcg by name for ensembl

# Can't be done with an apply due to API connection
toProcess <- as.data.frame(expand_grid(nameEnsembl.list, nameEnsembl.list))
colnames(toProcess) <- c("source_species", "target_species")
toProcess <- toProcess[toProcess$source_species != toProcess$target_species, ]

# Connection of the bioMart database
# /!\ If a specific version is needed, that must be indicated with "version=" otherwise, the current one is used
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")

# Carefull the ensembl API can reject the connection
# TODO : Register the rejected connection


for (i in 1:nrow(toProcess)){
    name_source <- toProcess$source_species[i]
    name_target <- toProcess$target_species[i]
    cat(i,"/",nrow(toProcess), " - ", name_source, " vs. ", name_target, "\n" , sep = "")
    
    # Connection to the specific species
    ensembl <- useDataset(paste0(name_source, "_gene_ensembl"), mart = ensembl)
    
    homology <- getBM(attributes = c("ensembl_gene_id",
                                     paste0(name_target, "_homolog_ensembl_gene"),
                                     paste0(name_target, "_homolog_orthology_type")),
                      mart = ensembl)
    # Output
    write.table(homology,
                paste0("results/", name_source,"-",name_target, "_homology.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)
}



