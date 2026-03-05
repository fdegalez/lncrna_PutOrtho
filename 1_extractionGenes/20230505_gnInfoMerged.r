# Header ---------------------------------------------------------------------------------------------
## Aim : Create a merged TSV file according to the previous script but adding the origin of the species.
## Date : 13/03/2023 (last update)
## Author(s): - Fabien DEGALEZ

# Working Directory ---------------------------------------------------------------------------
setwd("~/2_homology_2023/1_extractionGenes")

# Libraries -----------------------------------------------------------------------------------
library(stringi)
library(stringr)

geneInfo.listNames <- list.files("./results_gnInfo/", full.names = T, pattern = "gnInfo.tsv")

# Functions -----------------------------------------------------------------------------------
geneInfo.fct <- function(geneInfo.path){
    geneInfo.name <- str_remove(str_split(rev(str_split(geneInfo.path, "/", simplify = T))[1], "\\.", simplify = T)[1], "_gnInfo")
    cat(geneInfo.name, "\n")
    geneInfo <- read.delim(geneInfo.path, header = T, stringsAsFactors = F)
    geneInfo$species <- geneInfo.name
    return(geneInfo)
}

# Input --------------------------------------------------------------------------------------

# Script --------------------------------------------------------------------------------------

geneInfo.list <- NULL
geneInfo.list <- pblapply(geneInfo.listNames, 
                          geneInfo.fct)

geneInfo_merged <- do.call("rbind", geneInfo.list)
write.table(geneInfo_merged, paste0("results_gnInfo/0_allMerged_gnInfo.tsv"), quote = F, sep = '\t', row.names = F)
