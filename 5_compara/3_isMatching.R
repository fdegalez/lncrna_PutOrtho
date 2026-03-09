#!/usr/bin/env Rscript

###############################################################################
# Aim: Summarize Compara genomic alignment matches for each lncRNA
#
# Input:
#   output/MP_<species>.tsv
#   equivalenceNameEnsembl.xlsx
#
# Output:
#   output_isMatching/<species>_isMatching.tsv
###############################################################################

suppressPackageStartupMessages({
    library(readxl)
    library(stringr)
})

############################
# Input files
############################

nameEq <- read_excel("../0_metafiles/equivalenceName_ensembl.xlsx")

config.path <- "../data/config.txt"
config <- read.delim(config.path, header = TRUE, stringsAsFactors = FALSE)

############################
# Output directory
############################

dir.create("results_isMatching", showWarnings = FALSE)

############################
# Input Compara results
############################

MP.list <- list.files(
    "results_raw",
    full.names = TRUE,
)

############################
# Processing
############################

MP.file <- MP.list[1]
for (MP.file in MP.list){
    
    cat("Processing:", MP.file, "\n")
    
    dta <- read.delim(MP.file, header = FALSE, stringsAsFactors = FALSE)
    
    nameSp <- basename(MP.file)
    nameSp <- str_remove(nameSp, "MP_")
    nameSp <- str_remove(nameSp, "_compara.tsv")
    
    colnames(dta) <- c(
        "speciesName_source",
        "seq_id",
        "seq_region_start",
        "seq_region_end",
        "seq_region_strand",
        "species_scientificName",
        "species_shortName",
        "species_displayName",
        "genebuild",
        "assembly",
        "seqRegionName",
        "start",
        "end",
        "strand",
        "genes_inSlice"
    )
    
    ########################################
    # Remove self matches
    ########################################
    
    nameSp_short <- nameEq$species_shortName[
        which(nameSp == nameEq$species_customName)
    ]
    
    dta <- dta[dta$species_shortName != nameSp_short, ]
    
    ########################################
    # Count species matches per lncRNA
    ########################################
    
    nbSpeciesMatch_tot <- function(x){
        
        species <- unique(x$species_displayName)
        
        c(
            length(species),
            paste(sort(species), collapse=";")
        )
    }
    
    split_dta <- split(dta, dta$seq_id)
    
    tmp <- lapply(split_dta, nbSpeciesMatch_tot)
    
    df <- data.frame(
        lncRNA = names(tmp),
        nbSpeciesMatch_tot = sapply(tmp, "[", 1),
        SpeciesMatch_tot = sapply(tmp, "[", 2),
        stringsAsFactors = FALSE
    )
    
    df$nbSpeciesMatch_tot <- as.numeric(df$nbSpeciesMatch_tot)
    
    ########################################
    # Species set filtering
    ########################################
    
    toTest <- nameEq[nameEq$inStudiedSet == 1, ]
    toTest <- toTest[toTest$species_shortName != nameSp_short, ]
    
    for (name in toTest$species_displayName){
        
        toName <- nameEq$species_customName[
            match(name, nameEq$species_displayName)
        ]
        
        df[[paste0(toName, "_isMatching")]] <-
            as.numeric(grepl(name, df$SpeciesMatch_tot))
        
    }
    
    ########################################
    # Number of matches in the studied set
    ########################################
    
    df$nbSpeciesMatch_set <-
        apply(df[,4:ncol(df)], 1, sum)
    
    ########################################
    # Output
    ########################################
    
    write.table(
        df,
        paste0("results_isMatching/", nameSp, "_isMatching.tsv"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )
    
}
