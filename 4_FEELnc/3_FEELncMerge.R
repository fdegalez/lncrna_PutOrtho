#!/usr/bin/env Rscript

# Header --------------------------------------------------------------------------------------
## Aim : Merge FEELnc orthology results across all pairwise species comparisons
##       for each reference species.
## Date : 13/03/2023 (last update, revised)
## Author(s): - Fabien DEGALEZ

# Libraries -----------------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(stringr)
    library(pbapply)
})

# Input ---------------------------------------------------------------------------------------

config.path <- "../data/config.txt"
config <- read.delim(config.path, header = TRUE, stringsAsFactors = FALSE)

config.criteria.list <- c("inter1", "inter2", "open1", "open2", "custom")

# Functions -----------------------------------------------------------------------------------

uncollapse_lnc <- function(x){
    
    if(is.na(x[1])){
        return(data.frame(x))
    }
    
    x <- data.frame(t(x))
    x <- x[rep(1, str_count(x[1,1], ";") + 1), ]
    
    x[,1] <- as.character(str_split(x[1,1], ";", simplify = TRUE))
    
    return(x)
}

# Script --------------------------------------------------------------------------------------

for (config.criteria in config.criteria.list){
    
    cat("--", config.criteria, "\n")
    
    outdir <- file.path("results_orthoFeelnc_merged", config.criteria)
    
    unlink(outdir, recursive = TRUE)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    feelnc.list <- list.files(
        file.path("results_orthoFeelnc_agg", config.criteria),
        full.names = TRUE,
        pattern = "_lncConfigurationHomologyAggregated.tsv"
    )
    
    for (i in 1:nrow(config)){
        
        shortName <- config$shortName[i]
        ensemblName <- config$ensemblName[i]
        completeName <- config$completeName[i]
        
        cat(completeName, "\n")
        
        feelnc_interest.list <- feelnc.list[
            grepl(paste0("^", ensemblName, "-"), basename(feelnc.list))
        ]
        
        res <- data.frame(matrix(NA, nrow = 0, ncol = 1), stringsAsFactors = FALSE)
        colnames(res) <- paste0("lncRNA.", completeName)
        
        for (feelnc.path in feelnc_interest.list){
            
            feelnc.file <- read.delim(
                feelnc.path,
                header = TRUE,
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
            
            feelnc.file <- feelnc.file[, c(5,8,11)]
            
            test <- pbapply(feelnc.file, 1, uncollapse_lnc)
            test <- do.call(rbind.data.frame, test)
            
            res <- merge(
                res,
                test,
                by.x = paste0("lncRNA.", completeName),
                by.y = paste0("lncRNA.", completeName),
                all = TRUE
            )
        }
        
        tmp <- res[, grepl("orthology.type", colnames(res)), drop = FALSE]
        
        totalIndication <- function(x){
            sum(!is.na(x))
        }
        
        profilIndication <- function(x){
            
            if (all(is.na(x))){
                return(NA)
            }
            
            x <- x[!is.na(x)]
            
            tab <- data.frame(table(x), stringsAsFactors = FALSE)
            tab <- tab[order(tab$Freq), ]
            
            paste0(apply(tab, 1, paste0, collapse=":"), collapse=";")
        }
        
        nbTotalIndication <- pbapply(tmp, 1, totalIndication)
        profilIndication <- pbapply(tmp, 1, profilIndication)
        
        res$nbTotalIndication <- nbTotalIndication
        res$profilIndication <- profilIndication
        
        write.table(
            res,
            file.path(
                "results_orthoFeelnc_merged",
                config.criteria,
                paste0(ensemblName, "_feelncMerged.tsv")
            ),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE
        )
    }
}
