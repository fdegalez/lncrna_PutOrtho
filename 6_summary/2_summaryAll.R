#!/usr/bin/env Rscript

###############################################################################
# Summary per reference species
# Merge all pairwise summaries for a given species
###############################################################################

suppressPackageStartupMessages({
    library(stringr)
})

###############################################################################
# Load configuration
###############################################################################

config <- read.delim(
    "../data/config.txt",
    header = TRUE,
    stringsAsFactors = FALSE
)

sp_list <- config$completeName

dir.create("all", showWarnings = FALSE)

###############################################################################
# Main loop
###############################################################################


sp <- sp_list[1]
for (sp in sp_list){
    
    cat("Processing species:", sp, "\n")
    
    res <- NULL
    
    sp_toMerge <- sp_list[sp_list != sp]
    
    
    sp_toAdd <- sp_toMerge[1]
    for (sp_toAdd in sp_toMerge){
        
        cat("\t-", sp_toAdd, "\n")
        
        summary_path <- paste0(
            "pairwise//",
            sp, "_", sp_toAdd,
            "_summary.tsv"
        )
        
        if(!file.exists(summary_path)){
            warning("Missing summary file: ", summary_path)
            next
        }
        
        tmp <- read.delim(
            summary_path,
            header = TRUE,
            stringsAsFactors = FALSE
        )
        
        ########################################
        # Remove "_to_zero" cases
        ########################################
        
        col_ortho <- grep("orthology_type", colnames(tmp))
        
        if(length(col_ortho) > 0){
            idx <- grepl("_to_zero", tmp[,col_ortho])
            tmp[idx, grep("isTable_2", colnames(tmp))] <- 0
        }
        
        ########################################
        # Merge tables
        ########################################
        
        if(is.null(res)){
            
            res <- tmp
            
        } else {
            
            res <- merge(
                res,
                tmp,
                by = 1,
                all = TRUE
            )
            
        }
        
    }
    
    ########################################
    # Compute summary statistics
    ########################################
    
    for (numCol in grep("isTable", colnames(res))){
        res[is.na(res[, numCol]), numCol] <- 0 
    }
    
    res$sumIsTable <- apply(
        res[,grep("isTable_", colnames(res)), drop = FALSE],
        1,
        sum,
        na.rm = TRUE
    )
    
    res$sumIsTable_1 <- apply(
        res[,grep("isTable_1", colnames(res)), drop = FALSE],
        1,
        sum,
        na.rm = TRUE
    )
    
    res$sumIsTable_2 <- apply(
        res[,grep("isTable_2", colnames(res)), drop = FALSE],
        1,
        sum,
        na.rm = TRUE
    )
    
    res$sumIsTable_3 <- apply(
        res[,grep("isTable_3", colnames(res)), drop = FALSE],
        1,
        sum,
        na.rm = TRUE
    )
    
    ########################################
    # Write output
    ########################################
    
    write.table(
        res,
        paste0(
            "all//summary_",
            sp,
            "_all.tsv"
        ),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )
    
}

