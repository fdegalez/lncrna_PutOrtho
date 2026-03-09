#!/usr/bin/env Rscript

# Header --------------------------------------------------------------------------------------
## Aim : Infer pairwise putative lncRNA orthologies based on FEELnc configuration
##       using the custom configuration mode to handle genic/intergenic cases due to mis-annotation
## Author(s): - Fabien DEGALEZ

# Libraries -----------------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(stringr)
    library(pbapply)
    library(tidyr)
    library(dplyr)
})

# Variables -----------------------------------------------------------------------------------

config.path <- "../data/config.txt"
config <- read.delim(config.path, header = TRUE, stringsAsFactors = FALSE)
config.criteria <- "custom"
dist_interToGenic <- 5000

# Functions -----------------------------------------------------------------------------------

orthologous.fct <- function(pcg_id) {
    PCG_lnc_1 <- lnc_1[lnc_1$feelLncPcgGnId == pcg_id, , drop = FALSE]
    
    # Extraction of lncRNAs linked to the homologous PCG in species 2
    homologous_i <- unlist(homology.file[(homology.file[, 1] == pcg_id | homology.file[, 2] == pcg_id), 1:2])
    homologous_i <- setdiff(homologous_i, pcg_id)
    
    if (identical(homologous_i, character(0))) {
        return(NA)
    }
    
    PCG_lnc_2 <- lnc_2[lnc_2$feelLncPcgGnId %in% homologous_i, , drop = FALSE]
    
    if (nrow(PCG_lnc_1) == 0 || nrow(PCG_lnc_2) == 0) {
        return(NA)
    }
    
    PCG_lnc_1_sub <- PCG_lnc_1[, c(1, 3:5, 7), drop = FALSE]
    PCG_lnc_2_sub <- PCG_lnc_2[, c(1, 3:5, 7), drop = FALSE]
    
    res_i <- merge(
        PCG_lnc_1_sub,
        PCG_lnc_2_sub,
        by.x = "feelLncPcgClassName",
        by.y = "feelLncPcgClassName"
    )
    
    if (nrow(res_i) == 0) {
        return(NA)
    }
    
    res_i <- res_i[, c(4, 2, 1, 3, 5, 8, 6, 1, 7, 9), drop = FALSE]
    return(res_i)
}

whatOrthologies <- function(x) {
    nbLncRNA_species1 <- str_count(x[1], ";") + 1
    nbLncRNA_species2 <- str_count(x[2], ";") + 1
    
    if (nbLncRNA_species1 == 1 && nbLncRNA_species2 == 1) {
        return("one_to_one")
    } else if (nbLncRNA_species1 == 1 && nbLncRNA_species2 > 1) {
        return("one_to_many")
    } else if (nbLncRNA_species1 > 1 && nbLncRNA_species2 == 1) {
        return("many_to_one")
    } else if (nbLncRNA_species1 > 1 && nbLncRNA_species2 > 1) {
        if (nbLncRNA_species1 == nbLncRNA_species2) {
            return("many_to_many_strict")
        } else {
            return("many_to_many")
        }
    } else {
        return(NA)
    }
}

remove_duplicate <- function(x) {
    split1 <- str_split(x[1], ";", simplify = TRUE)
    toKeep <- !duplicated(as.character(split1))
    
    x[1] <- paste0(str_split(x[1], ";", simplify = TRUE)[toKeep], collapse = ";")
    x[2] <- paste0(str_split(x[2], ";", simplify = TRUE)[toKeep], collapse = ";")
    x[3] <- paste0(str_split(x[3], ";", simplify = TRUE)[toKeep], collapse = ";")
    
    return(c(x[1], x[2], x[3]))
}

custom_config_1 <- function(x, thr_distance) {
    config_lnc <- x["feelLncPcgClassName"]
    
    strand_lnc <- x["strand"]
    end_lnc <- as.numeric(x["end"])
    start_lnc <- as.numeric(x["start"])
    
    strand_pcg <- x["feelLncPcg_strand"]
    end_pcg <- as.numeric(x["feelLncPcg_end"])
    start_pcg <- as.numeric(x["feelLncPcg_start"])
    
    if (config_lnc == "lincDivg") {
        distance <- min(abs(end_lnc - start_pcg), abs(start_lnc - end_pcg))
        if (distance <= thr_distance) {
            return("lncgDivg")
        } else {
            return("lincDivg")
        }
    }
    
    if (config_lnc == "lincConv") {
        distance <- min(abs(end_lnc - start_pcg), abs(start_lnc - end_pcg))
        if (distance <= thr_distance) {
            return("lncgConv")
        } else {
            return("lincConv")
        }
    }
    
    if (grepl("lincSS", config_lnc)) {
        if (strand_lnc == "+") {
            if (end_lnc <= start_pcg) {
                if ((start_pcg - end_lnc) <= thr_distance) {
                    return("lncgSS.up")
                } else {
                    return("lincSS.up")
                }
            } else if (start_lnc >= end_pcg) {
                if ((start_lnc - end_pcg) <= thr_distance) {
                    return("lncgSS.dw")
                } else {
                    return("lincSS.dw")
                }
            }
        } else if (strand_lnc == "-") {
            if (start_lnc > end_pcg) {
                if ((start_lnc - end_pcg) <= thr_distance) {
                    return("lncgSS.up")
                } else {
                    return("lincSS.up")
                }
            } else if (start_pcg >= end_lnc) {
                if ((start_pcg - end_lnc) <= thr_distance) {
                    return("lncgSS.dw")
                } else {
                    return("lincSS.dw")
                }
            }
        }
    }
    
    if (grepl("lncgSS", config_lnc)) {
        if (strand_lnc == "+") {
            if (start_lnc < start_pcg) {
                return("lncgSS.up")
            } else if (start_lnc >= start_pcg) {
                return("lncgSS.dw")
            }
        } else if (strand_lnc == "-") {
            if (start_lnc > start_pcg) {
                return("lncgSS.up")
            } else if (start_lnc <= start_pcg) {
                return("lncgSS.dw")
            }
        }
    }
    
    if (grepl("lncgAS", config_lnc)) {
        if (strand_lnc == "+") {
            if (start_lnc <= end_pcg & start_lnc >= start_pcg) {
                return("lncgDivg")
            } else {
                return("lncgConv")
            }
        } else if (strand_lnc == "-") {
            if (end_lnc <= end_pcg & end_lnc >= start_pcg) {
                return("lncgDivg")
            } else {
                return("lncgConv")
            }
        }
    }
    
    return("ERROR")
}

# Input / Output --------------------------------------------------------------------------------

dir.create("results_orthoFeelnc", showWarnings = FALSE)
dir.create("results_orthoFeelnc_agg", showWarnings = FALSE)

unlink(file.path("results_orthoFeelnc", config.criteria), recursive = TRUE)
dir.create(file.path("results_orthoFeelnc", config.criteria), recursive = TRUE, showWarnings = FALSE)

unlink(file.path("results_orthoFeelnc_agg", config.criteria), recursive = TRUE)
dir.create(file.path("results_orthoFeelnc_agg", config.criteria), recursive = TRUE, showWarnings = FALSE)

toWork <- as.data.frame(expand_grid(config$completeName, config$completeName))
colnames(toWork) <- c("completeName_source", "completeName_target")
toWork$ensemblName_source <- config$ensemblName[match(toWork$completeName_source, config$completeName)]
toWork$ensemblName_target <- config$ensemblName[match(toWork$completeName_target, config$completeName)]
toWork <- toWork[toWork$ensemblName_source != toWork$ensemblName_target, ]

# Script --------------------------------------------------------------------------------------

k <- 1
for (k in 1:nrow(toWork)) {
    
    completeName_source <- as.character(toWork[k, 1])
    completeName_target <- as.character(toWork[k, 2])
    ensemblName_source <- as.character(toWork[k, 3])
    ensemblName_target <- as.character(toWork[k, 4])
    
    cat(ensemblName_source, "-", ensemblName_target, "\n")
    
    ## Import Homology file
    homology.path <- paste0(
        "../2_extractionOrthologyPCG//results/",
        ensemblName_source, "-", ensemblName_target, "_homology.tsv"
    )
    
    homology.file <- read.delim(homology.path, header = TRUE, stringsAsFactors = FALSE)
    homology.file <- homology.file[homology.file[, 3] == "ortholog_one2one", , drop = FALSE]
    
    ## Import Annotation file
    annot_1 <- read.delim(
        paste0("../1_extractionGenes/results_gnInfo/", completeName_source, "_gnInfo.tsv"),
        header = TRUE, stringsAsFactors = FALSE
    )
    
    annot_2 <- read.delim(
        paste0("../1_extractionGenes/results_gnInfo/", completeName_target, "_gnInfo.tsv"),
        header = TRUE, stringsAsFactors = FALSE
    )
    
    ## Import lncRNA configuration files
    lnc_1 <- read.delim(
        paste0("./results_geneLevel/", completeName_source, "_lncConfiguration_feelncclassifier.tsv"),
        header = TRUE, stringsAsFactors = FALSE
    )
    toKeep_lnc1 <- colnames(lnc_1)
    
    lnc_2 <- read.delim(
        paste0("./results_geneLevel/", completeName_target, "_lncConfiguration_feelncclassifier.tsv"),
        header = TRUE, stringsAsFactors = FALSE
    )
    toKeep_lnc2 <- colnames(lnc_2)
    
    ## lncRNA selection
    lnc_1 <- lnc_1[!is.na(lnc_1$feelLncPcgClassName) &
                       !grepl("unclassified", lnc_1$feelLncPcgClassName), , drop = FALSE]
    
    lnc_2 <- lnc_2[!is.na(lnc_2$feelLncPcgClassName) &
                       !grepl("unclassified", lnc_2$feelLncPcgClassName), , drop = FALSE]
    
    ## Add annotation information
    lnc_1 <- merge(lnc_1, annot_1, by.x = "gnId", by.y = "gene_id", all.x = TRUE)
    colnames(annot_1) <- paste0("feelLncPcg", "_", colnames(annot_1))
    lnc_1 <- merge(lnc_1, annot_1, by.x = "feelLncPcgGnId", by.y = "feelLncPcg_gene_id", all.x = TRUE)
    
    lnc_2 <- merge(lnc_2, annot_2, by.x = "gnId", by.y = "gene_id", all.x = TRUE)
    colnames(annot_2) <- paste0("feelLncPcg", "_", colnames(annot_2))
    lnc_2 <- merge(lnc_2, annot_2, by.x = "feelLncPcgGnId", by.y = "feelLncPcg_gene_id", all.x = TRUE)
    
    ## Configuration update before custom reclassification
    lnc_1$feelLncPcgClassName <- str_split(lnc_1$feelLncPcgClassName, "_", simplify = TRUE)[, 1]
    lnc_1$feelLncPcgClassName[lnc_1$feelLncPcgClassName == "lincSS"] <- "lncgSS"
    
    lnc_2$feelLncPcgClassName <- str_split(lnc_2$feelLncPcgClassName, "_", simplify = TRUE)[, 1]
    lnc_2$feelLncPcgClassName[lnc_2$feelLncPcgClassName == "lincSS"] <- "lncgSS"
    
    ## Custom configuration
    lnc_1$feelLncPcgClassName <- apply(lnc_1, 1, custom_config_1, dist_interToGenic, simplify = TRUE)
    lnc_2$feelLncPcgClassName <- apply(lnc_2, 1, custom_config_1, dist_interToGenic, simplify = TRUE)
    
    #print(completeName_source)
    #print(table(lnc_1$feelLncPcgClassName))
    
    ## Return to original columns only
    lnc_1 <- lnc_1[, toKeep_lnc1, drop = FALSE]
    lnc_2 <- lnc_2[, toKeep_lnc2, drop = FALSE]
    
    ### Creation of the final file
    colNames <- c(
        paste0("PCG.", completeName_source),
        paste0("lncRNA.", completeName_source),
        paste0("conf.", completeName_source),
        paste0("class.", completeName_source),
        paste0("distance.", completeName_source),
        paste0("PCG.", completeName_target),
        paste0("lncRNA.", completeName_target),
        paste0("conf.", completeName_target),
        paste0("class.", completeName_target),
        paste0("distance.", completeName_target)
    )
    
    toExtract <- lnc_1$feelLncPcgGnId %in% homology.file[, 1] | lnc_1$feelLncPcgGnId %in% homology.file[, 2]
    PCG_homologous <- unique(lnc_1[toExtract, "feelLncPcgGnId"])
    
    res <- pbsapply(PCG_homologous, orthologous.fct, USE.NAMES = TRUE)
    res <- res[!is.na(res)]
    
    if (length(res) == 0) {
        next
    }
    
    res <- do.call(rbind, res)
    
    colnames(res) <- colNames
    rownames(res) <- 1:nrow(res)
    
    write.table(
        res,
        paste0(
            "./results_orthoFeelnc/", config.criteria, "/",
            ensemblName_source, "-", ensemblName_target,
            "_lncConfigurationHomology.tsv"
        ),
        quote = FALSE, row.names = FALSE, sep = "\t", col.names = TRUE
    )
    
    ## Aggregated version
    res_agg <- aggregate(
        cbind(res[, 2], res[, 4], res[, 5],
              res[, 7], res[, 9], res[, 10]),
        list(res[, 1], res[, 3], res[, 6], res[, 8]),
        paste, collapse = ";"
    )
    
    colnames(res_agg) <- colnames(res[, c(1, 3, 6, 8, 2, 4, 5, 7, 9, 10)])
    
    res_agg[, c(8, 9, 10)] <- data.frame(t(pbapply(res_agg[, c(8, 9, 10)], 1, remove_duplicate)), stringsAsFactors = FALSE)
    res_agg[, c(5, 6, 7)] <- data.frame(t(pbapply(res_agg[, c(5, 6, 7)], 1, remove_duplicate)), stringsAsFactors = FALSE)
    
    res_agg[, paste0("orthology_type.", completeName_source, "-", completeName_target)] <-
        pbapply(res_agg[, c(5, 8)], 1, whatOrthologies)
    
    write.table(
        res_agg,
        paste0(
            "./results_orthoFeelnc_agg/", config.criteria, "/",
            ensemblName_source, "-", ensemblName_target,
            "_lncConfigurationHomologyAggregated.tsv"
        ),
        quote = FALSE, row.names = FALSE, sep = "\t", col.names = TRUE
    )
}
