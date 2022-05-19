library(stringr)
library(stringi)

## Import argument
args <- commandArgs(trailingOnly = T)
name1 <- args[1]
name2 <- args[2]
ensName1 <- args[3]
ensName2 <- args[4]
nammingConfig <- args[5]


## Functions
shorterConf_inter1 <- function(x){
    a <- str_split(x, "_", simplify = T)[1]
    return(a)
}

shorterConf_inter2 <- function(x){
    a <- str_split(x, "_", simplify = T)[1]
    return(str_sub(a, 1, 6))
}

shorterConf_open1 <- function(x){
    a <- str_split(x, "_", simplify = T)[1]
    a <- str_sub(a, 1, 6)
    if (a == "lincSS" | a == "lncgSS") {
        a <- "SS"
    } else if (a == "lincDi" | a == "lncgAS") {
        a <- "AS"
    } else if (a == "lincCo") {
        a <- "Conv"
    } 
}

shorterConf_open2 <- function(x){
    a <- str_split(x, "_", simplify = T)[1]
    a <- str_sub(a, 1, 6)
    if (a == "lincSS" | a == "lncgSS") {
        a <- "SS"
    } else if (a == "lincDi" | a == "lncgAS" | a == "lincCo") {
        a <- "AS"
    }
}

## Import generated files
# Homology file
if ((ensName1 != ensName2)) {
    cat("------------------------------------ \n")
    cat(paste("#", name1, "-", name2, "configuration homology", "#\n"))
    cat("------------------------------------ \n")
    
    ### Treatment of homology file (START)
    # File importation
    homology_file <- list.files(".", pattern = "_homology")
    homology_file <- read.delim(homology_file, header = T, stringsAsFactors = F)
    homology_file <- homology_file[homology_file[,3] == "ortholog_one2one", ]
    ### Treatment of homology file (END)
    
    ### Treatement of FEELnc file (START)
    # lncRNA configuration files
    lnc_1 <- read.delim(paste0(name1, "_lncConfiguration_feelncclassifier.tsv"),
                        header = T, stringsAsFactors = F)
    lnc_2 <- read.delim(paste0(name2, "_lncConfiguration_feelncclassifier.tsv"),
                        header = T, stringsAsFactors = F)
    
    ## lncRNA selection
    lnc_1 <- lnc_1[!is.na(lnc_1$feelLncPcgClassName) 
                   & !grepl("unclassified", lnc_1$feelLncPcgClassName), ]
    
    lnc_2 <- lnc_2[!is.na(lnc_2$feelLncPcgClassName) 
                   & !grepl("unclassified", lnc_2$feelLncPcgClassName), ]
    ### Treatement of FEELnc file (END)
    
    ### Configuration changement - if necessary (START)
    
    cat(paste0("Transformation of configuration in \"", nammingConfig , "\" ... "))
    if (nammingConfig == "inter1") {
        # intermediary 
        lnc_1$feelLncPcgClassName <- sapply(lnc_1$feelLncPcgClassName, shorterConf_inter1)
        lnc_2$feelLncPcgClassName <- sapply(lnc_2$feelLncPcgClassName, shorterConf_inter1)
    } else if (nammingConfig == "inter2") {
        # open - with analogs
        lnc_1$feelLncPcgClassName <- sapply(lnc_1$feelLncPcgClassName, shorterConf_inter2)
        lnc_2$feelLncPcgClassName <- sapply(lnc_2$feelLncPcgClassName, shorterConf_inter2)
    } else if (nammingConfig == "open1") {
        # open - with analogs
        lnc_1$feelLncPcgClassName <- sapply(lnc_1$feelLncPcgClassName, shorterConf_open1)
        lnc_2$feelLncPcgClassName <- sapply(lnc_2$feelLncPcgClassName, shorterConf_open1)
    } else if (nammingConfig == "open2") {
        # open - with analogs
        lnc_1$feelLncPcgClassName <- sapply(lnc_1$feelLncPcgClassName, shorterConf_open2)
        lnc_2$feelLncPcgClassName <- sapply(lnc_2$feelLncPcgClassName, shorterConf_open2)
    }
    cat("DONE\n")
    ### Configuration changement - if necessary (START)
    
    ### Creation of the final file (START)
    colNames <- c(paste0("PCG_", name1),
                  paste0("lncRNA_", name1),
                  paste0("conf_", name1),
                  paste0("class_", name1),
                  paste0("distance_", name1),
                  paste0("PCG_", name2),
                  paste0("lncRNA_", name2),
                  paste0("conf_", name2),
                  paste0("class_", name2),
                  paste0("distance_", name2))
    res <- data.frame(matrix(NA, ncol = length(colNames)),
                      stringsAsFactors = F)
    colnames(res) <- colNames
    res <- res[-1, ] 
    ### Creation of the final file (END)
    
    toExtract <- lnc_1$feelLncPcgGnId %in% homology_file[, 1] | lnc_1$feelLncPcgGnId %in% homology_file[, 2]
    PCG_homologous <- unique(lnc_1[toExtract, "feelLncPcgGnId"])
    
    
    cat("Creation of the configuration FEELNC homology file ... ")
    ## PCG associated
    i = PCG_homologous[1]
    for (i in PCG_homologous) {
        ## Extraction of sub-table for each PCG (START)
        ## La feature fixée est le PCG.
        # Extraction of the lncRNA linked to thse PCG in species 1 
        PCG_lnc_1 <- lnc_1[lnc_1$feelLncPcgGnId == i, ]
        
        # Exaction of the lncRNA linked to the homologous PCG in species 2 
        homologous_i <- unlist(homology_file[(homology_file[, 1] == i | homology_file[, 2] == i), 1:2])
        homologous_i <- setdiff(homologous_i, i)
        if (!(identical(homologous_i, character(0)))) {
            PCG_lnc_2 <- lnc_2[lnc_2$feelLncPcgGnId == homologous_i, ]
        } 
        ## Extraction of sub-table for each PCG (END)
        
        ## Crossed-table (START)
        PCG_lnc_1_sub <- PCG_lnc_1[, c(1,3:5,7)]
        PCG_lnc_2_sub <- PCG_lnc_2[, c(1,3:5,7)]
        
        res_i <- merge(PCG_lnc_1_sub, PCG_lnc_2_sub,
                       by.x = "feelLncPcgClassName", by.y = "feelLncPcgClassName")
        #Column rearragment
        res_i <- res_i[ , c(4,2,1,3,5,8,6,1,7,9)]
        #bind to res
        res <- rbind(res, res_i)
        ## Crossed-table (END)
    }
    colnames(res) <- colNames
    write.table(res, paste0(name1, "_", name2, "_lncConfigurationHomology.tsv"), quote = F, row.names = F, sep = "\t", col.names = T)
    cat("DONE\n")
    cat("------------------------------------ \n\n")
} else {
    cat("------------------------------------ \n")
    cat(paste("#", name1, "vs.", name2, "#\n"))
    cat("Same Specie, no homology file created !\n")
    cat("------------------------------------ \n\n")
}
