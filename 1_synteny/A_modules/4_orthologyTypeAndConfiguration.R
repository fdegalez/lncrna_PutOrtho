### Import packages
library(stringr)
library(stringi)

## Importing data
args <- commandArgs(trailingOnly = T)
file <- args[1]
nameSource <- args[2]
nameTarget <- args[3]
tableSource <- read.delim(file, header = T, stringsAsFactors = F)

cat("------------------------------------ \n")
cat(paste("#", nameSource, "-", nameTarget, "orthology type & strand orientation", "#\n"))
cat("------------------------------------ \n")


### Creation of orthology categories
orthology_categories <- function(x){
    nb_lncRNA_source <- str_count(x[1], ";") + 1
    nb_lncRNA_target <- str_count(x[10], ";") + 1
    
    if ((nb_lncRNA_source == 1) & (is.na(nb_lncRNA_target))){
        return("one_to_zero")
    } else if ((nb_lncRNA_source > 1) & (is.na(nb_lncRNA_target))){
        return("many_to_zero")
    } else if ((nb_lncRNA_source == 1) & (nb_lncRNA_target == 1)){
        return("one_to_one")
    } else if ((nb_lncRNA_source == 1) & (nb_lncRNA_target > 1)){
        return("one_to_many")
    } else if ((nb_lncRNA_source > 1) & (nb_lncRNA_target == 1)){
        return("many_to_one")
    } else if ((nb_lncRNA_source > 1) & (nb_lncRNA_target > 1)){
        return("many_to_many")
    } else {
        return("Error")
    }
}

cat("Creation of the orthology type ... ")

tableSource$orthology_type <- apply(tableSource, 1, orthology_categories)
cat("DONE\n")

## Creation of orientation conservation
## Table of possible configuration
strand_case_possibility <- data.frame(c("+++",
                                        "++-", "+-+", "-++",
                                        "+--", "-+-", "--+",
                                        "---"),
                                      c(1,
                                        2, 3, 4,
                                        4, 3, 2,
                                        1), 
                                      stringsAsFactors = F)
colnames(strand_case_possibility) <- c("combinaison", "group")


## Collapsing of the "many" case
strand_collapsing <- function(x){
    strand <- as.character(str_split(x, ";", simplify = T))
    strand_res <- c(strand[1])
    initialStrand <- strand[1]
    for (i in strand) {
        if (i != initialStrand)
            strand_res <- c(strand_res, i)
        initialStrand <- i 
    }
    return(paste0(strand_res, collapse = ";"))
}

## Orienation and group following the strand_case_possible table
strand_orientation <- function(x){
    ## Extracxtion of orientation from source species 
    strand_lncRNA_source <- strand_collapsing(x[4])
    strand_source <- paste0(c(x[5], strand_lncRNA_source, x[6] ), collapse = "")
    ## Extracxtion of orientation from target species 
    strand_lncRNA_target <- strand_collapsing(x[13])
    strand_target <- paste0(c(x[14], strand_lncRNA_target, x[15] ), collapse = "")
    ## Look for a matching group
    resSource <- strand_case_possibility[match(strand_source, strand_case_possibility$combinaison),]
    resTarget <- strand_case_possibility[match(strand_target, strand_case_possibility$combinaison),]
    
    if (nchar(strand_source) != nchar(strand_target)){
        orientation_group <- NA
        orientation_class <- "discordant_multi"
    } else {
        if ((nchar(strand_source) > 3) &  (nchar(strand_target) > 3)){
            if (strand_source == strand_target){
                orientation_group <- NA
                orientation_class <- "same_multi"
            } else {
                strand_source <- str_replace_all(strand_source, "\\+", "a")
                strand_source <- str_replace_all(strand_source, "\\-", "\\+")
                strand_source <- str_replace_all(strand_source, "a", "\\-")
                if (strand_source == strand_target){
                    orientation_group <- NA
                    orientation_class <- "reverse_multi"
                } else {
                    orientation_group <- NA
                    orientation_class <- "discordant_multi"
                }
            }
        } else if  ((nchar(strand_source) == 3) &  (nchar(strand_target) == 3)){
            if (resSource$combinaison == resTarget$combinaison){
                orientation_group <- resSource$group
                orientation_class <- "same"
            } else if ((resSource$combinaison != resTarget$combinaison) & (resSource$group == resTarget$group)){
                orientation_group <- resSource$group
                orientation_class <- "reverse"
            } else {
                orientation_group <- NA
                orientation_class <- "discordant"
            }
        }
    }
    return(c(orientation_class, orientation_group))
}

cat("Creation of the strand conservation ... ")

## Split into two tables to deal with the case where at least one orthologuous lncRNA can be found
tableSource_orientation_strand <- tableSource[(tableSource$orthology_type != "one_to_zero") & (tableSource$orthology_type != "many_to_zero"), ]
tableSource_no_orientation_strand <- tableSource[(tableSource$orthology_type == "one_to_zero") | (tableSource$orthology_type == "many_to_zero"), ]
cat("DONE\n")

## Addition of the class and group orientation
toAdd <- data.frame(t(apply(tableSource_orientation_strand, 1, strand_orientation)), stringsAsFactors = F)
colnames(toAdd) <- c("orientation_class", "orientation_group")
tableSource_orientation_strand <- cbind(tableSource_orientation_strand, toAdd)
tableSource_no_orientation_strand$orientation_class <- NA
tableSource_no_orientation_strand$orientation_group <- NA

## Merginf of the two sub-table
tableSource <- rbind(tableSource_orientation_strand, tableSource_no_orientation_strand)

cat("Saving ... ")
write.table(tableSource, paste0("4_", nameSource , "_", nameTarget, "_orthology_strandConservation.tsv"), quote = F, sep = '\t', row.names = F)
cat("DONE\n")
cat("------------------------------------ \n\n")