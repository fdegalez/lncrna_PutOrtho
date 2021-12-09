## Library
library(stringr)
library(stringi)

## Importing data
args <- commandArgs(trailingOnly = T)

source <- args[1]
tableSource <- read.delim(source, header = T, stringsAsFactors = F)
target <- args[2]
tableTarget <- read.delim(target, header = T, stringsAsFactors = F)
homologyInfo <- args[3]
orthology <- read.delim(homologyInfo, header = T, stringsAsFactors = F)
nameSource <- args[4]
nameTarget <- args[5]

# source <- "2_tableGTF_human.tsv"
# tableSource <- read.delim(source, header = T, stringsAsFactors = F)
# target <- "2_tableGTF_mouse.tsv"
# tableTarget <- read.delim(target, header = T, stringsAsFactors = F)
# homologyInfo <- "homologyFile.tsv"
# orthology <- read.delim(homologyInfo, header = T, stringsAsFactors = F)
# nameSource <- "human"
# nameTarget <- "mouse"

cat("------------------------------------ \n")
cat(paste("#", nameSource, " - ", nameTarget, " orthology table creation", "#\n"))
cat("------------------------------------ \n")


## Homology file formatting
orthology <- orthology[orthology[, 3] == "ortholog_one2one", ]

correspondance <- data.frame(tableSource[, 1], tableSource[, 2] %in% orthology[, 1] ,
                             tableSource[, 3] %in% orthology[, 1], stringsAsFactors = F)
colnames(correspondance) <- c("lncRNA_source_id", "PCG_left_source_orthologous", "PCG_right_source_orthologous")
correspondance$PCG_orthologous_nb <- apply(correspondance[, 2:3], 1, sum)

# Number of cases with 0 / 1 / 2 PCG one_to_one
cat("Number of orthologous PCG one_to_one bounding the lncRNA : \n")
cat(paste0(sum(correspondance$PCG_orthologous_nb == 0), "/", nrow(correspondance), " : ", "no ortholog (0)\n"))
cat(paste0(sum(correspondance$PCG_orthologous_nb == 1), "/", nrow(correspondance), " : ", "one ortholog (1) \n"))
cat(paste0(sum(correspondance$PCG_orthologous_nb == 2), "/", nrow(correspondance), " : ", "two orthologs (2) \n"))
cat("------------------------------------ \n")

### Focus on lncRNA with two PCG orthologous one_to_one
lncRNA_concerned <- correspondance[correspondance$PCG_orthologous_nb == 2, "lncRNA_source_id"]
tableSource <- tableSource[tableSource[, 1] %in% lncRNA_concerned, ]

find_lncRNA_orth_byPCGcouple <- function(x){
    
    PCG_left <- x[1]
    PCG_right <- x[2]
    
    table_source_tmp <- tableSource[tableSource[,2] == PCG_left & tableSource[,3] == PCG_right, ] 
    ## Concatenation
    colnames_table_source <- colnames(table_source_tmp)
    table_source_toPaste <- table_source_tmp[, c(1,4,7,8,9)]
    table_source_toPaste <- t(apply(table_source_toPaste, 2, function(x){return(paste0(x, collapse = ";"))}))
    table_source_tmp <- cbind(table_source_tmp[1, c(2,3,5,6)], table_source_toPaste)
    table_source_tmp <- table_source_tmp[, c(5,1,2,6,3,4,7,8,9)]
    colnames(table_source_tmp) <- colnames_table_source
    
    ## Orthologous in the target species
    PCG_orthologous_left <- orthology[match(PCG_left, orthology[,1]), 2]
    PCG_orthologous_right<- orthology[match(PCG_right, orthology[,1]), 2]
    
    table_target_tmp <- tableTarget[((tableTarget[, 2] %in% PCG_orthologous_left) | (tableTarget[, 3] %in% PCG_orthologous_left)) & 
                                        ((tableTarget[, 2] %in% PCG_orthologous_right) | (tableTarget[, 3] %in% PCG_orthologous_right)), ]
    colnames_table_target <- colnames(table_target_tmp)
    if (nrow(table_target_tmp) > 0) {
    ## Concatenation
    table_target_toPaste <- table_target_tmp[, c(1,4,7,8,9)]
    table_target_toPaste <- t(apply(table_target_toPaste, 2, function(x){return(paste0(x, collapse = ";"))}))
    table_target_tmp <- cbind(table_target_tmp[1, c(2,3,5,6)], table_target_toPaste)
    table_target_tmp <- table_target_tmp[, c(5,1,2,6,3,4,7,8,9)]
    colnames(table_target_tmp) <- colnames_table_target
    } else {
        table_target_tmp <- data.frame(t(rep(NA, ncol(table_target_tmp))), stringsAsFactors = F)
        colnames(table_target_tmp) <- colnames_table_target
    }
    
    res <- cbind(table_source_tmp, table_target_tmp)
    res <- sapply(res[1,], as.character)
    return(res)
}

tableSource_PCGcouple <- unique(tableSource[, c(2,3)])
cat("Creation of the table ... ")
orthologous_lncRNA_cases <- data.frame(t(apply(tableSource_PCGcouple, 1, find_lncRNA_orth_byPCGcouple)), stringsAsFactors = F)
cat("DONE\n")
cat("Saving ... ")
write.table(orthologous_lncRNA_cases, paste0("3_", nameSource, "_", nameTarget, "_orthology.tsv"), quote = F, sep = '\t', row.names = F)
cat("DONE\n")
cat("------------------------------------ \n\n")
