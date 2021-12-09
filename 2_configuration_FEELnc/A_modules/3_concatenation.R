## Libraries
library(stringr)
library(stringi)


## Functions
collaspingInfoCol <- function(x){
    return(paste0(x, collapse = ";"))
}

## Improting data

name1<- "human"
name2 <- "mouse"

args <- commandArgs(trailingOnly = T)
name1 <- args[1]
name2 <- args[2]

cat("------------------------------------ \n")
cat(paste("#", name1, "vs.", name2, "concatenation", "#\n"))
cat("------------------------------------ \n")
dta <- read.delim(paste0(name1, "_", name2, "_lncConfigurationHomology.tsv"), header = T, stringsAsFactors = F)

PCG <- unique(dta[, 1])

res <- NULL

#gene <- PCG[2]
#gene <- "ENSG00000142611"
cat("Concatenation of lncRNA in the same ortholog group ... ")
for (gene in PCG) {
    sub <- dta[dta[, 1] == gene, ]
    if (nrow(sub) == 1) {
        res <- rbind(res, sub)
    } else {
        conf <- unique(sub[, 3])
        #conf_i <- conf[2]
        for (conf_i in conf) {
            sub_sub <- sub[sub[, 3] == conf_i, ]
            if (nrow(sub_sub) == 1) {
                res <- rbind(res, sub_sub)
            } else {
                sub_sub_merge <- unique(sub_sub[, c(1,3,6,8)])
                sub_sub_collapse_Left <- data.frame(t(apply(unique(sub_sub[, c(2,4,5)]), 2, collaspingInfoCol)), stringsAsFactors = F)
                sub_sub_collapse_Right <- data.frame(t(apply(unique(sub_sub[, c(7,9,10)]), 2, collaspingInfoCol)), stringsAsFactors = F)
                sub_sub_concat <- cbind(sub_sub_merge, sub_sub_collapse_Left, sub_sub_collapse_Right)
                sub_sub_concat <- sub_sub_concat[, c(1,5,2,6,7,3,8,4,9,10)]
                res <- rbind(res, sub_sub_concat)
            }
        }  
    }
} 
cat("DONE\n")
## Adding information about type
cat("Assignation of the orthology configuration ... ")
whatOrthologies <- function(x){
    nbLncRNA_species1 <- str_count(x[1], ";")+1
    nbLncRNA_species2 <- str_count(x[2], ";")+1
    if (nbLncRNA_species1 == 1 & nbLncRNA_species2 == 1){
        return("one_to_one")
    } else if (nbLncRNA_species1 == 1 & nbLncRNA_species2 > 1) {
        return("one_to_many")
    } else if (nbLncRNA_species1 > 1 & nbLncRNA_species2 == 1) {
        return("many_to_one")
    } else if (nbLncRNA_species1 > 1 & nbLncRNA_species2 > 1) {
        if (nbLncRNA_species1 == nbLncRNA_species2) {
            return("many_to_many_strict")
        } else {
            return("many_to_many")
        }
    }
}

res$orthology_type <- apply(res[, c(2,7)], 1 , whatOrthologies)
cat("DONE\n")
cat("------------------------------------ \n\n")
write.table(res, paste0(name1, "_", name2, "_lncConfigurationHomology_concatenated.tsv"), quote = F, row.names = F, sep = "\t", col.names = T)
