# Header ---------------------------------------------------------------------------------------------
## Aim : Create a TSV file with standardized column and information from the GTF files provided.
## Date : 13/03/2023 (last update)
## Author(s): - Fabien DEGALEZ


# Working Directory ---------------------------------------------------------------------------
setwd("~/2_homology_2023/1_extractionGenes")

# Libraries -----------------------------------------------------------------------------------
library(stringi)
library(stringr)


# Variables -----------------------------------------------------------------------------------


# Functions -----------------------------------------------------------------------------------
extract_field <- function(x, fieldName) {
    # Extract the fieldName (ex : gene_id) from the "attributes" column of a gtf
    regEx <- paste0(fieldName, "[^;]*")
    tmp  <- str_split(str_extract(x, regEx), " ", simplify = T)
    if (all(is.na(tmp[2:length(tmp)] == ""))) {
        return("")
    } else {
        tmp <- paste0(tmp[2:length(tmp)], collapse = " ")
        tmp <- gsub('"', '', tmp)
        return(tmp)
    }
} 

# Input --------------------------------------------------------------------------------------

GTF.list <- list.files("./results", full.names = T, pattern = "genesOnly.gtf")

# Script --------------------------------------------------------------------------------------

dir.create("results_gnInfo", showWarnings = FALSE)
for (GTF.path in GTF.list){
    # Extract the names from the GTF
    GTF.name <- str_split(rev(str_split(GTF.path, "/", simplify = T))[1], "\\.", simplify = T)[1]
    cat(GTF.name, "\n")
    GTF <- read.delim(GTF.path, header = F, stringsAsFactors = F)
    colnames(GTF) <- c("seqname", "source", "feature",
                       "start", "end", "score",
                       "strand", "frame", "attribute")
    ## Parsing of the attribute field of the GTF
    GTF$gene_id <- pbsapply(GTF$attribute, extract_field, "gene_id")
    GTF$gene_name <- pbsapply(GTF$attribute, extract_field, "gene_name")
    GTF$gene_biotype <- pbsapply(GTF$attribute, extract_field, "gene_biotype")
    GTF <- GTF[,c("gene_id", "gene_name", "gene_biotype",
                  "seqname",
                  "start", "end",
                  "strand")]
    GTF[GTF == ""] <- NA
    # Output
    write.table(GTF, paste0("results_gnInfo/",GTF.name,"_gnInfo.tsv"), quote = F, sep = '\t', row.names = F)
}



# Output --------------------------------------------------------------------------------------


