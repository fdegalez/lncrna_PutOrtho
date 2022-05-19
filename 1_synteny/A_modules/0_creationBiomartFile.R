### Import package
library(biomaRt)

### Test Variables
#nameSpeciesSource <- "hsapiens"
#nameSpeciesTarget <- "mmusculus"
#ShortNameSpeciesSource <- "human"
#ShortNameSpeciesTarget <- "mouse"

### Import variables
args <- commandArgs(trailingOnly = T)
nameSpeciesSource <- args[1]
nameSpeciesTarget <- args[2]
ShortNameSpeciesSource <- args[3]
ShortNameSpeciesTarget <- args[4]

##### Connexion to Mart
cat("------------------------------------ \n")
cat(paste("#", ShortNameSpeciesSource, " - ", ShortNameSpeciesTarget, "#\n"))
cat("------------------------------------ \n")
cat("Connexion to BioMart ... ")
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
#datasets <- listDatasets(ensembl)
#head(datasets)
ensembl <- useDataset(paste0(nameSpeciesSource, "_gene_ensembl"), mart = ensembl)
cat("DONE\n")

###### Query building
cat("Creation of the homology file ... ")
res <- getBM(attributes = c("ensembl_gene_id",
                          paste0(nameSpeciesTarget, "_homolog_ensembl_gene"),
                          paste0(nameSpeciesTarget, "_homolog_orthology_type")),
             mart = ensembl) 

write.table(res, "homologyFile.tsv", quote = F, row.names = F, sep = "\t", col.names = T)
cat("DONE\n")
cat("------------------------------------ \n\n")