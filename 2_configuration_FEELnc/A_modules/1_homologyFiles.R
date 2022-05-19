### Import package
library(biomaRt)

### Import variable
args <- commandArgs(trailingOnly = T)
ensName1 <- args[1]
ensName2 <- args[2]
name1 <- args[3]
name2 <- args[4]

if (ensName1 != ensName2){
##### Connexion to Mart
cat("------------------------------------ \n")
cat(paste("#", name1, "-", name2, "#\n"))
cat("------------------------------------ \n")
cat("Connexion to BioMart ... ")
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset(paste0(ensName1, "_gene_ensembl"), mart = ensembl)
cat("DONE\n")

###### Query building
cat("Creation of the homology file ... ")
res <- getBM(attributes = c("ensembl_gene_id",
                            paste0(ensName2, "_homolog_ensembl_gene"),
                            paste0(ensName2, "_homolog_orthology_type")),
             mart = ensembl) 

colnames(res)[1] <- paste0(ensName1, "_homolog_ensembl_gene")

toName <- paste0(name1, "_", name2, "_homology.tsv")
write.table(res, toName, quote = F, row.names = F, sep = "\t", col.names = T)
cat("DONE\n")
cat("------------------------------------ \n\n")
} else {
  cat("------------------------------------ \n")
  cat(paste("#", name1, "vs.", name2, "#\n"))
  cat("Same Specie, no homology file created !")
  cat("------------------------------------ \n\n")
}
