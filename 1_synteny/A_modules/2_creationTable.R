### Import libraries
library(stringr)
library(stringi)

## Variables that can be change : 
# RegEx use for identification and selection of lncRNA & PCG
lncRegEx <- "antisense|lincRNA|lncRNA|sense_exonic|sense_intronic"
pcgRegEx <- "protein_coding"

## Importing function to Parse
source("../A_modules/2_Rfunction_gtfAttributeParser.R")

## Importing data
args <- commandArgs(trailingOnly = T)
#file <- "/home/fabien/1_homologyLNCPCG_toPublish_V1/1_synteny/human_comparedTo_mouse/human_geneExtracted.gtf"
#nameSmpl <- "human"
file <- args[1]
nameSmpl <- args[2]


cat("------------------------------------ \n")
cat(paste("#", nameSmpl, " table creation", "#\n"))
cat("------------------------------------ \n")
# Improtation of the GTF containing only genes
GTF <- read.delim(file, header = F, stringsAsFactors = F)
colnames(GTF) <- c("seqname", "source", "feature",
                   "start", "end", "score",
                   "strand", "frame", "attribute")
## Parsing of the attribute field of the GTF
GTF <- data.frame(GTF, gtfAttributeParser(GTF), stringsAsFactors = FALSE)


# Extraction of strand, position to add information later
strandInfo <- GTF[, c("gene_id", "strand")]
posInfo <- GTF[, c("gene_id", "start")]

## Creation of a field containing a simple biotype to identify each features
simpleBiotype <- rep(NA,nrow(GTF))
simpleBiotype[grep(lncRegEx, GTF$gene_biotype)] <- "lnc"
simpleBiotype[grep(pcgRegEx, GTF$gene_biotype)] <- "pcg"
GTF$simpleBiotype <- simpleBiotype
## Selection of PCG and lncRNA only
GTF <- GTF[grep("lnc|pcg", GTF$simpleBiotype), ]

## Extraction of lncRNA id
listLNC <- GTF[ GTF$simpleBiotype == "lnc", ]

# Fonction to find each PCG which bounded the lncRNA
pcgIdEachSide <- function(y){
    x <- y["gene_id"]
    # lncRNA position
    chrLNC <- y["seqname"]
    startLNC <- as.numeric(y["start"])
    endLNC <- as.numeric(y["end"])
    strandLNC <- y["strand"]
    ## Left PCG
    leftPCG <- GTF[GTF$seqname == chrLNC & GTF$start < startLNC & GTF$simpleBiotype == "pcg" , ]
    if (nrow(leftPCG) > 0) {
        leftPCG <- leftPCG[, c("seqname", "start", "end", "strand", "gene_id")]
        leftPCG$diff <- startLNC - leftPCG$start
        leftPCG <- leftPCG[leftPCG$diff == min(leftPCG$diff), ]
        leftPCG <- leftPCG[1, ]
    } else {
        leftPCG <- data.frame(t(rep(NA, 6)), stringsAsFactors = F)
        colnames(leftPCG) <- c("seqname", "start", "end", "strand", "gene_id", "diff")
    }
    
    ## Right PCG
    rightPCG <- GTF[GTF$seqname == chrLNC & GTF$start >= startLNC & GTF$simpleBiotype == "pcg", ]
    if (nrow(rightPCG) > 0) {
        rightPCG <- rightPCG[, c("seqname", "start", "end", "strand", "gene_id")]
        rightPCG$diff <- abs(startLNC - rightPCG$start)
        rightPCG <- rightPCG[rightPCG$diff == min(rightPCG$diff), ]
        rightPCG <- rightPCG[1, ]
    }else {
        rightPCG <- data.frame(t(rep(NA, 6)), stringsAsFactors = F)
        colnames(rightPCG) <- c("seqname", "start", "end", "strand", "gene_id", "diff")
    }
    
    ## Configuration type of the lncRNA compared to the tw PCGs
    if (any(!is.na(leftPCG)) & any(!is.na(rightPCG))) {
        LNC_start <- startLNC
        LNC_end <- endLNC
        leftPCG_start <- leftPCG$start
        leftPCG_end <- leftPCG$end
        rightPCG_start <- rightPCG$start
        rightPCG_end <- rightPCG$end
        
        # Case 1 
        if  (((leftPCG_end >= LNC_start) & (LNC_start >= leftPCG_start)) &
             ((leftPCG_end >= LNC_end) & (LNC_end >= leftPCG_start))){
            LNC_conf <- "contained_PCG_Left"
            # Case 2
        } else if (((leftPCG_end >= LNC_start) & (LNC_start >= leftPCG_start)) &
                   ((leftPCG_end <= LNC_end) & (LNC_end <= rightPCG_start))){
            LNC_conf <- "overlap_PCG_Left"
            # Case 3
        } else if (((leftPCG_end >= LNC_start) & (LNC_start >= leftPCG_start)) &
                   ((rightPCG_end >= LNC_end) & (LNC_end >= rightPCG_start))){
            LNC_conf <- "overlap_PCG_both"
            # Case 4 
        } else if (((rightPCG_start >= LNC_start) & (LNC_start >= leftPCG_end)) &
                   ((rightPCG_start >= LNC_end) & (LNC_end >= leftPCG_end))){
            LNC_conf <- "between"
            # Case 5
        } else if (((rightPCG_start >= LNC_start) & (LNC_start >= leftPCG_end)) &
                   ((rightPCG_start <= LNC_end) & (LNC_end <= rightPCG_end))){
            LNC_conf <- "overlap_PCG_Right"
            # Case 6 
        } else if (((rightPCG_start <= LNC_start) & (LNC_start <= rightPCG_end)) &
                   ((rightPCG_start <= LNC_end) & (LNC_end <= rightPCG_end))){
            LNC_conf <- "contained_PCG_Right"
            # Case 7 
        } else if (((leftPCG_start <= LNC_start) & (LNC_start <= leftPCG_end)) &
                   ((LNC_end >= rightPCG_end))){
            LNC_conf <- "overlap_PCG_left_AND_PCG_right_contained"
            # Case 8
        } else if (((leftPCG_end <= LNC_start) & (LNC_start <= rightPCG_start)) &
                  ((LNC_end >= rightPCG_end))){
            LNC_conf <- "PCG_right_contained"
        } else {
            LNC_conf <- "Error/other"
        }
    } else {
        LNC_conf <- NA
    }
    ### Faire le return
    return(c(x, leftPCG$gene_id, rightPCG$gene_id,
             strandLNC, leftPCG$strand, rightPCG$strand,
             LNC_conf, leftPCG$diff, rightPCG$diff))
}

cat("Creation of the table ... ")
boundedLNC <- data.frame(t(apply(listLNC, 1, pcgIdEachSide)), stringsAsFactors = F)
colnames(boundedLNC) <- c(paste0("lncRNA_", nameSmpl), paste0("PCG_left_", nameSmpl), paste0("PCG_right_", nameSmpl),
                          paste0("lncRNA_strand_", nameSmpl), paste0("PCG_left_strand_", nameSmpl), paste0("PCG_right_strand_", nameSmpl),
                          paste0("lncRNA_conf_", nameSmpl), paste0("lncRNA_PCG_left_ditsance_", nameSmpl), paste0("lncRNA_PCG_right_ditsance_", nameSmpl))
cat("DONE\n")

## Save
cat("Saving ... ")
write.table(boundedLNC, paste0("2_tableGTF_",nameSmpl,".tsv"), quote = F, sep = '\t', row.names = F)
cat("DONE\n")
cat("------------------------------------ \n\n")