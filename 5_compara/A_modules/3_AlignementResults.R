#setwd("~/1_homologyLNCPCG/mappingAddition/")

### Variables (start)
alignmentSource <- "human"
alignmentOn <- "mouse"
sizeBlockAuthorized <- 500 # in bp
## Import argument
args <- commandArgs(trailingOnly = T)
alignmentSource <- args[1]
alignmentOn <- args[2]
sizeBlockAuthorized <- args[3] # in bp

regEx_LNC <- "lncRNA|lincRNA|antisense|sense_intronic|sense_exonic|sense_overlapping"
### Variables (end)

### Loading packages (start)
library(stringi)
library(stringr)
### Loading packages (end)

### Loading functions (start)
# GTF attributeParser
source("../../A_modules/3_Rfunction_gtfAttributeParser.R") #ToChange

cat("------------------------------------ \n")
cat(paste("#", alignmentSource, "vs.", alignmentOn, " matching comparison", "#\n"))
cat("------------------------------------ \n")

# Extraction of genetic information from string as "chr:start-end:strand"
extractionGenomic <- function(x){
    # ex 1:16324642-16325166:1
    temp <- str_split(x[2], ":", simplify = T)
    temp2 <- str_split(temp[2], "-", simplify = T)
    chr <- temp[1]
    strand <- temp[3]
    start <- as.numeric(temp2[1])
    end <- as.numeric(temp2[2])
    if (strand == 1) {strand = "+"
    } else if (strand == -1) { strand = "-"}
    # ex : c( ,1 ,16324642; 16325166, +)
    return(c(x[1], chr, start, end , strand))
}
### (end)

### Necessary files (start)
# Short name no longe available in Ensembl/Compara API
# Use of a equivalence file for names concerning the Amniotes groups, used in the analysis.
cat("Loading of the equivalence file ... ")
equivalenceName <- read.delim("../../A_modules/3_nameEquivalenceEnsemblAmniotes.txt", header = F,
                              sep= "\t", stringsAsFactors = F, strip.white = T)
equivalenceName <- equivalenceName[, 1:2]
colnames(equivalenceName) <- c("shortName", "longName")
equivalenceName$shortName <- tolower(equivalenceName$shortName)
cat("DONE\n")
# GTF of the "alignementOn" species to extract biotype and gnId
cat("Reading the GTF target file ... ")
gtf <- read.delim(paste0(alignmentOn, "_geneExtracted.gtf"), header = F, stringsAsFactors = T)
colnames(gtf) <- c("seqname", "source", "feature",
                   "start", "end", "score",
                   "strand", "frame", "attribute")

gtf <- data.frame(gtf, gtfAttributeParser(gtf), stringsAsFactors = FALSE)
cat("DONE\n")
# List of all files (1 by lncRNA mapped) from the "alignementSource" species
allFiles <- list.files(paste0("./", alignmentSource, "_alignement_MP_63amniotes/"), pattern = "_alignementAmniotes.txt")
### Necessary files (end)

### Beginning of the loop (start)
# Creation of the results file
resultsMatching <- data.frame((matrix(NA, ncol = 4, nrow = 1)), stringsAsFactors = F)
colnames(resultsMatching) <- c("gnId_source", "isMatching", "gnId_target", "gnBiotype_target")

# Two blocks but near
#file <- "ENSG00000223479_human_alignementAmniotes.txt"
# Two blocks not on the same chr
#file <- "ENSG00000237568_human_alignementAmniotes.txt"
cat("Comparison of the alignement zone of all the lncRNAs from the source species ... ")
for (file in allFiles){
    
    # Opening of the file/lncRNA match
    dta_file <- read.delim(paste0("./", alignmentSource, "_alignement_MP_63amniotes/",file), header = F, stringsAsFactors = F, sep = "\t")
    colnames(dta_file) <- c("speciesLongName", "assembly", "genCoord")
    
    # Extraction of the information of the "alignmentSource" species
    infoSpeciesSource <- dta_file[1,]
    strand_SpeciesSource <- unlist(extractionGenomic(infoSpeciesSource[1, c(1,3)]))[5]
    colnames(infoSpeciesSource) <- c("nameFile", "gnId", "genCoord")
    dta_file <- dta_file[-1, ]
    
    # Deletion of non matching species
    dta_file <- dta_file[dta_file$genCoord != 'undefined', ]
    
    # Transform the long scientific name into a short common name
    alignementResults <- merge(dta_file, equivalenceName, by.x = "speciesLongName", by.y = "longName", all.x = T)
    alignementResults <- alignementResults[, c("shortName", "genCoord")]
    
    # Focus on the results of the "alignementOn" species
    alignementResults <- alignementResults[alignementResults$shortName %in% alignmentOn, ]
    
    ## If there is an alignement
    if (nrow(alignementResults) > 0) {
        isMatching <- 1
        # Extraction of the genomic position
        alignementResults <- data.frame(t(apply(alignementResults, 1 , extractionGenomic)), stringsAsFactors = F)
        colnames(alignementResults) <- c("shortName", "chr", "start", "end", "strand")
        alignementResults$start <- as.numeric(alignementResults$start)
        alignementResults$end <- as.numeric(alignementResults$end)
        # If the strand_SpeciesSource == "-", all the strand matching must be inverted
        # (see perl API ensembl program, and probleme betwwen gene and slice)
        if (strand_SpeciesSource == "-") {
            alignementResults$strand[alignementResults$strand == "-"] <- "tmpStr"
            alignementResults$strand[alignementResults$strand == "+"] <- "-"
            alignementResults$strand[alignementResults$strand == "tmpStr"] <- "+"
        }
        
        # If there is more than one block : 
        if (nrow(alignementResults) > 1){
            # a) If it is on very different part of the genome (different chromosome or size between blocks
            # too important (> sizeBlockAuthorized)
            if ((length(unique(alignementResults$chr)) > 1)) {
                toBind <- c(infoSpeciesSource$gnId, -1, NA, NA)
            } else {
                # b) If they are on the same part of the genome (< sizeBlockAuthorized)
                # we create a "false" unique bloc with the extreme start/end from all the blocks
                start_shifted <- alignementResults$start[2:nrow(alignementResults)] 
                end_shifted <- alignementResults$end[1:(nrow(alignementResults)-1)]
                
                toCompareToSizeBlock <- start_shifted - end_shifted
                toCompareToSizeBlock <- toCompareToSizeBlock > sizeBlockAuthorized
                
                if (any(toCompareToSizeBlock)) {
                    toBind <- c(infoSpeciesSource$gnId, -1, NA, NA)
                } else {
                    alignementResults[1, "start"] <- min(alignementResults$start)
                    alignementResults[1, "end"] <- max(alignementResults$end)
                    alignementResults <- alignementResults[1, ]
                }
            }
        }
    } 
    else {
        # If there is no alignement, only kept the lncRNA id from the "alignmentSource" species
        # and isMatching = 0.  
        toBind <- c(infoSpeciesSource$gnId, 0, NA, NA)
    }
    if (nrow(alignementResults) == 1){
        # Extraction of the gnIds and gnBiotypes from the "alignmentOn" species
        chr_toSearch <- alignementResults[1, "chr"]
        start_toSearch <- as.numeric(alignementResults[1, "start"])
        end_toSearch <- as.numeric(alignementResults[1, "end"])
        strand_toSearch <- alignementResults[1, "strand"]
        
        
        cond1 <- gtf$strand == strand_toSearch
        cond2 <- gtf$seqname %in% chr_toSearch
        cond3 <- gtf$start <= start_toSearch &  gtf$end >= start_toSearch
        cond4 <- gtf$end <= start_toSearch &  gtf$end >= end_toSearch
        
        matching <- gtf[cond1 & cond2 & (cond3 | cond4), ]
        
        # If there is no features in the "alignmentOn" species
        if (nrow(matching) == 0) {
            gnId_target <- NA
            gnBiotype_target <- NA
        } 
        # If there is onr or more features in the "alignmentOn" species (collapsing)
        else {
            gnId_target <- paste0(matching$gene_id, collapse = ";")
            gnBiotype_target <- paste0(matching$gene_biotype, collapse = ";")
        }
        
        toBind <- c(infoSpeciesSource$gnId, isMatching, gnId_target, gnBiotype_target)
    }
    # Binding the sub_results to the global table
    resultsMatching <- rbind(resultsMatching, toBind)
    colnames(resultsMatching) <- c("gnId_source", "isMatching", "gnId_target", "gnBiotype_target")
}
resultsMatching <- resultsMatching[-1, ]

#Addition of a column to find lncRNA easily
resultsMatching$isLncRNA[grepl(regEx_LNC, resultsMatching$gnBiotype_target, perl = T)] <- 1

cat("DONE\n")
# Writing to the file
cat("Writing of the file ... ")
write.table(resultsMatching, 
            paste0(alignmentSource,"_alignedTo_", alignmentOn,".tsv"),
            quote = F, sep = "\t", row.names = F)
cat("DONE\n")
cat("------------------------------------ \n\n")
### Beginning of the loop (end)
