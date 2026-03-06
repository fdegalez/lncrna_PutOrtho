#!/usr/bin/env Rscript

# ==============================================================================
# Script: creationTableLncRNAbetweenPCG.R
#
# Aim:
#   For each species, identify the closest PCG upstream and downstream
#   of each lncRNA.
#
# Input:
#   ../1_extractionGenes/results/*genesOnly.gtf
#
# Output:
#   results_table/*_lncRNAbetweenPcg.tsv
# ==============================================================================


# Libraries --------------------------------------------------------------------

suppressPackageStartupMessages({
    library(stringr)
    library(pbapply)
})


# Biotype definitions ----------------------------------------------------------

PCG.regex <- c("protein_coding")

lncRNA.regex <- c(
    "lncRNA","lincRNA",
    "antisense",
    "sense_overlapping",
    "sense_intronic"
)


# Functions --------------------------------------------------------------------

extract_attr <- function(attr, field){
    
    pattern <- paste0("(^|;\\s*)", field, "\\s+([^;]+)")
    m <- regexec(pattern, attr)
    res <- regmatches(attr, m)
    
    sapply(res, function(x){
        if(length(x)>=3) trimws(x[3]) else NA
    })
    
}



pcgIdEachSide <- function(y, GTF){
    
    chrLNC   <- y["seqname"]
    startLNC <- as.numeric(y["start"])
    endLNC   <- as.numeric(y["end"])
    
    
    leftPCG <- GTF[
        GTF$seqname==chrLNC &
            GTF$start < startLNC &
            GTF$simpleBiotype=="pcg", ]
    
    if(nrow(leftPCG)>0){
        
        leftPCG$diff <- startLNC-leftPCG$start
        leftPCG <- leftPCG[leftPCG$diff==min(leftPCG$diff),][1,]
        
    }else{
        
        leftPCG <- data.frame(gene_id=NA,strand=NA,diff=NA)
        
    }
    
    
    rightPCG <- GTF[
        GTF$seqname==chrLNC &
            GTF$start >= startLNC &
            GTF$simpleBiotype=="pcg", ]
    
    if(nrow(rightPCG)>0){
        
        rightPCG$diff <- abs(startLNC-rightPCG$start)
        rightPCG <- rightPCG[rightPCG$diff==min(rightPCG$diff),][1,]
        
    }else{
        
        rightPCG <- data.frame(gene_id=NA,strand=NA,diff=NA)
        
    }
    
    
    return(c(
        y["gene_id"],
        leftPCG$gene_id,
        rightPCG$gene_id,
        y["strand"],
        leftPCG$strand,
        rightPCG$strand,
        leftPCG$diff,
        rightPCG$diff
    ))
    
}



# Input ------------------------------------------------------------------------

dir.create("results_table", showWarnings=FALSE)

gtf_files <- list.files(
    "../1_extractionGenes/results/",
    pattern="genesOnly.gtf$",
    full.names=TRUE
)


# Main loop --------------------------------------------------------------------
gtf_path <- gtf_files[1]
for(gtf_path in gtf_files){
    
    species <- str_split(str_remove(basename(gtf_path), "_genesOnly.gtf"), "\\.", simplify = T)[,1]
    
    cat("Processing:",species,"\n")
    
    GTF <- read.delim(gtf_path,header=FALSE,stringsAsFactors=FALSE)
    
    colnames(GTF) <- c(
        "seqname","source","feature",
        "start","end","score",
        "strand","frame","attribute"
    )
    
    
    GTF$gene_id <- extract_attr(GTF$attribute,"gene_id")
    
    GTF$gene_biotype <- extract_attr(GTF$attribute,"gene_biotype")
    
    
    simpleBiotype <- rep(NA,nrow(GTF))
    
    simpleBiotype[GTF$gene_biotype %in% lncRNA.regex] <- "lnc"
    simpleBiotype[GTF$gene_biotype %in% PCG.regex]    <- "pcg"
    
    GTF$simpleBiotype <- simpleBiotype
    
    
    GTF <- GTF[GTF$simpleBiotype %in% c("lnc","pcg"),]
    
    listLNC <- GTF[GTF$simpleBiotype=="lnc",]
    
    
    boundedLNC <- data.frame(
        t(pbapply(
            listLNC,
            1,
            function(row) pcgIdEachSide(row,GTF)
        )),
        stringsAsFactors=FALSE
    )
    
    
    colnames(boundedLNC) <- c(
        paste0("lncRNA.",species),
        paste0("PCG_left.",species),
        paste0("PCG_right.",species),
        paste0("lncRNA_strand.",species),
        paste0("PCG_left_strand.",species),
        paste0("PCG_right_strand.",species),
        paste0("lncRNA_PCG_left_distance.",species),
        paste0("lncRNA_PCG_right_distance.",species)
    )
    
    
    write.table(
        boundedLNC,
        file=paste0("results_table/",species,"_lncRNAbetweenPcg.tsv"),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )
    
}
