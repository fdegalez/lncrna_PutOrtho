#!/usr/bin/env Rscript

# ==============================================================================
# Script: syntenyBySpecies.R
# Aim:
#   Identify putative lncRNA orthologs between species based on flanking PCGs
# ==============================================================================


# Libraries --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(stringr)
  library(pbapply)
})


# Parameters -------------------------------------------------------------------

config_path <- "../data/config.txt"

results_table_dir <- "results_table"
results_synteny_dir <- "results_synteny"

dir.create(results_synteny_dir, showWarnings = FALSE)


# Orientation table ------------------------------------------------------------

strand_case_possibility <- data.frame(
  combinaison = c("+++","++-","+-+","-++","+--","-+-","--+","---"),
  group = c(1,2,3,4,4,3,2,1),
  stringsAsFactors = FALSE
)


# Functions --------------------------------------------------------------------

test_PCG <- function(source,target){

  if(source==target) return("same")

  tmp <- str_replace_all(source,"\\+","A")
  tmp <- str_replace_all(tmp,"\\-","+")
  tmp <- str_replace_all(tmp,"A","-")

  if(tmp==target) return("reverse")

  return("discordant")

}

find_lncRNA_orth_byPCGcouple <- function(x){
  
  # Source PCGs
  PCG_left  <- as.character(x[1])
  PCG_right <-  as.character(x[2])
  
  # Source table: lncRNAs bounded by this exact PCG pair
  table_source_tmp <- source_table[
    source_table[,2] %in% PCG_left &
      source_table[,3] %in% PCG_right,
  ]
  
  colnames_table_source <- colnames(table_source_tmp)
  
  # Orthologous PCGs in target species
  PCG_orthologous_left  <- homology[match(PCG_left,  homology[,1]), 2]
  PCG_orthologous_right <- homology[match(PCG_right, homology[,1]), 2]
  
  # Target table: lncRNAs associated with the orthologous PCG pair
  table_target_tmp <- target_table[
    ((target_table[,2] %in% PCG_orthologous_left)  | (target_table[,3] %in% PCG_orthologous_left)) &
      ((target_table[,2] %in% PCG_orthologous_right) | (target_table[,3] %in% PCG_orthologous_right)),
  ]
  
  colnames_table_target <- colnames(table_target_tmp)
  
  # PCG strand information
  PCG_source <- paste0(c(unique(table_source_tmp[,5]), unique(table_source_tmp[,6])), collapse = "")
  PCG_target <- paste0(c(unique(table_target_tmp[,5]), unique(table_target_tmp[,6])), collapse = "")
  
  # Source lncRNAs split by strand
  table_source_toPaste <- table_source_tmp[, c(1,4,7,8)]
  
  table_source_toPaste_forward <- table_source_toPaste[table_source_toPaste[,2] == "+", , drop = FALSE]
  table_source_toPaste_reverse <- table_source_toPaste[table_source_toPaste[,2] == "-", , drop = FALSE]
  
  if (nrow(table_source_toPaste_forward) > 0) {
    table_source_toPaste_forward <- t(apply(table_source_toPaste_forward, 2, function(z) paste0(z, collapse = ";")))
    table_source_toPaste_forward <- cbind(table_source_tmp[1, c(2,3,5,6)], table_source_toPaste_forward)
    table_source_toPaste_forward <- table_source_toPaste_forward[, c(5,1,2,6,3,4,7,8)]
    colnames(table_source_toPaste_forward) <- colnames_table_source
  } else {
    table_source_toPaste_forward <- data.frame(t(rep(NA, ncol(table_source_tmp))), stringsAsFactors = FALSE)
    colnames(table_source_toPaste_forward) <- colnames_table_source
  }
  
  if (nrow(table_source_toPaste_reverse) > 0) {
    table_source_toPaste_reverse <- t(apply(table_source_toPaste_reverse, 2, function(z) paste0(z, collapse = ";")))
    table_source_toPaste_reverse <- cbind(table_source_tmp[1, c(2,3,5,6)], table_source_toPaste_reverse)
    table_source_toPaste_reverse <- table_source_toPaste_reverse[, c(5,1,2,6,3,4,7,8)]
    colnames(table_source_toPaste_reverse) <- colnames_table_source
  } else {
    table_source_toPaste_reverse <- data.frame(t(rep(NA, ncol(table_source_tmp))), stringsAsFactors = FALSE)
    colnames(table_source_toPaste_reverse) <- colnames_table_source
  }
  
  # Target lncRNAs split by strand
  table_target_toPaste <- table_target_tmp[, c(1,4,7,8), drop = FALSE]
  
  table_target_toPaste_forward <- table_target_toPaste[table_target_tmp[,4] == "+", , drop = FALSE]
  table_target_toPaste_reverse <- table_target_toPaste[table_target_tmp[,4] == "-", , drop = FALSE]
  
  if (nrow(table_target_toPaste_forward) > 0) {
    table_target_toPaste_forward <- t(apply(table_target_toPaste_forward, 2, function(z) paste0(z, collapse = ";")))
    table_target_toPaste_forward <- cbind(table_target_tmp[1, c(2,3,5,6)], table_target_toPaste_forward)
    table_target_toPaste_forward <- table_target_toPaste_forward[, c(5,1,2,6,3,4,7,8)]
    colnames(table_target_toPaste_forward) <- colnames_table_target
  } else {
    table_target_toPaste_forward <- data.frame(t(rep(NA, ncol(table_target_tmp))), stringsAsFactors = FALSE)
    colnames(table_target_toPaste_forward) <- colnames_table_target
  }
  
  if (nrow(table_target_toPaste_reverse) > 0) {
    table_target_toPaste_reverse <- t(apply(table_target_toPaste_reverse, 2, function(z) paste0(z, collapse = ";")))
    table_target_toPaste_reverse <- cbind(table_target_tmp[1, c(2,3,5,6)], table_target_toPaste_reverse)
    table_target_toPaste_reverse <- table_target_toPaste_reverse[, c(5,1,2,6,3,4,7,8)]
    colnames(table_target_toPaste_reverse) <- colnames_table_target
  } else {
    table_target_toPaste_reverse <- data.frame(t(rep(NA, ncol(table_target_tmp))), stringsAsFactors = FALSE)
    colnames(table_target_toPaste_reverse) <- colnames_table_target
  }
  
  # Symmetry case
  if (nrow(table_target_tmp) > 0) {
    if (length(unique(table_target_tmp[,3])) == 1 &&
        !is.na(PCG_orthologous_left) &&
        PCG_orthologous_left == unique(table_target_tmp[,3])) {
      PCG_target <- paste(rev(strsplit(PCG_target, "")[[1]]), collapse = "")
    }
  }
  
  strand_couple <- test_PCG(PCG_source, PCG_target)
  
  # Merge source / target according to PCG orientation
  if (strand_couple == "reverse") {
    res_forward <- cbind(table_source_toPaste_forward, table_target_toPaste_reverse)
    res_reverse <- cbind(table_source_toPaste_reverse, table_target_toPaste_forward)
  } else {
    res_forward <- cbind(table_source_toPaste_forward, table_target_toPaste_forward)
    res_reverse <- cbind(table_source_toPaste_reverse, table_target_toPaste_reverse)
  }
  
  res <- rbind(res_forward, res_reverse)
  res$PCG_strandCouple <- strand_couple
  
  return(res)
}

orthology_categories <- function(x){

  nb_source <- str_count(x[1],";")+1
  nb_target <- str_count(x[9],";")+1

  if(nb_source==1 & is.na(nb_target)) return("one_to_zero")
  if(nb_source>1 & is.na(nb_target)) return("many_to_zero")
  if(nb_source==1 & nb_target==1) return("one_to_one")
  if(nb_source==1 & nb_target>1) return("one_to_many")
  if(nb_source>1 & nb_target==1) return("many_to_one")

  if(nb_source>1 & nb_target>1){

    if(nb_source==nb_target){
      return("many_to_many_sameNumber")
    }else{
      return("many_to_many_DifferentNumber")
    }

  }

  return("Error")

}




# Load config ------------------------------------------------------------------

config <- read.delim(config_path, stringsAsFactors = FALSE)

pairs <- expand.grid(
  config$ensemblName,
  config$ensemblName,
  stringsAsFactors = FALSE
)

colnames(pairs) <- c("source","target")

pairs$source_full <- config$completeName[match(pairs$source,config$ensemblName)]
pairs$target_full <- config$completeName[match(pairs$target,config$ensemblName)]

pairs <- pairs[pairs$source!=pairs$target,]


# Main loop --------------------------------------------------------------------

i<-1
for(i in seq_len(nrow(pairs))){

  cat(i,"/",nrow(pairs),":",pairs$source[i],"-",pairs$target[i],"\n")

  source <- pairs$source[i]
  target <- pairs$target[i]

  source_full <- pairs$source_full[i]
  target_full <- pairs$target_full[i]


  # Load data -------------------------------------------------------------

  source_table <- read.delim(
    paste0(results_table_dir,"/",source_full,"_lncRNAbetweenPcg.tsv"),
    stringsAsFactors = FALSE
  )

  target_table <- read.delim(
    paste0(results_table_dir,"/",target_full,"_lncRNAbetweenPcg.tsv"),
    stringsAsFactors = FALSE
  )

  homology <- read.delim(
    paste0("../2_extractionOrthologyPCG/results/",source,"-",target,"_homology.tsv"),
    stringsAsFactors = FALSE
  )


  homology <- homology[homology[,3]=="ortholog_one2one",]


  # Identify lncRNAs with orthologous PCGs ------------------------------------

  correspondance <- data.frame(
    source_table[,1],
    source_table[,2] %in% homology[,1],
    source_table[,3] %in% homology[,1]
  )

  colnames(correspondance) <- c(
    "lncRNA_source_id",
    "PCG_left_source_orthologous",
    "PCG_right_source_orthologous"
  )


  lncRNA_concerned <- correspondance[
    correspondance$PCG_left_source_orthologous &
    correspondance$PCG_right_source_orthologous,
    "lncRNA_source_id"
  ]


  source_table <- source_table[source_table[,1] %in% lncRNA_concerned,]

  source_table <- source_table[source_table[,2]!=source_table[,3],]


  PCG_couples <- unique(source_table[,c(2,3)])


  # Orthology inference --------------------------------------------------------

  results <- pbapply(
    PCG_couples,
    1,
    find_lncRNA_orth_byPCGcouple
  )

  results <- do.call(rbind,results)


  results <- results[!apply(results,1,function(x) all(is.na(x))),]

  results <- results[!is.na(results[,1]),]
  results$type <- pbapply(results,1,orthology_categories)


  write.table(
    results,
    paste0(
      results_synteny_dir,"/",
      source,"-",target,"_synteny.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

}
