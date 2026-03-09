#!/usr/bin/env Rscript

###############################################################################
# Pairwise summary of lncRNA orthology detected by the 3 methods
###############################################################################

suppressPackageStartupMessages({
    library(stringr)
    library(pbapply)
})

###############################################################################
# Load configuration
###############################################################################

config <- read.delim(
    "../data/config.txt",
    header = TRUE,
    stringsAsFactors = FALSE
)

sp_list <- config$completeName

dir.create("pairwise/", showWarnings = FALSE)

###############################################################################
# Function
###############################################################################

detect_methods <- function(lnc_id, use_align){
    
    ########################################
    # Method 1 : synteny
    ########################################
    
    isTable_1 <- 0
    orthology_type_synteny <- NA
    lncRNA_synteny_id <- NA
    
    
    
    tmp <- file_synteny[grepl(lnc_id, file_synteny[,1]), ]
    
    if(nrow(tmp) > 0){
        isTable_1 <- 1
        orthology_type_synteny <- tmp$type[1]
        lncRNA_synteny_id <- tmp[,9][1]
    }
    
    ########################################
    # Method 2 : configuration
    ########################################
    
    isTable_2 <- 0
    orthology_type_configuration <- NA
    configuration_type <- NA
    lncRNA_configuration_id <- NA
    
    tmp <- file_configuration[grepl(lnc_id, file_configuration[,5]), ]
    
    if(nrow(tmp) > 0){
        isTable_2 <- 1
        orthology_type_configuration <- tmp[,11][1]
        configuration_type <- tmp[,4][1]
        lncRNA_configuration_id <- tmp[,8][1]
    }
    
    ########################################
    # Method 3 : alignment
    ########################################
    
    isTable_3 <- NA
    
    if(use_align == 1){
        if(lnc_id %in% file_alignment[,1]){
            isTable_3 <- 1
        } else {
            isTable_3 <- 0
        }
    }
    
    return(c(
        lnc_id,
        isTable_1,
        orthology_type_synteny,
        lncRNA_synteny_id,
        isTable_2,
        orthology_type_configuration,
        configuration_type,
        lncRNA_configuration_id,
        isTable_3
    ))
}

###############################################################################
# Pairwise loop
###############################################################################


sp1 <- sp_list[1]

for(sp1 in sp_list){
    
    sp2 <- sp_list[sp_list != sp1][1]
    
    for(sp2 in sp_list[sp_list != sp1]){
        
        cat(sp1, "-", sp2, "\n")
        
        sp1_ens <- config$ensemblName[config$completeName == sp1]
        sp2_ens <- config$ensemblName[config$completeName == sp2]
        
        ########################################
        # Input paths
        ########################################
        
        synteny_path <- paste0(
            "../3_synteny/results_synteny/",
            sp1_ens,"-",sp2_ens,"_synteny.tsv"
        )
        
        configuration_path <- paste0(
            "../4_FEELnc/results_orthoFeelnc_agg/custom/",
            sp1_ens,"-",sp2_ens,
            "_lncConfigurationHomologyAggregated.tsv"
        )
        
        alignment_path <- paste0(
            "../5_compara/results_isMatching/",
            sp1,"_isMatching.tsv"
        )
        
        ########################################
        # Check required files
        ########################################
        
        if(!file.exists(synteny_path) || !file.exists(configuration_path)){
            warning("Missing files for ", sp1,"-",sp2)
            next
        }
        
        ########################################
        # Load files
        ########################################
        
        file_synteny <- read.delim(
            synteny_path,
            header = TRUE,
            stringsAsFactors = FALSE
        )
        
        table(file_synteny$type)
        file_synteny <- file_synteny[!grepl("zero", file_synteny$type),]
        
        file_configuration <- read.delim(
            configuration_path,
            header = TRUE,
            stringsAsFactors = FALSE
        )
        
        ########################################
        # Alignment file
        ########################################
        
        use_align <- 0
        
        if(file.exists(alignment_path)){
            
            file_alignment <- read.delim(
                alignment_path,
                header = TRUE,
                stringsAsFactors = FALSE
            )
            
            if(any(grepl(sp2, colnames(file_alignment)))){
                
                num_col <- grep(sp2, colnames(file_alignment))
                
                file_alignment <- file_alignment[
                    file_alignment[,num_col] == 1,
                    c(1,num_col)
                ]
                
                use_align <- 1
            }
        }
        
        ########################################
        # Build lncRNA list
        ########################################
        
        list_meth1 <- unique(unlist(str_split(file_synteny[,1],";")))
        list_meth1 <- list_meth1[list_meth1!=""]
        
        list_meth2 <- unique(unlist(str_split(file_configuration[,5],";")))
        list_meth2 <- list_meth2[list_meth2!=""]
        
        if(use_align==1){
            
            list_meth3 <- unique(file_alignment[,1])
            
            lncRNA_list <- unique(
                c(list_meth1,list_meth2,list_meth3)
            )
            
        } else {
            
            lncRNA_list <- unique(
                c(list_meth1,list_meth2)
            )
            
        }
        
        ########################################
        # Detection
        ########################################
    
        res <- pbsapply(
            lncRNA_list,
            detect_methods,
            use_align = use_align
        )
        
        
        res <- data.frame(t(pbsapply(
            lncRNA_list,
            detect_methods,
            use_align = use_align
        )), stringsAsFactors = F)
        
        ########################################
        # Column names
        ########################################
        
        colnames(res) <- c(
            paste0("lncRNAid_",sp1),
            paste0("isTable_1_",sp2),
            paste0("orthology_type_synteny_",sp2),
            paste0("lncRNA_synteny_id_",sp2),
            paste0("isTable_2_",sp2),
            paste0("orthology_type_configuration_",sp2),
            paste0("configuration_type_configuration_",sp2),
            paste0("lncRNA_configuration_id_",sp2),
            paste0("isTable_3_",sp2)
        )
        
        ########################################
        # Write output
        ########################################
        
        write.table(
            res,
            paste0(
                "pairwise/",
                sp1,"_",sp2,
                "_summary.tsv"
            ),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE
        )
        
    }
}

