# Header ---------------------------------------------------------------------------------------------
## Aim : Create FASTA sequence based on coordinates
## Date : 05/04/2024 (last update)
## Author(s): -Fabien DEGALEZ


# Working Directory ---------------------------------------------------------------------------
setwd("~/2_homology_2023/getFasta")

# Libraries -----------------------------------------------------------------------------------


# Variables -----------------------------------------------------------------------------------
bed_file <- "test.bed"
fasta_file <- "/home/fabien/Rocco/0_referenceFiles/goat_CHI/ARS1/01_genome/Capra_hircus.ARS1.dna.toplevel.fa"
output_name <- "test.fa"
strand_option <- F

# Functions -----------------------------------------------------------------------------------

## DECOMPRESS FASTA
#commandToSend <- paste("gunzip", fasta_file)
#system(command = commandToSend)


# Input --------------------------------------------------------------------------------------


# Script --------------------------------------------------------------------------------------
if (strand_option == T){
    strand <- "-s"
} else {
    strand <- ""
}


commandToSend <- paste("bedtools getfasta -fi", fasta_file, "-bed", bed_file, "-name -fo", output_name, strand)
system(command = commandToSend)


