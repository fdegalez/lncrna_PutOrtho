#!/bin/bash

# Mandatory input: GTF and fasta files
bed_file=$1
fasta_file=$2
output_name=$3
strand_option=$4
# Optional parameters


if [[ $fasta_file == *.gz ]]; then
  echo "Le fichier est compressé. Il va être décompressé."
  gunzip -c $fasta_file > temp.fasta
  mv temp.fasta ${fasta_file%%.*}
  fasta_file=${fasta_file%%.*}
fi


# Output the fasta sequence using bedtools getFasta function
bedtools getfasta $strand_option $output_name -fi $fasta_file -bed $bed_file $output_name

bedtools getfasta -fi $fasta_file -bed $bed_file -name -fo $output_name