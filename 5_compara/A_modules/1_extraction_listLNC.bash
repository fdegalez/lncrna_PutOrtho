#!/bin/bash
name=$1
GTF=$2

regEx_lncRNA='gene_biotype\s\"(lncRNA|lincRNA|sense_intronic|sense_exonic|antisense)\"'

grep -v "#" $GTF | awk '{if ($3 == "gene") print $0}' | grep -P $regEx_lncRNA | grep -P -o "gene_id[^;]*" | cut -d " " -f2 | sed "s/\"//g" > 1_listLNC_${name}.list