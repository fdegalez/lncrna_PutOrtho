#!/bin/bash

#####
##### Extraction lncRNA & mRNA + FEELncClassifier + geneLevel annotation
#####

## Aim of the proram :
# Extraction of the correct lncRNA features from the GTF and application
# of the FEELnc_classifier program

## Inputs :
GTF=$1 #Necessary to work
name=$2 #Suffix which will be use to name the output files

## Outputs : 
# {name}_feelncclassifier.log - general statistics on the number of interactions
# {name}_classes_feelncclassifier.txt - tabulated-format file with all the interactions
# {name}_lncConfiguration_feelncclassifier.tsv - tabulated file with all the configuration at the gene level

#####
#####
#####

## Configuration of the FEELnc environment

## To change if necessary : 
cd /home/fabien/00_SOFTWARE/FEELnc #FEELNC directory
FEELnc_classifier="/home/fabien/00_SOFTWARE/FEELnc/scripts/FEELnc_classifier.pl" #FEELNC_classifier path

export FEELNCPATH=${PWD}
export PERL5LIB=$PERL5LIB:${FEELNCPATH}/lib/
export PATH=$PATH:${FEELNCPATH}/scripts/
export PATH=$PATH:${FEELNCPATH}/utils/
export PATH=$PATH:${FEELNCPATH}/bin/LINUX/

cd -

## Extraction the lncRNA features from the GTF
##To change if necessary

regEx_lncRNA='gene_biotype\s\"(lncRNA|lincRNA|sense_intronic|sense_exonic|antisense)\"'
regEx_mRNA='gene_biotype\s\"protein_coding\"'

#regEx_lncRNA='gene_biotype\s(lncRNA|lincRNA|sense_intronic|sense_exonic|antisense)'
#regEx_mRNA='gene_biotype\sprotein_coding'

#lncRNA
echo "------------------------------------"
printf "Extraction of lncRNA ... "
zgrep -v "#" $GTF | grep -P $regEx_lncRNA > ${name}_LNCextracted.tmp.gtf
printf "DONE\n"
#mRNA 
printf "Extraction of mRNA ... "
zgrep -v "#" $GTF | grep -P $regEx_mRNA > ${name}_mRNAextracted.tmp.gtf
printf "DONE\n"

#lncRNA relative to mRNA
printf "FEELnc_classifier ... \n "
$FEELnc_classifier \
-i ${name}_LNCextracted.tmp.gtf \
-a ${name}_mRNAextracted.tmp.gtf \
-l ${name}_feelncclassifier.log > ${name}_classes_feelncclassifier.txt
#rm *.tmp.gtf
printf "FEELnc_classifier ... DONE"

#geneLevelAnnotation
printf "Creation of the gene level annotation of lncRNA ... "
Rscript ../../A_modules/0_FEELnc_tpLevel2gnLevelClassification.R ${name} ${GTF}
printf "DONE\n"
echo "------------------------------------"

rm *tmp* 