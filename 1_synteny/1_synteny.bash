#!/bin/bash

### Exemple of command
### 1_synteny.bash "human" "mouse" "hsapiens" "mmusculus" "/home/fabien/References/Human_GRCh38.p13/Ensembl/Homo_sapiens.GRCh38.101.gtf" "/home/fabien/References/Mouse_GRCm38.p6/Ensembl/Mus_musculus.GRCm38.102.gtf"
 

## To work 

#shortNameSource="human"
#shortNameTarget="mouse"
#ensemblNameSource="hsapiens"
#ensemblNameTarget="mmusculus"
#GTFsource="/home/fabien/References/Human_GRCh38.p13/Ensembl/Homo_sapiens.GRCh38.101.gtf"
#GTFtarget="/home/fabien/References/Mouse_GRCm38.p6/Ensembl/Mus_musculus.GRCm38.102.gtf"

##
shortNameSource=$1
shortNameTarget=$2
ensemblNameSource=$3
ensemblNameTarget=$4
GTFsource=$5
GTFtarget=$6

## Directory Creation

mkdir "${shortNameSource}_comparedTo_${shortNameTarget}"
cd "${shortNameSource}_comparedTo_${shortNameTarget}"

## Creation of the homology file between the two species 

echo $ensemblNameSource $ensemblNameTarget $shortNameSource $shortNameTarget
Rscript ../A_modules/0_creationBiomartFile.R $ensemblNameSource $ensemblNameTarget $shortNameSource $shortNameTarget

## Extraction of the gene from the GTF of each species
# • Species Source
../A_modules/1_extractGeneGTF.bash $GTFsource > ${shortNameSource}_geneExtracted.gtf
# • Species Target
../A_modules/1_extractGeneGTF.bash $GTFtarget > ${shortNameTarget}_geneExtracted.gtf

## Creation of the table (lncRNA position compared to the two nearest PCGs) for each species
# • Species Source
Rscript ../A_modules/2_creationTable.R ${shortNameSource}_geneExtracted.gtf ${shortNameSource}
# • Species Target
Rscript ../A_modules/2_creationTable.R ${shortNameTarget}_geneExtracted.gtf ${shortNameTarget}

## Creation of the homology file
Rscript ../A_modules/3_LNChomologyDetection.R \
2_tableGTF_${shortNameSource}.tsv \
2_tableGTF_${shortNameTarget}.tsv \
homologyFile.tsv \
${shortNameSource} \
${shortNameTarget} 

## Concatenation of the homology file and adding orthology type and configuration of strand
Rscript ../A_modules/4_orthologyTypeAndConfiguration.R \
3_${shortNameSource}_${shortNameTarget}_orthology.tsv \
${shortNameSource} \
${shortNameTarget} 

cd ..
mv "${shortNameSource}_comparedTo_${shortNameTarget}" B_results