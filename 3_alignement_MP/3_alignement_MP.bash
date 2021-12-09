#!/bin/bash

### Exemple of command
### 3_alignement_MP.bash "human" "mouse" "/home/fabien/References/Human_GRCh38.p13/Ensembl/Homo_sapiens.GRCh38.101.gtf" "/home/fabien/References/Mouse_GRCm38.p6/Ensembl/Mus_musculus.GRCm38.102.gtf"

#shortNameSource="human"
#shortNameTarget="mouse"
#GTFsource="/home/fabien/References/Human_GRCh38.p13/Ensembl/Homo_sapiens.GRCh38.101.gtf"
#GTFtarget="/home/fabien/References/Mouse_GRCm38.p6/Ensembl/Mus_musculus.GRCm38.102.gtf"
#sizeBlockAuthorized=500

shortNameSource=$1
shortNameTarget=$2
GTFsource=$3
GTFtarget=$4
sizeBlockAuthorized=$5


## The alignement of lncRNAs of the main sepcies takes a lot of time, if the directory already exit, we just copy it. 
# For species 1 
sp1_file=`find . -maxdepth 2  -name  ${shortNameSource}_alignement_MP_63amniotes | head -1`
rp1=`realpath $sp1_file`

# Directory Creation
mkdir "${shortNameSource}_comparedTo_${shortNameTarget}"
cd "${shortNameSource}_comparedTo_${shortNameTarget}"


## Extraction of lncRNA of the species to map them 
bash ../../A_modules/1_extraction_listLNC.bash ${shortNameSource} ${GTFsource}

## Alignement of the lncRNAs of the species
if [ -z "$sp1_file" ]
then
      echo "The ${shortNameSource}_alignement_MP_63amniotes directory doesn't exist and will be created."
      mkdir ${shortNameSource}_alignement_MP_63amniotes
      perl ../../A_modules/2_alignementCompara_MP.pl ${shortNameSource} 1_listLNC_${shortNameSource}.list
else
      echo "The ${shortNameSource}_alignement_MP_63amniotes directory already exist and has been copied."
	  cp -r $rp1 ./
fi

## Gene extraction for the homology alignement
bash ../../A_modules/3_extractGeneGTF.bash $GTFtarget > ${shortNameTarget}_geneExtracted.gtf
Rscript ../../A_modules/3_AlignementResults.R ${shortNameSource} ${shortNameTarget} ${sizeBlockAuthorized}

cd ..