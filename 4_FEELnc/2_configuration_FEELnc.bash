#!/bin/bash

## Example : 

### 2_configuration_FEELnc.bash "human" "mouse" "hsapiens" "mmusculus" "/home/fabien/References/Human_GRCh38.p13/Ensembl/Homo_sapiens.GRCh38.101.gtf" "/home/fabien/References/Mouse_GRCm38.p6/Ensembl/Mus_musculus.GRCm38.102.gtf" "open1"

## To work 


######################################################################################
## Example : 
# 2_configuration_FEELnc.bash \
#"human" \
#"mouse" \
#"hsapiens" \
#"mmusculus" \
#"Ensembl/Homo_sapiens.GRCh38.101.gtf" \
#"Ensembl/Mus_musculus.GRCm38.102.gtf" \
#"open2"

## To debug (if needed) :
#shortNameSource="human"
#shortNameTarget="mouse"
#ensemblNameSource="hsapiens"
#ensemblNameTarget="mmusculus"
#GTFsource="/home/fabien/References/Human_GRCh38.p13/Ensembl/Homo_sapiens.GRCh38.101.gtf"
#GTFtarget="/home/fabien/References/Mouse_GRCm38.p6/Ensembl/Mus_musculus.GRCm38.102.gtf"
#configNamming="open"
######################################################################################

## Variables
shortNameSource=$1
shortNameTarget=$2
ensemblNameSource=$3
ensemblNameTarget=$4
GTFsource=$5
GTFtarget=$6
configNamming=$7 #Degree of precision for the configuration : 
					# strict : lincSSup != lincSSdw != lncgSS
					# inter : (linSSup = lincSSdw) != lngSS
					# open : lincSSup = lincSSdw = lncgSS

######################################################################################

## FELLnc takes a lot of time, if the file already exit, we just copy it. 
# For species 1 
sp1_file=`find . -name  ${shortNameSource}_lncConfiguration_feelncclassifier.tsv | head -1`
# For species 2
sp2_file=`find . -name  ${shortNameTarget}_lncConfiguration_feelncclassifier.tsv | head -1`

## Creation of the FeelNC file for both species
if ! [ -z "$sp1_file" ]
then
      rp1=`realpath -q  $sp1_file`
fi
# For species 2
if ! [ -z "$sp2_file" ]
then
      rp2=`realpath -q  $sp2_file`
fi


# Directory Creation
mkdir "${shortNameSource}_comparedTo_${shortNameTarget}"
cd "${shortNameSource}_comparedTo_${shortNameTarget}"

## Creation of the FeelNC file for both species
if [ -z "$sp1_file" ]
then
      printf "The ${shortNameSource} file doesn't exist and will be created. \n"
	  bash ../../A_modules/0_FEELnc_classifier_adapted.bash ${GTFsource} ${shortNameSource}
else
      printf "The ${shortNameSource} file already exist and has been copied. \n"
	  cp $rp1 ./
fi
# For species 2
if [ -z "$sp2_file" ]
then
      printf "The ${shortNameTarget} file doesn't exist and will be created. \n"
	  bash ../../A_modules/0_FEELnc_classifier_adapted.bash ${GTFtarget} ${shortNameTarget}
else
      printf "The ${shortNameTarget} file already exist and has been copied. \n"
	  cp $rp2 ./
fi

## Creation of the homology file
Rscript ../../A_modules/1_homologyFiles.R  ${ensemblNameSource} ${ensemblNameTarget} ${shortNameSource} ${shortNameTarget}

## Creation of the homology configuration files
Rscript ../../A_modules/2_configurationHomology.R ${shortNameSource} ${shortNameTarget} ${ensemblNameSource} ${ensemblNameTarget} ${configNamming}

## Concatenation of the different cases
Rscript ../../A_modules/3_concatenation.R ${shortNameSource} ${shortNameTarget}

cd ..