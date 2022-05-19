#!/bin/bash

### Command usage :
#bash 1_generationMulti_synteny.bash [CONFIG FILE ABSOLUTE PATH]

######################################################################################
## Common usage :
# 1_synteny.bash \
# "sp1_Name" \
# "sp2_Name" \
# "sp1_ensName" \
# "sp2_ensName" \
# "sp1_GTFpath" \
# "sp2_GTFpath"
 
## Example : 
# 1_synteny.bash \
# "human" \
# "mouse" \
# "hsapiens" \
# "mmusculus" \
# "Ensembl/Homo_sapiens.GRCh38.101.gtf" \
# "Ensembl/Mus_musculus.GRCm38.102.gtf"
######################################################################################

## Aim of this program : From the config.txt file containing all the information by species we want to use in the analysis
# generate a file that will be launch to create all the possible intersection. 

######################################################################################

## Global Variables : 
RED=$(tput setaf 1)
NORMAL=$(tput sgr0)
start=`date`
clear

######################################################################################

# The config file is the input of the program (first argument)
config=$1 ## IN ABSOLUTE PATH !

# Number of lines in the file (= number of species + 1 )
nb=`wc -l $config | cut -d" " -f1`
((nb_sp=nb-1))
nb_done=0 

# Displayer method
printf "#####################################\n"
printf "#####################################\n"
printf "#####################################\n\n"
printf "###       SYNTENY - Method 1      ###\n\n"
printf "#####################################\n"
printf "#####################################\n"
printf "#####################################\n\n\n\n\n"
sleep 5

# Displayer nb of species
printf "##################################### \n\n"
printf "Launching of the crossed analysis for ${RED} $nb_sp species ${NORMAL} \n\n"
printf "##################################### \n\n\n"

# Move to the associated directory
cd 1_synteny

# Establishment of all crossing
# For two species, the program is launched in both configuration (sp1→sp2 / sp2→sp1)
# because depending of the target specie, annotation can change (e.g one_to_many → many_to_one, or one_to_zero cases)

for count1 in  `seq 2 $nb`
do
sp1_Name=`sed -n "${count1}p" $config | cut -f1 -d$'\t'`
sp1_ensName=`sed -n "${count1}p" $config | cut -f2 -d$'\t'`
sp1_GTFpath=`sed -n "${count1}p" $config | cut -f3 -d$'\t'`
for count2 in  `seq 2 $nb`
do

sp2_Name=`sed -n "${count2}p" $config | cut -f1 -d$'\t'`
sp2_ensName=`sed -n "${count2}p" $config | cut -f2 -d$'\t'`
sp2_GTFpath=`sed -n "${count2}p" $config | cut -f3 -d$'\t'`


if [ "$sp1_ensName" != "$sp2_ensName" ]; then
printf "##################################### \n"
((nb_done=nb_done+1))
printf "Crossing n° ${nb_done} : ${sp1_Name} & ${sp2_Name} \n"

bash ./1_synteny.bash \
${sp1_Name} \
${sp2_Name} \
${sp1_ensName} \
${sp2_ensName} \
${sp1_GTFpath} \
${sp2_GTFpath}
printf "##################################### \n\n"
fi
done
done

cd .. 
end=`date`


# Displayer end time
printf "##################################### \n\n"
printf "Total of crossing :  ${nb_done} ${NORMAL}\n"
printf "Starting time : ${RED} ${start} ${NORMAL}\n"
printf "Ending time :${RED}  ${end} ${NORMAL}\n"
printf "##################################### \n\n"