#!/bin/bash

### The command must be used like that :
#bash 2_generationMulti_configurationFEELnc.bash [CONFIG FILE ABSOLUTE PATH] [configNamming]

## Common usage :
# 2_configuration_FEELnc.bash \
#"sp1_Name" \
#"sp2_Name" \
#"sp1_ensName" \
#"sp2_ensName" \
#"sp1_GTFpath" \
#"sp2_GTFpath" \
#"configNamming"

## Example : 
### 2_configuration_FEELnc.bash \
#"human" \
#"mouse" \
#"hsapiens" \
#"mmusculus" \
#"/home/fabien/References/Human_GRCh38.p13/Ensembl/Homo_sapiens.GRCh38.101.gtf" \
#"/home/fabien/References/Mouse_GRCm38.p6/Ensembl/Mus_musculus.GRCm38.102.gtf" \
#"open2"

## Aim of this program : From the config.txt file containing all the information by species we want to use in the analysis
# generate a file that will be launch to create all the possible intersection. 

# The config file is the input of the program (first argument)
config=$1 ## IN ABSOLUTE PATH !
configNamming=$2

# Time
start=`date`

# Number of lines in the file (= number of species + 1 )
nb=`wc -l $config | cut -d" " -f1`
nb_done=0 
echo "Launching of the crossed analysis ... "

cd 2_configuration_FEELnc/B_results

for count1 in  `seq 2 $nb`
do
echo "--------------------------"
sp1_Name=`sed -n "${count1}p" $config | cut -f1 -d$'\t'`
sp1_ensName=`sed -n "${count1}p" $config | cut -f2 -d$'\t'`
sp1_GTFpath=`sed -n "${count1}p" $config | cut -f3 -d$'\t'`
let count1p1=$count1+1 
for count2 in  `seq 2 $nb`
do
sp2_Name=`sed -n "${count2}p" $config | cut -f1 -d$'\t'`
sp2_ensName=`sed -n "${count2}p" $config | cut -f2 -d$'\t'`
sp2_GTFpath=`sed -n "${count2}p" $config | cut -f3 -d$'\t'`


if [ "$sp1_Name" != "$sp2_Name" ]; then
    ((nb_done=nb_done+1))
    echo "Crossing n° ${nb_done} : ${sp1_Name} & ${sp2_Name}"

    bash ../2_configuration_FEELnc.bash \
    ${sp1_Name} \
    ${sp2_Name} \
    ${sp1_ensName} \
    ${sp2_ensName} \
    ${sp1_GTFpath} \
    ${sp2_GTFpath} \
    ${configNamming}

fi
done
done

cd .. 
end=`date`

echo "--------------------------"
echo "Total of crossing :  ${nb_done}"
echo "Starting time : ${start}"
echo "Ending time : ${end}"
echo "--------------------------"