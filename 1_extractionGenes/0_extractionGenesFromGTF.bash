#!/bin/bash

# Can be used both by changing the path or in commmand line with options
configFile_path="/home/fabien/2_homology_2023/data/config.txt"
configFile_path="$1"

# Clear
rm -r results
mkdir results

for GTF in $(tail -n+2 ${configFile_path} | cut -f4)
do
echo $GTF
ouptut_name=`echo $GTF | rev | cut -d'/' -f-1 | rev | sed 's/.gtf/_genesOnly.gtf/g'`
grep -v "#" $GTF | awk -F "\t" '{
if ($3 == "gene") print $0
}' > results/${ouptut_name}
done


