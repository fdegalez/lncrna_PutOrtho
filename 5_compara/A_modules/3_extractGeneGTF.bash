#!/bin/bash

GTF="$1"
grep -v "#" $GTF | awk -F "\t" '{
if ($3 == "gene") print $0
}'