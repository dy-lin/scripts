#!/bin/bash

# /projects/btl/dlin/scripts/Rscripts/summarize.R

# should produce the summarized.csv file
# Parse this file to create file for plotting
# Gene, TPM, Treatment, Condition
condition=$1
#condition="wpw"
summary=summarized.csv
# summary=/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/summarized.csv
count=2

# write the file first
echo "Gene,TPM,Treatment,Condition" > final_summary.csv
for column in $(head -n1 $summary | tr "," " " | cut -d " " -f 2-)
do
	awk -v par="$condition" -v var="$column" -v bar="$count" -F "," 'BEGIN{OFS=","}{print $1, $bar, var, par }' <(tail -n +2 $summary) >> final_summary.csv
	count=$((count+1))
done


