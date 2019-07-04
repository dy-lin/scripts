#!/bin/bash
PROGRAM=$(basename $0)
if [[ "$#" -eq 0 ]]
then
	echo "$PROGRAM <FASTQ files>" 1>&2
	exit 1
fi
sum=0
for file in $@
do
	if [[ "$file" == *.gz ]]
	then
	#	count=$(wc -l $file | awk '{print $1}')
	#	temp=$((count/4))
		temp=$(zgrep -c '^@HS' $file)
		sum=$((sum+temp))
		echo "$(basename $file): $temp reads"
	elif [[ "$file" == *.fastq || "$file" == *.fq ]]
	then
	#	count=$(wc -l $file | awk '{print $1}')
	#	temp=$((count/4))
		temp=$(grep -c '^@HS' $file)
		sum=$((sum+temp))
		echo "$(basename $file): $temp reads"
	else
		echo "$(basename $file): Wrong file type."
	fi
done

echo "Total: $sum reads"
		
