#!/bin/bash
PROGRAM=$(basename $0)
if [[ "$#" -eq 0 ]]
then
	echo "$PROGRAM <FASTQ files>" 1>&2
	exit 1
fi
sum=0
# If the first argument is a 'key', i.e. NOT a file, then use grep. If it is a file, use wc 
if [[ -e "$1" ]]
then
	for file in $@
	do
		if [[ -e "$file" ]]
		then
			if [[ "$file" == *.fastq || "$file" == *.fq  ]]
			then
				count=$(wc -l $file | awk '{print $1}')
				temp=$((count/4))
				echo "$(basename $file): $temp reads"
				sum=$((sum+temp))
			elif [[ "$file" == *.gz  ]]
			then
				count=$(pigz -d -c $file | wc -l)
				temp=$((count/4))
				sum=$((sum+temp))
				echo "$(basename $file): $temp reads"
			else
				echo "$(basename $file): Wrong file type." 1>&2
			fi
		else
			echo "$(basename $file) does not exist." 1>&2
		fi
	done
else
	key=$1
	shift
	key=^$key

	for file in $@
	do
		if [[  -e "$file" ]]
		then
			if [[ "$file" == *.gz ]]
			then
				temp=$(zgrep -c $key $file)
				sum=$((sum + temp))
			elif [[ "$file" == *.fastq || "$file" == *.fq ]]
			then
				temp=$(grep -c $key $file)
				sum=$((sum+temp))
			else
				echo "$(basename $file): Wrong file type." 1>&2
			fi
		else
			echo "$(basename $file) does not exist." 1>&2
		fi
	done
fi
echo "Total: $sum reads"
		
