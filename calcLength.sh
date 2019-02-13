#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)
# Extract intron lengths given a GFF file
# Prints results to screen, use bash redirection to save file to a tsv
intron=false
gethelp=false
while getopts :hi opt
do
	case $opt in
		h) gethelp=true;;
		i) intron=true;;
		\?)	echo "$PROGRAM: Invalid option $opt" 1>&2; exit 1;;
	esac
done

shift $((OPTIND-1))
if [[ "$#" -ne 1 || "$gethelp" = true ]]
then
	echo "USAGE: $(basenmae $0) <GFF file>" 1>&2
	echo "DESCRIPTION: Takes a GFF file, extracts feature length." 1>&2
	echo -e "OPTIONS:\n\t-h\tShow help menu\n\t-i\tGet intron length" 1>&2
	exit 1
fi
gff=$1

IFS=$'\n'

for line in $(cat $gff)
do
	feature=$(echo $line | awk -F "\t" '{print $3}')
	if [ "$feature" == "mRNA" ]
	then
		name=$(echo $line | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//')
		begin=$(echo $line | awk -F "\t" '{print $4}')
		end=$(echo $line | awk -F "\t" '{print $5}')
		length=$((end-begin+1))
		echo -e "$name\t$length"
	fi
	if [ "$feature" == "intron" ]
	then
		begin=$(echo $line | awk -F "\t" '{print $4}')
		end=$(echo $line | awk -F "\t" '{print $5}')
		length=$((end-begin+1))
		echo -e "$name:intron\t$length"
	fi
done
	
