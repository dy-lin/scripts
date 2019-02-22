#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)
# Extract intron lengths given a GFF file
# Prints results to screen, use bash redirection to save file to a tsv
intron=false
gethelp=false
table=false
suppress=false
while getopts :hist opt
do
	case $opt in
		h) gethelp=true;;
		i) intron=true;;
		s) suppress=true;;
		t) table=true;;
		\?)	echo "$PROGRAM: Invalid option $opt" 1>&2; exit 1;;
	esac
done

shift $((OPTIND-1))
if [[ "$#" -lt 1 || "$gethelp" = true ]]
then
	echo "USAGE: $(basenmae $0) <GFF file>" 1>&2
	echo "DESCRIPTION: Takes a GFF file, extracts feature length." 1>&2
	echo -e "OPTIONS:\n\t-h\tShow help menu\n\t-i\tGet intron lengths\n\t-s\tSuppress default output of transcript lengths\n\t-t\tPrint as table" 1>&2
	exit 1
fi

IFS=$'\n'
for gff in "$@"
do

	if [ -e temp.tsv ]
	then
		rm temp.tsv
	fi
	if [[ "$intron" = true && "$(grep -c '	intron	' $gff)" -eq 0 ]]
	then
		echo "$(basename $gff):"
		if [[ "$suppress" = true ]]
		then
			echo "Default output of transcript lengths has been suppressed."
		fi
		echo "There are no intron annotations in this GFF."
		echo
		continue
#		if [[ "$suppress" = true ]]
#		then
#			if [[ "$#" -ne 1 ]]
#			then
#				cat $gff
#				continue
#			fi
#		fi
	fi
	echo "$(basename $gff):" >> temp.tsv
	for line in $(awk '!/^#/' $gff)
	do
		feature=$(echo $line | awk -F "\t" '{print $3}')
		if [[ "$feature" == "mRNA" ]]
		then
			name=$(echo $line | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//')
			if [[ "$suppress" = false ]]
			then
				begin=$(echo $line | awk -F "\t" '{print $4}')
				end=$(echo $line | awk -F "\t" '{print $5}')
				length=$((end-begin+1))
				echo -e "$name\t$length" >> temp.tsv
			fi
		fi
		if [[ "$intron" = true && "$feature" == "intron" ]]
		then
			name=$(echo $line | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | sed 's/^ID=//')
			begin=$(echo $line | awk -F "\t" '{print $4}')
			end=$(echo $line | awk -F "\t" '{print $5}')
			length=$((end-begin+1))
			echo -e "$name\t$length" >> temp.tsv
		fi
	done
	if [[ "$#" -ne 1 ]]
	then
		head -n1 temp.tsv
	fi
	if [[ "$table" = true ]]
	then
		sort -t$'\t' -k1 <(tail -n +2 temp.tsv) | column -t -s$'\t'
	else
		sort -t$'\t' -k1 <(tail -n +2 temp.tsv)
	fi
	if [[ "$#" -ne 1 ]]
	then
		echo
	fi
done
if [ -e temp.tsv ]
then
	rm temp.tsv
fi
