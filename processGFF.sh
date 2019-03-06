#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)
gethelp=false
introns=false

while getopts :hi opt
do 
	case $opt in
		h) gethelp=true;;
		i) introns=true;;
		\?) echo "$PROGRAM: Invalid option $opt" 1>&2; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ "$#" -ne 1 ||  "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM <GFF file>" 1>&2
	echo "DESCRIPTON: Takes a (maker generated) GFF file and adds introns and genes." 1>&2
	echo -e "OPTIONS:\n\t-h\tShow help menu\n\t-i\tAdd introns" 1>&2
	exit 1
fi

present=false
gff=$1
IFS=$'\n'
# in the case that genes are already in the file, and introns are selected-- nothing would happen, exit
if [[ "$(awk '/\tintron\t/' $gff | wc -l)" -gt 0 ]]
then
	introns=false
	present=true
fi

if [[ "$(awk '/\tgene\t/' $gff | wc -l)" -gt 0 ]]
then
	genes=false
else
	genes=true
fi

if [[ "$genes" = true ]] 
then
	if [[ "$introns" = true ]]
	then
		echo "Annotating genes and introns..." 1>&2
	else
		if [[ "$present" = true ]]
		then
			echo "Introns have already been annotated. Annotating genes..." 1>&2
		else
			echo "Annotating genes..." 1>&2
		fi
	fi
else
	if [[ "$introns" = true ]]
	then
		echo "Genes have already been annotated. Annotating introns..." 1>&2
	else
		if [[ "$present" = true ]]
		then
			echo "Genes and introns have already been annotated." 1>&2
			exit 1
		else
			echo "Genes have already been annotated." 1>&2
			exit 1
		fi
	fi

fi

if [ -e "temp.gff" ]
then
	rm temp.gff
fi

# if genes=true; means file has no genes-- annotate!
if [[ "$genes" = true ]]
then
	# take transcript line information and transform into gene information
	for line in $(cat $gff)
	do
		feature=$(echo $line | awk -F "\t" '{print $3}')
		if [[ "$feature" == "mRNA" ]]
		then
			# The parent attribute of mRNA is the ID of gene
			transcriptID=$(echo $line | awk -F "\t" '{print $9}' | awk -F "ID=" '{print $2}' | awk -F ";" '{print $1}')
			ID=$(echo $line | awk -F "\t" '{print $9}' | awk -F "Parent=" '{print $2}' | awk -F ";" '{print $1}')
			if [[ -z "$ID" ]]
			then
				if [[ -e "temp.gff" ]]
				then
					rm temp.gff
				fi
				echo "$transcriptID is missing the Parent attribute." 1>&2
				exit 1
			fi
			geneline=$(echo "$line" | sed 's/\tmRNA\t/\tgene\t/' | sed "s/ID=.\+$/ID=$ID/")
			if [[ "$transcriptID" == *-mRNA-1 ]]
			then
				echo $geneline >> temp.gff
			fi
			echo $line >> temp.gff
		else
			echo $line >> temp.gff
		fi	
	done
	mv temp.gff $gff
fi

if [[ "$introns" = true ]]
then
	src=$(awk '!/^#/' $gff | head -n 1 | awk '{print $2}')
	gt gff3 -addintrons -retainids -checkids -setsource $src $gff > temp.gff
	# Fix intron IDs
	mv temp.gff $gff
	for line in $(cat $gff)
	do
		feature=$(echo $line | awk -F "\t" '{print $3}')
		if [[ "$feature" == "mRNA" ]]
		then
			count=0
		fi
		if [[ "$feature" == "intron" ]]
		then
			count=$((count+1))
			parent=$(echo $line | awk -F "\t" '{print $9}' | awk -F "=" '{print $2}')
			if [[ "$count" -eq 1 ]]
			then
				newline=$(echo $line | sed "s/Parent=.\+$/ID=$parent:intron;Parent=$parent/")
			else
				newline=$(echo $line | sed "s/Parent=.\+$/ID=$parent:intron$count;Parent=$parent/")
			fi
			echo $newline >> temp.gff
		else
			echo $line >> temp.gff
		fi
	done
	mv temp.gff $gff
fi


