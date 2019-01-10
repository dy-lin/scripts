#!/bin/bash
#set -eu -o pipefail
PROGRAM=$(basename $0)
fileout=false
output=reads.in
gethelp=false
while getopts :hfo: opt
do
	case $opt in
		h) gethelp=true;;
		f) fileout=true;;
		o) output=$OPTARG; fileout=true;;
		\?) echo "$PROGRAM: invalid option: $OPTARG" >&2; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ ( "$#" -ne 1 && "$#" -ne 2 ) || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM <file list of reads>" >&2
	echo "DESCRIPTION: Takes a file containing a list of reads (or two separate read pair files) and outputs them 1 read pair (2 reads) per line." >&2
	echo -e "OPTIONS:\n\t-f\t\tWrite output to default file named 'reads.in'\n\t-o <FILENAME>\tWrite output to <FILENAME>" >&2
	exit 1
fi
if [[ "$1" == "/dev/fd/63" ]]
then
	cp $1 copy.in
	file=copy.in
else
	file=$1
fi
if [[ "$#" -eq 1 ]]
then
	total=$(sort -u $file | wc -l)
	if [[ "$total" -eq 1 && "$(head -n 1 $file | awk '{print NF}')" -gt 2 ]]
	then
		cat $file | tr ' ' '\n' | sed '/^$/d' > tempfile
		file=tempfile
		total=$(sort -u $file | wc -l)
	fi
	if [[ $(( $total % 2 )) -ne 0 ]]
	then
		echo "Some reads are missing their pair." >&2
		exit 1
	fi
	for line in $(sort -u $file)
	do
		if [[ $(echo $line | grep -c '_[0-9]_1_001[_.]') -gt 0 ]]
		then
			pair=$(echo $line | sed 's/_1_001/_2_001/')
			if [ "$fileout" = true ]
			then
				echo "$line $pair" > $output
			else
				echo "$line $pair" 
			fi
		elif [[ $(echo $line | grep -c '_R1[_.]') -gt 0 ]]
		then
			pair=$(echo $line | sed 's/_R1/_R2/')
			if [ "$fileout" = true ]
			then
				echo "$line $pair" > $output
			else
				echo "$line $pair" 
			fi
		elif [[ $(echo $line | grep -c '[_-]first') -gt 0 ]]
		then
			pair=$(echo $line | sed 's/first/second/')
			if [ "$fileout" = true ]
			then
				echo "$line $pair" > $output
			else
				echo "$line $pair" 
			fi
		else
			continue
		fi
	done
fi
	if [[ -e tempfile ]]
	then
		rm tempfile
	fi
	if [[ -e copy.in ]]
	then
		rm copy.in
	fi
if [[ "$#" -eq 2 ]]
then
	if [ "$fileout" = true ]
	then
		mergeReads.py $1 $2 > $output
	else
		mergeReads.py $1 $2 
	fi
fi

if [ "$fileout" = true ]
then
	echo "Output written to $output." >&2
fi
