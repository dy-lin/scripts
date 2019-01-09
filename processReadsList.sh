#!/bin/bash
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
		\?) echo "$PROGRAM: invalid option: $OPTARG" 1>&2; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ "$#" -ne 1 || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM <file list of reads>" 1>&2
	echo "DESCRIPTION: Takes a file containing a list of reads (1 per line) and outputs them 1 read pair (2 reads) per line." 1>&2
	echo -e "OPTIONS:\n\t\t-f\t\tWrite output to default file named 'reads.in'\n\t\t-o <FILENAME>\tWrite output to <FILENAME>" 1>&2
	exit 1
fi

file=$1
total=$(sort -u $file | wc -l)

if [[ $(( $total % 2 )) -ne 0 ]]
then
	echo "Some reads are missing their pair." 1>&2
	exit 1
else
	count=$((total/2))
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

if [ "$fileout" = true ]
then
	echo "Output written to $output." 1>&2
fi
