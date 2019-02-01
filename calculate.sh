#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)
average=false
median=false
mode=false

while getopts :amo opt
do
	case $opt in
		a) average=true;;
		m) median=true;;
		o) mode=true;;
		\?) echo "$PROGRAM: invalid option: $opt" >&2; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $PROGRAM <FILE>" 1>&2
	echo "DESCRIPTION: Takes a file with numbers, and prints the desired calculation result to STDOUT." 1>&2
	echo -e "OPTIONS: At least ONE option must be selected.\n\t-a\tmean (average) \n\t-m\tmedian\n\t-o\tmode"
	exit 1
fi

if [[ "$average" = false && "$median" = false && "$mode" = false ]]
then
	echo "USAGE: $PROGRAM <FILE>" 1>&2
	echo "DESCRIPTION: Takes a file with numbers, and prints the mean/median/mode (depending on flags) to STDOUT." 1>&2
	echo -e "OPTIONS: At least ONE option must be selected.\n\t-a\tmean (average) \n\t-m\tmedian\n\t-o\tmode"
	exit 1
fi
file=$1
cp $file temp.out
file=temp.out
linecount=$(wc -l $file | awk '{print $1}')
if [[ "$median" = true ]]
then
	# if numbers all on one line, e.g. an array
	if [[ "$linecount" -eq 1 ]]
	then
		echo "Array detected! Calculating median..." 1>&2
		cat $file | tr ' ' '\n' | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print "Median: " (a[x-1]+a[x])/2; else print "Median: " a[x-1]; }'
	else
		echo "Calculating median..." 1>&2
		sort -n $file | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print "Median: " (a[x-1]+a[x])/2; else print "Median: " a[x-1]; }'
	fi
fi
if [[ "$average" = true ]]
then
	if [[ "$linecount" -eq 1 ]]
	then
		echo "Array detected! Calculating mean..." 1>&2
		cat $file | awk '{print}' | tr ' ' '\n' | awk '{sum+=$1} END{ print "Mean: " sum/NR}'
	else
		echo "Calculating mean..." 1>&2
		awk '{sum+=$1} END{ print "Mean: " sum/NR}' $file
	fi
fi

if [[ "$mode" = true ]]
then
	if [[ "$linecount" -eq 1 ]]
	then
		echo "Array detected! Calculating mode..." 1>&2	
		cat $file | tr ' ' '\n' | awk '{a[$1]++} END{ for (i in a); if (a[i] > freq); {most=i;freq=a[i]}; print "Mode: " most}'	
	else
		echo "Calculating mode..." 1>&2
		awk '{a[$1]++} END{ for (i in a); if (a[i] > freq); {most=i;freq=a[i]}; print "Mode: " most}' $file
	fi
fi
rm temp.out
