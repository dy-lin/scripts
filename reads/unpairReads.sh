#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)
gethelp=false
verbose=false
while getopts :hv opt
do
	case $opt in
		h) gethelp=true;;
		v) verbose=true;;
		\?) echo "$PROGRAM: invalid option: $OPTARG" 1>&2; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ "$#" -ne 1 || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM <list of read pairs>" 1>&2
	echo "DESCRIPTION: Takes a list of reads and splits the read pairs into two separate files." 1>&2
	echo -e "OPTIONS:\n\t-h\tLoad help menu\n\t-v\tVerbose (prints contents of file on screen)" 1>&2
	exit 1
fi

file=$1

if [[ -e "reads1.in" ]]
then
	rm reads1.in
fi

if [[ -e "reads2.in" ]]
then
	rm reads2.in
fi
# Copy original file to a temp so as to not modify original file
cp $file temp.txt
file=temp.txt

# if all reads are on the same line
if [[ "$(head -n 1 $file | awk '{print NF}')" -gt 2 ]]
then
	echo "$(head -n 1 $file | awk '{print NF}') reads on the same line detected..." 1>&2
	pairReads.sh $file > temp.out
	cp temp.out $file
	rm temp.out
fi

# if all read pairs are individually on separate line
if [[ "$(head -n 1 $file | awk '{print NF}')" -eq 1 ]]
then
	echo "$(wc -l $file | awk '{print $1}') lines with one read per line detected..." 1>&2
	pairReads.sh $file > temp.out
	cp temp.out $file
	rm temp.out
fi

# if there are pairs on each line
if [[ "$(head -n 1 $file | awk '{print NF}')" -eq 2 ]]
then
	echo "$(wc -l $file | awk '{print $1}') lines with one read pair (2 reads) per line detected..." 1>&2
	while IFS=' ' read -r reads1 reads2
	do
		echo $reads1 >> reads1.in
		echo $reads2 >> reads2.in
	done < $file
fi

if [[ "$verbose" = true ]]
then
	echo -e "\nreads1.in:" 1>&2
	cat reads1.in
	echo -e "\nreads2.in:" 1>&2
	cat reads2.in
else
	echo "...DONE!" 1>&2
fi
if [[ -e "$file" ]]
then
	rm "$file"
fi
