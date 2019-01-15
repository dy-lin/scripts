#!/bin/bash
set -eu -o pipefail
# Find identical sequences when sequence IDs are different
partial=false
gethelp=false
delete=false
PROGRAM=$(basename $0)
while getopts :hpd opt
do
	case $opt in
		d) delete=true;;
		h) gethelp=true;;
		p) partial=true;;
		\?) echo "$PROGRAM: invalid option: $opt" >&2; exit 1;;
	esac
done

shift $((OPTIND-1))

if [[ ( "$#" -ne 1 && "$#" -ne 2 ) || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM [-dhp] <FASTA file>" >&2
	echo "DESCRIPTION: Checks for duplicates in sequences (not sequence IDs). If given one FASTA file, checks for duplicate sequences within itself. If given two FASTA files, checks for duplicates among the two." >&2
	echo -e "OPTIONS:\n\t-d\tDelete the matches\n\t-h\tShow help menu\n\t-p\tInclude partial matches" >&2
	exit 1
fi

if [ "$#" -eq 2 ]
then
	file1=$1
	file2=$2
fi

if [ "$#" -eq 1 ]
then
	file1=$1
	file2=$1
	if [[ "$partial" = true ]]
	then
		echo "----------------------------------------------"
		echo "Identical and Partial Match sequences in $(basename $file1):"
	else
		echo "----------------------------------------------"
		echo "Identical sequences in $(basename $file1):"
	fi
		echo -e "----------------------------------------------\n"
fi
if [[ "$delete" = true ]]
then
	filename1=${file1%.*}
	ext1=${file1##*.}
	output1=${filename1}.duplicatesRemoved.${ext1}
	cp $file1 $output1
	if [ "$#" -eq 2 ]
	then
		filename2=${file2%.*}
		ext2=${file2##*.}
		output2=${filename2}.duplicatesRemoved.${ext2}
		cp $file2 $output2
	fi
fi
# If working with two files, print duplicate sequences' seqIDs from each file, adding a tab before each seqID
# If working with one file, print all duplicate sequences' seqIDs

for line in $(awk '!/^>/ {print}' $file1 | sort -u)
do
	if [[ "$partial" = false ]]
	then
		if [[ "$#" -eq 2 ]]
		then
			if [[ "$(grep -cw "$line" $file2)" -gt 0 ]]
			then
				echo "Sequence(s) from $(basename $file1):"
				grep -wB1 "$line" $file1 | awk '/^>/ {print}' | sed 's/^>/\t>/'
				
				if [[ "$delete" = true ]]
				then
					if [[ "$(grep -cw "$line" $file1)" -gt 1 ]]
					then
						for sequence in $(grep -wB1 "$line" $file1 | awk '/^>/ {print}' | tail -n +2)
						do
							sed -i "/$sequence/,+1 d" $output1
						done
					fi
				fi
				
				echo -e "\nIdentical sequence(s) from $(basename $file2):"
				grep -wB1 "$line" $file2 | awk '/^>/ {print}' | sed 's/^>/\t>/'
				
				if [[ "$delete" = true ]]
				then
					for sequence in $(grep -wB1 "$line" $file2 | awk '/^>/ {print}')
					do
						sed -i "/$sequence/,+1 d" $output2
					done
				fi
				
				echo -e "\n----------------------------------------------\n"
			fi
		else
			if [[ "$(grep -wc "$line" $file2)" -gt 1 ]]
			then
				grep -wB1 "$line" $file1 | awk '/^>/ {print}'
				if [[ "$delete" = true ]]
				then
					for sequence in $(grep -wB1 "$line" $file1 | awk '/^>/ {print}')
					do
						sed -i "/$sequence/,+1 d" $output1
					done
				fi
				echo -e "\n----------------------------------------------\n"
			fi	
		fi
	# if partial = false
	else
		if [[ "$#" -eq 2 ]]
		then
			if [[ "$(grep -c "$line" $file2)" -gt 0 ]]
			then
				echo "Sequence(s) from $(basename $file1):"
				grep -B1 "$line" $file1 | awk '/^>/ {print}' | sed 's/^>/\t>/'
				if [[ "$delete" = true ]]
				then
					if [[ "$(grep -c "$line" $file1)" -gt 1 ]]
					then
						for sequence in $(grep -B1 "$line" $file1 | awk '/^>/ {print}' | tail -n +2)
						do
							sed -i "/$sequence/,+1 d" $output1
						done
					fi
				fi
				echo -e "\nIdentical and Partial Match sequence(s) from $(basename $file2):"
				grep -B1 "$line" $file2 | awk '/^>/ {print}' | sed 's/^>/\t>/'
				if [[ "$delete" = true ]]
				then
					for sequence in $(grep -B1 "$line" $file2 | awk '/^>/ {print}')
					do
						sed -i "/$sequence/,+1 d" $output2
					done
				fi
				echo -e "\n----------------------------------------------\n"
			fi
		else
			if [[ "$(grep -c "$line" $file2)" -gt 1 ]]
			then
				grep -B1 "$line" $file1 | awk '/^>/ {print}'
				if [[ "$delete" = true ]]
				then
					for sequence in $(grep -B1 "$line" $file1 | awk '/^>/ {print}')
					do
						sed -i "/$sequence/,+1 d" $output1
					done
				fi
				echo -e "\n----------------------------------------------\n"
			fi
		fi
	fi
done
if [[ "$delete" = true && "$#" -eq 2 ]]
then
	cat $output1 $output2 > ${filename1}.${filename2}.merged.${ext1}
fi
