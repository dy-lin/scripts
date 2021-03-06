#!/bin/bash
set -eu -o pipefail
# Find identical sequences when sequence IDs are different
partial=false
gethelp=false
delete=false
merge=false
inplace=false
verbose=false
PROGRAM=$(basename $0)
while getopts :dhimpv opt
do
	case $opt in
		d) delete=true;;
		h) gethelp=true;;
		i) inplace=true;delete=true;;
		m) merge=true;delete=true;;
		p) partial=true;;
		v) verbose=true;;
		\?) echo "$PROGRAM: invalid option: $opt" >&2; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ ( "$#" -ne 1 && "$#" -ne 2 ) || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM [-dhimpv] <FASTA file>" >&2
	echo "DESCRIPTION: Checks for duplicates in sequences (not sequence IDs). If given one FASTA file, checks for duplicate sequences within itself. If given two FASTA files, checks for duplicates among the two." >&2
	echo -e "OPTIONS:\n\t-d\tDelete the matches (changes documented in a new file)\n\t-h\tShow help menu\n\t-i\tDelete matches in place (implies -d)\n\t-m\tMerge the files (implies -d)\n\t-p\tInclude partial matches\n\t-v\tVerbose (print matches to STDOUT)" >&2
	exit 1
fi

if [[ "$merge" = true && "$#" -ne 2 ]]
then
	echo -e "Two FASTA files must given in order to merge two files into one.\n" >&2
	echo "USAGE: $PROGRAM [-dhimpv] <FASTA file>" >&2
	echo "DESCRIPTION: Checks for duplicates in sequences (not sequence IDs). If given one FASTA file, checks for duplicate sequences within itself. If given two FASTA files, checks for duplicates among the two." >&2
	echo -e "OPTIONS:\n\t-d\tDelete the matches (changes documented in a new file)\n\t-h\tShow help menu\n\t-i\tDelete matches in place (implies -d)\n\t-m\tMerge the two files (implies -d)\n\t-p\tInclude partial matches\n\t-v\tVerbose (print matches to STDOUT)" >&2
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
		if [[ "$verbose" = true ]]
		then
			echo "------------------------------------------------------------------"
			echo "Identical and Partial Match sequences in $(basename $file1):"
		fi
	else
		if [[ "$verbose" = true ]]
		then
			echo -e "------------------------------------------------------------------"
			echo "Identical sequences in $(basename $file1):"
		fi
	fi
		if [[ "$verbose" = true ]]
		then
			echo -e "------------------------------------------------------------------"
		fi
fi
count=0
dir1=$(dirname $file1)
basefile1=$(basename $file1)
if [[ "$delete" = true ]]
then
	filename1=${basefile1%.*}
	ext1=${basefile1##*.}
	if [[ "$partial" = true ]]
	then
		output=${dir1}/${filename1}.partial.unique.${ext1}
	else
		output1=${dir1}/${filename1}.unique.${ext1}
	fi
	cp $file1 $output1
	if [[ "$partial" = false && -e "${dir1}/${filename1}.deletions.${ext1}" ]]
	then
		rm "${dir1}/${filename1}.deletions.${ext1}"
	fi
	if [[ "$partial" = true && -e "${dir1}/${filename1}.partial.deletions.${ext1}" ]]
	then
		rm "${dir1}/${filename1}.partial.deletions.${ext1}"
	fi
	if [ "$#" -eq 2 ]
	then
		dir2=$(dirname $file2)
		basefile2=$(basename $file2)
		filename2=${basefile2%.*}
		ext2=${basefile2##*.}
		if [[ "$partial" = true ]]
		then
			output2=${dir2}/${filename2}.partial.unique.${ext2}
		else
			output2=${dir2}/${filename2}.unique.${ext2}
		fi
		cp $file2 $output2
	fi
fi
# If working with two files, print duplicate sequences' seqIDs from each file, adding a tab before each seqID
# If working with one file, print all duplicate sequences' seqIDs
IFS=$'\n'
delcount=0
keepcount=0
for line in $(awk '!/^>/ {print}' $file1 | sort -u)
do
	if [[ "$partial" = false ]]
	then
		if [[ "$#" -eq 2 ]]
		then
			if [[ "$(grep -cw "$line" $file2)" -gt 0 ]]
			then
				if [[ "$verbose" = true ]]
				then
					echo "Sequence(s) from $(basename $file1):"
					grep -wB1 "$line" $file1 | awk '/^>/ {print}' | sed 's/^>/\t>/'
				fi
					hits=$(grep -wc "$line" $file1)
					count=$((count+hits))
				if [[ "$delete" = true ]]
				then
					if [[ "$(grep -cw "$line" $file1)" -gt 1 ]]
					then
						for sequence in $(grep -wB1 "$line" $file1 | awk '/^>/ {print}' | tail -n +2)
						do
							# escape backslashes, forward slashes, and square brackets
							sequence=$(echo $sequence | sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
							sed -i "/$sequence/,+1 d" $output1
						done
					fi
				fi
				
				if [[ "$verbose" = true ]]
				then
					echo -e "\nIdentical sequence(s) from $(basename $file2):"
					grep -wB1 "$line" $file2 | awk '/^>/ {print}' | sed 's/^>/\t>/'
				fi
				if [[ "$delete" = true ]]
				then
					for sequence in $(grep -wB1 "$line" $file2 | awk '/^>/ {print}')
					do
						sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
						sed -i "/$sequence/,+1 d" $output2
					done
				fi
				if [[ "$verbose" = true ]]
				then
					echo -e "\n------------------------------------------------------------------\n"
				fi
			fi
		# if $# = 1
		else
			if [[ "$(grep -wc "$line" $file2)" -gt 1 ]]
			then
				if [[ "$verbose" = true ]]
				then
					grep -wB1 "$line" $file1 | awk '/^>/ {print}'
				fi
				hits=$(grep -wc "$line" $file1)
				count=$((count + hits))
				if [[ "$delete" = true ]]
				then
					removed=false
					# Put all matches to sequence in deletions, but needs to keep one, so the kept one must be removed from the deletions
					# the kept one will either be unique if removed=true, or if removed=false, just remove the last two lines of the file
					grep -wB1 --no-group-separator "$line" $file1 >> ${dir1}/${filename1}.deletions.${ext1}
					for sequence in $(grep -wB1 "$line" $file1 | awk '/^>/ {print}')
					do
						sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
						if [[ "$(grep -cw "$sequence" $file1)" -eq 1 ]]
						then
							if [[ "$removed" = false ]]
							then
								removed=true
								kept="$sequence"
								continue
							fi
						fi
						sed -i "/$sequence/,+1 d" $output1
					done
					# Remove the one kept from the deletion file if one was kept
					if [[ "$removed" = true ]]
					then
						sed -i "/$kept/,+1 d" ${dir1}/${filename1}.deletions.${ext1}
					fi
					# If at the end of the loop, removed is still false, then no deletion occurred....
					# Remove all but one-- but sed will replace all -- add last one back (since none of sequence IDs are unique)
					if [[ "$removed" = false ]]
					then
						for sequence in $(grep -wB1 "$line" $file1 | awk '/^>/ {print}')
						do
							sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
							# Remove all matched sequences
							sed -i "/$sequence/,+1 d" $output1

						done
						# added back the last two as to keep one of each duplicate
						grep -wB1 --no-group-separator "$line" $file1 | tail -n 2 >> $output1
						# all were added to deletion list, but one was added back, need to remove that one from the deletion list, remove last two lines of deletion list
						head -n -2 ${dir1}/${filename1}.deletions.${ext1} > temp.out
						mv temp.out ${dir1}/${filename1}.deletions.${ext1}
					fi
						
					delcount=$((delcount+hits-1))
					keepcount=$((keepcount+1))
				fi
				if [[ "$verbose" = true ]]
				then
					echo -e "\n------------------------------------------------------------------\n"
				fi	
			fi
		fi
	# if partial = true
	else
		if [[ "$#" -eq 2 ]]
		then
			if [[ "$(grep -c "$line" $file2)" -gt 0 ]]
			then
				if [[ "$verbose" = true ]]
				then
					echo "Sequence(s) from $(basename $file1):"
					grep -B1 "$line" $file1 | awk '/^>/ {print}' | sed 's/^>/\t>/'	
				fi
					hits=$(grep -c "$line" $file1)
					count=$((count+hits))
				if [[ "$delete" = true ]]
				then
					if [[ "$(grep -c "$line" $file1)" -gt 1 ]]
					then
						for sequence in $(grep -B1 "$line" $file1 | awk '/^>/ {print}' | tail -n +2)
						do
							sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
							sed -i "/$sequence/,+1 d" $output1
						done
					fi
				fi

				if [[ "$verbose" = true ]]
				then
					echo -e "\nIdentical and Partial Match sequence(s) from $(basename $file2):"
					grep -B1 "$line" $file2 | awk '/^>/ {print}' | sed 's/^>/\t>/'
				fi
				if [[ "$delete" = true ]]
				then
					for sequence in $(grep -B1 "$line" $file2 | awk '/^>/ {print}')
					do
						sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g') 
						sed -i "/$sequence/,+1 d" $output2
					done
				fi
				if [[ "$verbose" = true ]]
				then
					echo -e "\n------------------------------------------------------------------\n"
				fi
			fi
		else
			if [[ "$(grep -c "$line" $file2)" -gt 1 ]]
			then
				if [[ "$verbose" = true ]]
				then
					grep -B1 "$line" $file1 | awk '/^>/ {print}'
				fi
				hits=$(grep -c "$line" $file1)
				count=$((count + hits))
				if [[ "$delete" = true ]]
				then
					removed=false
					grep -B1 --no-group-separator "$line" $file1 >> ${dir1}/${filename1}.deletions.${ext1}
					for sequence in $(grep -B1 "$line" $output1 | awk '/^>/ {print}' | tail -n +2)
					do
						sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
						if [[ "$(grep -c "$sequence" $file1)" -eq 1 ]]
						then
							if [[ "$removed" = false ]]
							then
								removed=true
								kept="$sequence"
								continue
							fi
						fi
						sed -i "/$sequence/,+1 d" $output1
					
					done
					if [[ "$removed" = true ]]
					then
						sed -i "/$kept/,+1 d" ${dir1}/${filename1}.deletions.${ext1}
					fi
					if [[ "$removed" = false ]]
					then
						for sequence in $(grep -B1 "$line" $file1 | awk '/^>/ {print}')
						do
							sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
							sed -i "/$sequence/,+1 d" $output1
						done
						grep -B1 --no-group-separator "$line" $file1 | tail -n 2 >> $output1
						head -n -2 ${dir1}/${filename1}.deletions.${ext1} > temp.out
						mv temp.out ${dir1}/${filename1}.deletions.${ext1}
					fi
					delcount=$((delcount+hits-1))
					keepcount=$((keepcount+1))
				fi
				if [[ "$verbose" = true ]]
				then
					echo -e "\n------------------------------------------------------------------\n"
				fi
			fi
		fi
	fi
done

# If Partial is true, then sequences from file2 must be checked against sequences in file1
revcount=0
if [[ "$partial" = true && "$#" -eq 2 ]]
then
	for line in $(awk '!/^>/ {print}' $file2 | sort -u)
	do
		duplicate=false
		# If sequence from file2 is found in file1
		if [[ "$(grep -c "$line" $file1)" -gt 0 ]]
		then
			if [[ "$verbose" = true ]]
			then
				for item in $(grep "$line" $file2)
				do
					if [[ "$(grep -wc "$line" $file1)" -gt 0 ]]
					then
						duplicate=true
						break
					fi
				done

				# Print those sequence names from file2, only if they don't match entirely in file1 (no reciprocal duplicates for identical sequences)
				if [[ "$duplicate" = false ]]
				then
					echo "Sequence(s) from $file2:"
					grep -B1 "$line" $file2 | awk '/^>/ {print}' | sed 's/^>/\t>/'
				fi
			fi
			if [[ "$delete" = true ]]
			then
				# if delete = true, remove all but the first sequence from file2
				for sequence in $(grep -B1 "$line" $output2 | awk '/^>/ {print}' | tail -n +2)
				do
					sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g'| sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
					sed -i "/$sequence/,+1 d" $output2
				done
			fi
			if [[ "$verbose" = true ]]
			then
				# Print those sequence names from file 1
				if [[ "$duplicate" = false ]]
				then
					echo "Identical and Partial Match sequence(s) from $file1:"
					grep -B1 "$line" $file1 | awk '/^>/ {print}' | sed 's/^>/\t>/'
				fi
			fi
			hits=$(grep -c "$line" $file2)
			revcount=$((revcount+hits))
			if [[ "$delete" = true ]]
			then
				# if delete = true, remove all of those from file1
				for sequence in $(grep -B1 "$line" $output1 | awk '/^>/ {print}')
				do
					sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
					sed -i "/$sequence/,+1 d" $output1
				done
			fi
			if [[ "$verbose" = true ]]
			then
				if [[ "$duplicate" = false ]]
				then
					echo -e "\n------------------------------------------------------------------\n"
				fi
			fi
		fi
	done
fi		

if [[ "$verbose" = true && "$#" -eq 1 && "$count" -eq 0 ]]
then
	echo -e "\tNONE"
fi

# Check for duplicate sequences within file 1 and remove them
if [[ "$delete" = true ]]
then
	for line in $(awk '!/^>/ {print}' $output1 | sort -u)
	do
		if [[ "$partial" = true ]]
		then
			if [[ "$(grep -c "$line" $output1)" -gt 1 ]]
			then
				removed=false
				if [[ "$#" -eq 1 ]]
				then
					grep -B1 --no-group-separator "$line" $file1 >> ${dir1}/${filename1}.deletions.${ext1}
				fi
				for sequence in $(grep -B1 "$line" $output1 | awk '/^>/ {print}')
				do
					sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g') 
					if [[ "$(grep -c "$sequence" $file1)" -eq 1 ]]
					then
						if [[ "$removed" = false ]]
						then
							removed=true
							kept="$sequence"
							continue
						fi
					fi
					sed -i "/$sequence/,+1 d" $output1
					delcount=$((delcount+1))
				done
				if [[ "$removed" = true && "$#" -eq 1 ]]
				then
					sed -i "/$kept/,+1 d" ${dir1}/${filename1}.deletions.${ext1}
				fi
				if [[ "$removed" = false ]]
				then
					for sequence in $(grep -B1 "$line" $file1 | awk '/^>/ {print}')
					do
						sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
						sed -i "/$sequence/,+1 d" $output1
					done
					grep -B1 --no-group-separator "$line" $file1 | tail -n 2 >> $output1
					if [[ "$#" -eq 1 ]]
					then
						head -n -2 ${dir1}/${filename1}.deletions.${ext1} > temp.out
						mv temp.out ${dir1}/${filename1}.deletions.${ext1}
					fi
				fi
				keepcount=$((keepcount+1))
			fi
		else
			if [[ "$(grep -wc "$line" $output1)" -gt 1 ]]
			then
				removed=false
				if [[ "$#" -eq 1 ]]
				then
					grep -wB1 --no-group-separator "$line" $file1 >> ${dir1}/${filename1}.deletions.${ext1}
				fi
				for sequence in $(grep -wB1 "$line" $output1 | awk '/^>/ {print}' | tail -n +2)
				do

					sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
					if [[ "$(grep -cw "$sequence" $file1)" -eq 1 ]]
					then
						if [[ "$removed" = false ]]
						then
							removed=true
							kept="$sequence"
							continue
						fi
					fi
					sed -i "/$sequence/,+1 d" $output1
					delcount=$((delcount+1))
				done
				if [[ "$removed" = true && "$#" -eq 1 ]]
				then
					sed -i "/$kept/,+1 d" ${dir1}/${filename1}.deletions.${ext1}
				fi
				if [[ "$removed" = false ]]
				then
					for sequence in $(grep -wB1 "$line" $file1 | awk '/^>/ {print}')
					do
						sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
						sed -i "/$sequence/,+1 d" $output1
					done
					grep -wB1 --no-group-separator "$line" $file1 | tail -n 2 >> $output1
					if [[ "$#" -eq 1 ]]
					then
						head -n -2 ${dir1}/${filename1}.deletions.${ext1} > temp.out
						mv temp.out ${dir1}/${filename1}.deletions.${ext1}
					fi
				fi

			fi
		fi
	done
fi
if [[ "$verbose" = false && "$#" -eq 1 ]]
then
	ogcount=$(awk '/^>/ {print}' $file1 | wc -l)
	if [[ "$partial" = true ]] 
	then
		echo "There are $count sequence(s) that are entirely or partially in $basefile1."
		#echo "Number of identical and partial match sequence(s) in $basefile1: $count"
		if [[ "$count" -ne 0 ]]
		then
			if [[ "$delete" = true ]]
			then
				echo "Of those identical and partial match sequence(s), $keepcount were kept and its $delcount duplicates were removed."
				echo "After removal, there are now $((ogcount-delcount)) unique sequences (originally a total of $ogcount sequences)."
				if [[ "$inplace" = true ]]
				then
					echo "Changes made directly to $basefile1." >&2
					if [[ -e "${dir1}/${filename1}.partial.deletions.${ext1}" ]]
					then
						echo "Deleted sequences saved in ${filename1}.partial.deletions.${ext1}" >&2
					fi
				else 
					echo "Writing modified $basefile1 to $(basename $output1)." >&2
					if [[ -e "${dir1}/${filename1}.partial.deletions.${ext1}" ]]
					then
						echo "Deleted sequences saved in ${filename1}.partial.deletions.${ext1}" >&2
					fi
				 fi
			fi
		fi
	else
		echo "There are $count identical sequence(s) within $basefile1."
	#	echo "Number of identical sequence(s) in $basefile: $count"
		if [[ "$count" -ne 0 ]]
		then
			if [[ "$delete" = true ]]
			then
				echo "Of those identical sequence(s), $keepcount were kept, and its $delcount duplicates were removed."
				echo "After removal, there are now $((ogcount-delcount)) unique sequences (originally a total of $ogcount sequences)."
				if [[ "$inplace" = true ]]
				then
					echo "Changes made directly to $basefile1." >&2
					if [[ -e "${dir1}/${filename1}.deletions.${ext1}" ]]
					then
						echo "Deleted sequences saved in ${filename1}.deletions.${ext1}" >&2
					fi	
				else
					echo "Writing modified $basefile1 to $(basename $output1)." >&2
					if [[ -e "${dir1}/${filename1}.deletions.${ext1}" ]]
					then
						echo "Deleted sequences saved in ${filename1}.deletions.${ext1}" >&2
					fi	
				fi
			fi
		fi
	fi
fi





# Check for duplicate sequences in file 2 and remove them
# changed output2 to file2 except for in sed
if [[ "$delete" = true && "$#" -eq 2 ]]
then
	for line in $(awk '!/^>/ {print}' $file2 | sort -u)
	do
		if [[ "$partial" = true ]]
		then
			if [[ "$(grep -c "$line" $file2)" -gt 1 ]]
			then
				for sequence in $(grep -B1 "$line" $output2 | awk '/^>/ {print}' | tail -n +2)
				do
					sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
					sed -i "/$sequence/,+1 d" $output2
				done
			fi
		else
			if [[ "$(grep -wc "$line" $file2)" -gt 1 ]]
			then
				for sequence in $(grep -wB1 "$line" $output2 | awk '/^>/ {print}' | tail -n +2)
				do
					sequence=$(echo $sequence |sed 's~\\~\\\\\\~g' | sed 's~\/~\\/~g' | sed 's~\[~\\[~g' | sed 's~\]~\\]~g')
					sed -i "/$sequence/,+1 d" $output2
				done
			fi
		fi
	done
fi

# Find partial and identical counts within file1 and file2 individually
dupcount1=0
dupcount2=0
if [[ "$partial" = true ]]
then
	for line in $(awk '!/^>/ {print}' $file1 | sort -u)
	do
		hits=$(grep -c $line $file1)
		if [[ "$hits" -gt 1 ]]
		then
			dupcount1=$((dupcount1+hits))
		fi
	done
	
	for line in $(awk '!/^>/ {print}' $file2 | sort -u)
	do
		hits=$(grep -c $line $file2)
		if [[ "$hits" -gt 1 ]]
		then
			dupcount2=$((dupcount2+hits))
		fi
	done
fi

if [[ "$verbose" = false && "$#" -eq 2 ]]
then
	if [[ "$partial" = false ]]
	then
		echo "There are $count identical sequence(s) from $basefile1 in $basefile2."
		#	echo "Number of identical sequence(s) from $basefile1 in $basefile2: $count"
		num1=$(awk '!/^>/ {print}' $file1 | sort | uniq -c | sed 's/^\t*[[:space:]]*//' | awk '{if($1>1) print $1}')
		num2=$(awk '!/^>/ {print}' $file2 | sort | uniq -c | sed 's/^\t*[[:space:]]*//' | awk '{if($1>1) print $1}')
		sum1=0
		for num in $num1
		do
			sum1=$((sum1+num))
		done
		sum2=0
		for num in $num2
		do
			sum2=$((sum2+num))
		done
		echo "There are $sum1 identical sequence(s) within $basefile1."
		echo "There are $sum2 identical sequence(s) within $basefile2."
		#		echo "Number of identical sequence(s) within $basefile1: $sum1"
		#		echo "Number of identical sequence(s) within $basefile2: $sum2"
	else
		echo "There are $count sequence(s) in $basefile1 that are entirely or partially in $basefile2."
		#		echo "Number of identical and partial match sequence(s) from $basefile1 in $basefile2: $count"
		#		echo -e "Number of identical and partial match sequence(s) from $basefile2 in $basefile1: $revcount\n"
		echo -e "There are $revcount sequence(s) in $basefile2 that are entirely or partially in $basefile1.\n"
		num1=$(awk '!/^>/ {print}' $file1 | sort | uniq -c | sed 's/^\t*[[:space:]]*//' | awk '{if($1>1) print $1}')
		num2=$(awk '!/^>/ {print}' $file2 | sort | uniq -c | sed 's/^\t*[[:space:]]*//' | awk '{if($1>1) print $1}')
		sum1=0
		for num in $num1
		do
			sum1=$((sum1+num))
		done
		sum2=0
		for num in $num2
		do
			sum2=$((sum2+num))
		done
		if [[ "$dupcount1" -gt "$sum1" ]]
		then
			partialcount1=$((dupcount1-sum1))
		else
			partialcount1=0
		fi
		if [[ "$dupcount2" -gt "$sum2" ]]
		then
			partialcount2=$((dupcount2-sum2))
		else
			partialcount2=0
		fi
		echo "There are $sum1 identical sequence(s) within $basefile1."
		#		echo "Number of identical sequence(s) within $basefile1: $sum1"
		echo -e "There are $partialcount1 partially matched sequence(s) within $basefile1.\n"
		#		echo -e "Number of partial match sequence(s) within $basefile1: $partialcount1\n"
		echo "There are $sum2 identical sequence(s) within $basefile2."
		#		echo "Number of identical sequence(s) within $basefile2: $sum2"
		echo -e "There are $partialcount2 partially matched sequence(s) within $basefile2.\n"
		#		echo "Number of partial match sequence(s) within $basefile2: $partialcount2"
	fi
fi


if [[ "$inplace" = true ]]
then
	if [[ "$#" -eq 2 ]]
	then
		mv $output1 $file1
		mv $output2 $file2
	else
		mv $output1 $file1
	fi
fi

if [[ "$merge" = true && "$#" -eq 2 ]]
then
	if [[ "$inplace" = true ]] 
	then
		echo "Changes made directly to $basefile1 and $basefile2." >&2
		if [[ "$partial" = true ]]
		then
			cat $file1 $file2 > ${dir1}/${filename1}.${filename2}.partial.merged.${ext1}
			echo "Merged $basefile1 and $basefile2 into ${filename1}.${filename2}.partial.merged.${ext1}." >&2
		else
			cat $file1 $file2 > ${dir1}/${filename1}.${filename2}.merged.${ext1}
			echo "Merged $basefile1 and $basefile2 into ${filename1}.${filename2}.merged.${ext1}." >&2
		fi	
	else
		echo "Writing modified $basefile1 to $(basename $output1)." >&2
		echo "Writing modified $basefile2 to $(basename $output2)." >&2
		if [[ "$partial" = true ]]
		then
			cat $output1 $output2 > ${dir1}/${filename1}.${filename2}.partial.merged.${ext1}
			echo "Merged $(basename $output1) and $(basename $output2) into ${filename1}.${filename2}.partial.merged.${ext1}." >&2
		else
			cat $output1 $output2 > ${dir1}/${filename1}.${filename2}.merged.${ext1}
			echo "Merged $(basename $output1) and $(basename $output2) into ${filename1}.${filename2}.merged.${ext1}." >&2
		fi
	fi
fi

if [[ "$merge" = false && "$#" -eq 2 && "$delete" = true ]]
then
	if [[ "$inplace" = true ]]
	then
		echo "Changes made directly to $basefile1 and $basefile2." >&2
	else
		echo "Writing modified $basefile1 to $(basename $output1)." >&2
		echo "Writing modified $basefile2 to $(basename $output2)." >&2
	fi
fi

if [[ "$count" -eq 0 && "$#" -eq 1 && "$delete" = true && -e "$output1" ]]
then
	rm $output1
fi
