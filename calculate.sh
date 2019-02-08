#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)
average=false
median=false
mode=false
gethelp=false
while getopts :amoh opt
do
	case $opt in
		a) average=true;;
		h) gethelp=true;;
		m) median=true;;
		o) mode=true;;
		\?) echo "$PROGRAM: invalid option: $opt" >&2; exit 1;;
	esac
done
shift $((OPTIND-1))
if [[ "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM [OPTIONS] <FILE>" 1>&2
	echo "DESCRIPTION: Takes a file with numbers, and prints the desired calculation result to STDOUT." 1>&2
	echo -e "OPTIONS: At least ONE of -a, -m or -o must be selected.\n\t-a\tmean (average) \n\t-m\tmedian\n\t-o\tmode\n\t-h\tShow help menu" 1>&2 
	exit 1
fi

pipe=false
if [ -p /dev/stdin ]
then
	pipe=true
	if [[ "$#" -ne 0 ]]
	then
		echo "USAGE: $PROGRAM [OPTIONS] <FILE>" 1>&2
		echo "DESCRIPTION: Takes a file with numbers, and prints the desired calculation result to STDOUT." 1>&2
		echo -e "OPTIONS: At least ONE of -a, -m, or -o must be selected.\n\t-a\tmean (average) \n\t-m\tmedian\n\t-o\tmode\n\t-h\tShow help menu" 1>&2
		exit 1
	fi
	cat > temp.out
	file=temp.out
else
	if [[ "$#" -ne 1 ]]
	then
		echo "USAGE: $PROGRAM [OPTIONS] <FILE>" 1>&2
		echo "DESCRIPTION: Takes a file with numbers, and prints the desired calculation result to STDOUT." 1>&2
		echo -e "OPTIONS: At least ONE of -a, -m, or -o must be selected.\n\t-a\tmean (average) \n\t-m\tmedian\n\t-o\tmode\n\t-h\tShow help menu" 1>&2
		exit 1
	fi
	if [[ "$1" == "/dev/fd/63" ]]
	then
		cp $1 temp.out
		file=temp.out
	else
		file=$1
	fi
fi
if [[ -e "$file" ]]
then
	filesize=$(ls -l $file | awk '{print $5}')
else
	if [[ "$pipe" == false ]]
	then
		echo "ERROR: File does not exist." 1>&2
	else
		echo "ERROR: STDIN is empty." 1>&2
	fi
		echo "USAGE: $PROGRAM [OPTIONS] <FILE>" 1>&2
		echo "DESCRIPTION: Takes a file with numbers, and prints the desired calculation result to STDOUT." 1>&2
		echo -e "OPTIONS: At least ONE of -a, -m, or -o must be selected.\n\t-a\tmean (average) \n\t-m\tmedian\n\t-o\tmode\n\t-h\tShow help menu" 1>&2
		if [[ -e "temp.out" ]]
		then
			rm temp.out
		fi

		exit 1
		
fi
if [[ "$filesize" -le 1 ]]
then
	if [ ! -s "$file" ]
	then
		if [[ "$pipe" == false ]]
		then
			echo "ERROR: File does not exist." 1>&2
		else
			echo "ERROR: STDIN is empty." 1>&2
		fi
	else
		if [[ "$pipe" == false ]]
		then
			echo "ERROR: File is empty." 1>&2
		else
			echo "ERROR: STDIN is empty." 1>&2
		fi
	fi
	echo "USAGE: $PROGRAM [OPTIONS] <FILE>" 1>&2
	echo "DESCRIPTION: Takes a file with numbers, and prints the desired calculation result to STDOUT." 1>&2		
	echo -e "OPTIONS: At least ONE of -a, -m or -o must be selected.\n\t-a\tmean (average) \n\t-m\tmedian\n\t-o\tmode\n\t-h\tShow help menu" 1>&2
	if [[ -e "temp.out" ]]
	then
		rm temp.out
	fi
	exit 1
fi


if [[ "$average" = false && "$median" = false && "$mode" = false ]]
then
	echo "ERROR: At least one calculation has to be selected." 1>&2
	echo "USAGE: $PROGRAM [OPTIONS] <FILE>" 1>&2
	echo "DESCRIPTION: Takes a file with numbers, and prints the mean/median/mode (depending on flags) to STDOUT." 1>&2
	echo -e "OPTIONS: At least ONE of -a, -m or -o must be selected.\n\t-a\tmean (average) \n\t-m\tmedian\n\t-o\tmode\n\t-h\tShow help menu" 1>&2
	if [[ -e "temp.out" ]]
	then
		rm temp.out
	fi
	exit 1
fi

# Check if input are numbers!
number=true
for num in $(cat $file)
do
	if [[ "$num" != ?(-)+([0-9])?(.)*([0-9]) ]]
	then
		number=false
		break
	fi
done

if [[ "$number" = false ]]
then
	echo "ERROR: Input must contain numbers only." 1>&2
	echo "USAGE: $PROGRAM [OPTIONS] <FILE>" 1>&2
	echo "DESCRIPTION: Takes a file with numbers, and prints the mean/median/mode (depending on flags) to STDOUT." 1>&2
	echo -e "OPTIONS: At least ONE of -a, -m or -o must be selected.\n\t-a\tmean (average) \n\t-m\tmedian\n\t-o\tmode\n\t-h\tShow help menu" 1>&2 
	if [[ -e "temp.out" ]]
	then
		rm temp.out
	fi
	exit 1
fi

linecount=$(wc -l $file | awk '{print $1}')
if [[ "$median" = true ]]
then
	# if numbers all on one line, e.g. an array
	if [[ "$linecount" -eq 1 ]]
	then
		echo "Array detected! Calculating median..." 1>&2
		cat $file | tr ' ' '\n' | sort -n | awk ' {a[i++]=$1} END {x=int((i+1)/2); if (x < (i+1)/2) print "Median: " (a[x-1]+a[x])/2; else print "Median: " a[x-1]}'
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
		cat $file | awk '{print}' | tr ' ' '\n' | awk '{sum+=$1} END{ if (NR !=0 ) {print "Mean: " sum/NR} else { print "There are no numbers given!"; exit 1}}'
	else
		echo "Calculating mean..." 1>&2
		awk '{sum+=$1} END{ if (NR != 0) {print "Mean: " sum/NR} else { print "There are no numbers given!"; exit 1}}' $file
	fi
		if [[ "$?" -eq 1 ]]
		then
			if [[ -e "temp.out" ]]
			then
				rm temp.out
			fi
			exit 1
		fi
fi

if [[ "$mode" = true ]]
then
	if [[ "$linecount" -eq 1 ]]
	then
		echo "Array detected! Calculating mode..." 1>&2	
		cat $file | tr ' ' '\n' | awk '{a[$1]++} END{ for ( i in a ) { if (a[i] > freq ) {most=i; freq=a[i]} } for (i in a) {for (j in a) { if (a[i] != a[j]) {same=1}}} if(same == 1) {print "Mode:", most} else {print "All values occur at the same frequency."}}'
	else
		echo "Calculating mode..." 1>&2
		awk '{a[$1]++} END{ for ( i in a ) { if (a[i] > freq ) {most=i; freq=a[i]} } for (i in a) {for (j in a) { if (a[i] != a[j]) {same=1}}} if(same == 1) {print "Mode:", most} else {print "All values occur at the same frequency."; exit 1}}' $file
	fi
		if [[ "$?" -eq 1 ]]
		then
			if [[ -e "temp.out" ]]
			then
				rm temp.out
			fi
			exit 1
		fi
	
fi

if [[ -e "temp.out" ]]
then
	rm temp.out
fi
