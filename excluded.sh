#!/bin/bash

#Used to 'cross-reference' two files

# Specific flags accepted only-- if other arguments-- exit
# Print only missing items if -m
# Print only found items if -f
# Print both if -mf or -fm
# Print nothing (except for STDERR messages) if no flags
if [[ "$#" -eq 4 && "$1" == "-m" && "$1" == "-f" ]]
then
	missing=true
	found=true
	shift shift
elif [[ "$#" -eq 4 && "$1" == "-f" && "$1" == "-m" ]]
then
	missing=true
	found=true
	shift
	shift
elif [[ "$#" -eq 3 && "$1" == "-m" ]]
then
	missing=true
	found=false
	shift
elif [[ "$#" -eq 3 && "$1" == "-f" ]]
then
	missing=false
	found=true
	shift
elif [[ "$#" -eq 3 ]] && [[ "$1" == "-mf" || "$1" == "-fm" ]]
then
	missing=true
	found=true
	shift
elif [[ "$#" -eq 2 ]] && [[ "${1:0:1}" != "-" && "${2:0:1}" != "-" ]]
then
	missing=false
	found=false
else
	echo "USAGE: $(basename $0) [-mf] <keyword file> <searched file>"
	echo -e "OPTIONS:\n\t- [-m] print missing items\n\t- [-f] print found items"
	exit 1
fi

keywords=$1
file=$2
IFS=$'\n'
num_missing=0
num_found=0
for keyword in $(cat $keywords)
do
	# If keyword is a FASTA header, extract sequence name, else take first word as sequence name
	if [ "${keyword:0:1}" == ">" ]
	then
		protein=$(echo $keyword | awk -F ">" '{print $2}' | awk '{print $1}')
	else
		protein=$(echo $keyword | awk '{print $1}')
	fi
	count=$(grep -c $protein $file)
	if [ "$count" -eq 0 ]
	then
		num_missing=$((num_missing+1))
		if [[ "$missing" = true ]]
		then
			echo $protein
		fi
	else
		num_found=$((num_found+1))
		if [[ "$found" = true ]]
		then
			echo $protein
		fi
	fi
done

total=$((num_missing + num_found))

# Print messages to STDERR so that STDOUT can be redirected/piped
if [[ "$missing" = true && "$found" = false ]] 
then
	echo "${num_missing}/${total} are missing." 1>&2
fi
if [[ "$found" = true  && "$missing" = false ]]
then
	echo "${num_found}/${total} have been found." 1>&2
fi
if [[ "$missing" = true && "$found" = true ]]
then
	echo "All items have been found." 1>&2
fi
if [[ "$missing" = false && "$found" = false ]]
then
	echo "None of the items have been found." 1>&2
fi
