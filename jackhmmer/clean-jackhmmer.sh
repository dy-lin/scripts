#!/bin/bash 
set -euo pipefail
force=false 
gethelp=false 
PROGRAM=$(basename $0)
while getopts :hf opt
do
	case $opt in
		h) gethelp=true;;
		f) force=true;;
		\?) echo "ERROR: $PROGRAM: invalid option: $opt" 1>&2; exit 1;;
	esac
done

shift $((OPTIND-1))
if [[ "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM [directory]" 1>&2
	echo "DESCRIPTION: Removes JACKHMMER and BLAST output files in the specified directory. If no directory is given, the current working directory is used."  1>&2
	echo -e "OPTIONS:\n\t-h\tShow help menu\n\t-f\tForce (Default: asks for confirmation before removing any files)" 1>&2
	exit 1
fi
if [[ "$#" -eq 0 ]]
then
	dir=$(pwd)
	if [[ "$force" = false ]]
	then
		while read -p "Cleaning current directory --  Are you sure? (Y/N): " answer
		do
			case $answer in
				Y) break;;
				y) break;;
				Yes) break;;
				yes) break;;
				N) echo "Aborted." 1>&2; exit 1;;
				n) echo "Aborted." 1>&2; exit 1;;
				no) echo "Aborted." 1>&2; exit 1;;	
				No) echo "Aborted." 1>&2; exit 1;;
				\?) echo "Invalid answer." 1>&2;;
			esac
		done
	fi
	for i in ${dir}/guide-*
	do
		if [[ -e "$i" ]]
		then
			rm $i
		fi
	done

	for i in ${dir}/jackhmmer_bs*
	do
		if [[ -e "$i" ]]
		then
			rm $i
		fi
	done
else
	for path in $@
	do
		dir=$(echo $path | sed 's/\/$//')
		if [[ ! -d "$dir" ]]
		then
			echo "$dir does not exist." 1>&2
			continue
		fi
		if [[ "$force" = false ]]
		then
			while read -p "Cleaning ${dir} --  Are you sure? (Y/N): " answer
			do
				case $answer in
					Y) break;;
					y) break;;
					Yes) break;;
					yes) break;;
					N) echo "Aborted." 1>&2; exit 1;;
					n) echo "Aborted." 1>&2; exit 1;;
					no) echo "Aborted." 1>&2; exit 1;;	
					No) echo "Aborted." 1>&2; exit 1;;
					\?) echo "Invalid answer." 1>&2;;
				esac
			done
		fi
		for i in ${dir}/guide-*
		do
			if [[ -e "$i" ]]
			then
				rm $i
			fi
		done

		for i in ${dir}/jackhmmer_bs*
		do
			if [[ -e "$i" ]]
			then
				rm $i
			fi
		done
	done
fi
