#!/bin/bash
set -euo pipefail
PROGRAM=$(basename $0)

function show_help() {
# check if symlinks need to be absolute paths for second argument?
	echo "USAGE #1: $PROGRAM <Symlinked item> <Symlink location/name>" 1>&2
	echo "USAGE #2: $PROGRAM <Symlink list>" 1>&2
	echo "DESCRIPTION: Takes the file or directory to link and links them to the specified directory. If a single argument is given, it will be treated as a list, where each row will be executed like USAGE #1." 1>&2
	echo -e "OPTIONS:\n\t-f\tforce\n\t-h\tshow help menu\n\t-t\ttesting (print instead of symlink)"
	exit 1
}

if [[ "$HOSTNAME" == dlin02* ]]
then
	readlink=$(which greadlink)
else
	readlink=$(which readlink)
fi

testing=false
force=false
while getopts :fht opt
do
	case $opt in
		f) force=true;;
		h) show_help;;
		t) testing=true;;
		?) echo "ERROR: $PROGRAM: Invalid option $opt" 1>&2; show_help; exit 1;;
	esac
done

shift $((OPTIND-1))

if [[ "$#" -eq 0 || "$#" -gt 2 ]]
then
	show_help
fi

if [[ "$#" -eq 1  ]]
then
	file=$1
	# read the file and do the symlinks
	while read symlink loc
	do
		if [[ "$force" = true ]]
		then
			if [[ "$testing" = true ]]
			then
				echo "$($readlink -f $symlink) $($readlink -f $loc)"
			else
				ln -fs $($readlink -f $symlink) $($readlink -f $loc)
			fi
		else
			if [[ "$testing" = true ]]
			then
				echo "$($readlink -f $symlink) $($readlink -f $loc)"
			else
				ln -s $($readlink -f $symlink) $($readlink -f $loc)
			fi
		fi
	done < $file 
elif [[ "$#" -eq 2 ]]
then
	if [[ "$force" = true ]]
	then
		if [[ "$testing" = true ]]
		then
			echo "$($readlink -f $symlink) $($readlink -f $loc)"
		else
			ln -fs $($readlink -f $symlink) $($readlink -f $loc)
		fi
	else
		if [[ "$testing" = true ]]
		then
			echo "$($readlink -f $symlink) $($readlink -f $loc)"
		else
			ln -s $($readlink -f $symlink) $($readlink -f $loc)
		fi
	fi
fi
