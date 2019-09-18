#!/bin/bash
set -euo pipefail
PROGRAM=$(basename $0)

if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <Files/symlinks to copy>" 1>&2
	echo "DESCRIPTION: Takes files/symlinks and copies to them to the current directory."
	exit 1
fi

for file in $@
do
	if [ -L "$file" ]
	then
		cp "$(readlink -f $file)" $(basename $file)
	else
		cp $file .
	fi
done

# Symlink those original files to this new directory?
