#!/bin/bash
set -euo pipefail
PROGRAM=$(basename $0)

if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <File(s) to link>" 1>&2
	echo "DESCRIPTION: Takes file(s) to link and links them to the current directory." 1>&2
	exit 1
fi

for file in $@
do
	file_path=$(readlink -f $file)
	filename=$(basename $file)
	ln -s $file_path $filename
done


