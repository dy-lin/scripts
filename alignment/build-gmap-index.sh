#!/bin/bash
PROGRAM=$(basename $0)
# use the gmap binary specified in the config file
# source $(echo $(dirname "$0") | sed 's/bin//')/gmap_config.txt
if [[ "$#" -ne 2 ]]
then
	echo "USAGE: $PROGRAM <GMAP-index path> <FASTA file>" 1>&2
	exit 1
fi

db=$(dirname $1)
dbname=$(basename $1)
infile=$2

gmap_build -D $db -d $dbname --sort=none $infile > gmap-build.log 2>&1

### EOF ###
