#!/bin/bash
PROGRAM=$(basename $0)
# use the gmap binary specified in the config file
# source $(echo $(dirname "$0") | sed 's/bin//')/gmap_config.txt
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <desired DB directory> <DB name> <FASTA file> <logfile>" 1>&2
	exit 1
fi

db=$1; shift
dbname=$1; shift
infile=$1; shift
logfile=$1

gmap_build -D $db -d $dbname --sort=none $infile > $logfile 2>&1

### EOF ###
