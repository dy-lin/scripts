#!/bin/bash

# use the gmap binary specified in the config file
#source $(echo $(dirname "$0") | sed 's/bin//')/gmap_config.txt
PROGRAM=$(basename $0)
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <DB directory path> <DB name> <threads> <format> <output prefix> <FASTA file> <logfile> <splice>" 1>&2
	echo -e "Formats:\n\tA\talignment (blast-like format)\n\t2\tGFF format" 1>&2
	exit 1
fi

db=$1; shift
dbname=$1; shift
threads=$1; shift
fmt=$1; shift
prefix=$1; shift
infile=$1; shift
logfile=$1; shift
splice=$1 # set to N by main script if doing cDNA:transcript alignments; blank otherwise

if [ -z "${splice}" ]; then
    /projects/btl/lcoombe/bin/gmap-2017-11-15_hpce/gmap-2017-11-15/bin/gmapl -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        -t $threads \
        -f $fmt \
        -x 20 \
        -O \
        --split-output=$prefix \
        -n 10000 \
        $infile > $logfile 2>&1
else
    /projects/btl/lcoombe/bin/gmap-2017-11-15_hpce/gmap-2017-11-15/bin/gmapl -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        --nosplicing \
        -t $threads \
        -f $fmt \
        -x 20 \
        -O \
        --split-output=$prefix \
        -n 10000 \
        $infile > $logfile 2>&1
fi

### EOF ###
