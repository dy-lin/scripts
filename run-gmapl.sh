#!/bin/bash

# use the gmap binary specified in the config file
#source $(echo $(dirname "$0") | sed 's/bin//')/gmap_config.txt
PROGRAM=$(basename $0)
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <GMAP-index path> <threads> <output prefix> <FASTA file> [splice]" 1>&2
	exit 1
fi
db=$(dirname $1)
dbname=$(basename $1); shift
threads=$1; shift
prefix=$1; shift
infile=$1; shift
splice=$1 # set to N by main script if doing cDNA:transcript alignments; blank otherwise

if [ -z "${splice}" ]; then
    /projects/btl/lcoombe/bin/gmap-2017-11-15_hpce/gmap-2017-11-15/bin/gmapl -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        -t $threads \
        -f 2 \
        -x 20 \
        -O \
        -n 10000 \
        $infile > ${prefix}.gmapl.gff 2> ${prefix}.gmapl.log
else
    /projects/btl/lcoombe/bin/gmap-2017-11-15_hpce/gmap-2017-11-15/bin/gmapl -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        --nosplicing \
        -t $threads \
        -f 2 \
        -x 20 \
        -O \
        -n 10000 \
        $infile > ${prefix}.gmapl.gff 2> ${prefix}.gmapl.log
fi

if [ -z "${splice}" ]; then
    /projects/btl/lcoombe/bin/gmap-2017-11-15_hpce/gmap-2017-11-15/bin/gmapl -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        -t $threads \
        -A \
        -x 20 \
        -O \
        -n 10000 \
        $infile > ${prefix}.gmapl.aln 2>> ${prefix}.gmapl.log
else
    /projects/btl/lcoombe/bin/gmap-2017-11-15_hpce/gmap-2017-11-15/bin/gmapl -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        --nosplicing \
        -t $threads \
        -A \
        -x 20 \
        -O \
        -n 10000 \
        $infile > ${prefix}.gmapl.aln 2>> ${prefix}.gmapl.log
fi

### EOF ###
