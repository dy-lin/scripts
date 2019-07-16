#!/bin/bash

# use the gmap binary specified in the config file
#source $(echo $(dirname "$0") | sed 's/bin//')/gmap_config.txt

PROGRAM=$(basename $0)
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <GMAP-index path> <threads> <output prefix> <FASTA file> [splice]" 1>&2
	echo "Note: To specify an output directory, use <output prefix>." 1>&2
	exit 1
fi
db=$(dirname $1)
dbname=$(basename $1);shift
threads=$1; shift
if [[ -d "$1" ]]
then
	# If a entire directory path is given as a prefix, use a default prefix.
	prefix=output
	dir=$1; shift
else
	prefix=$(basename $1)
	dir=$(dirname $1); shift
fi
infile=$1; shift
splice=$1 # set to N by main script if doing cDNA:transcript alignments; blank otherwise

if [ -z "${splice}" ]; then
    /projects/btl/lcoombe/bin/gmap-2017-11-15_hpce/gmap-2017-11-15/bin/gmap -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        -t $threads \
        -A \
        -x 20 \
        -O \
        -n 10000 \
        $infile > ${dir}/${prefix}.gmap.aln 2> ${dir}/${prefix}.gmap.log
	if [[ "$?" -eq 9 ]]
	then
		echo "Genome is too large to run gmap. Running gmapl..." 1>&2
		rm *.gmap.*
		run-gmapl.sh $db/$dbname $threads ${dir}/$prefix $infile
	fi
else
    /projects/btl/lcoombe/bin/gmap-2017-11-15_hpce/gmap-2017-11-15/bin/gmap -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        --nosplicing \
        -t $threads \
        -A \
        -x 20 \
        -O \
        -n 10000 \
        $infile > ${dir}/${prefix}.gmap.aln 2> ${dir}/${prefix}.gmap.log
	if [[ "$?" -eq 9 ]]
	then
		echo "Genome is too large to run gmap. Running gmapl..." 1>&2
		rm *.gmap.*
		run-gmapl.sh $db/$dbname $threads $dir/$prefix $infile $splice
	fi
fi

if [ -z "${splice}" ]; then
    /projects/btl/lcoombe/bin/gmap-2017-11-15_hpce/gmap-2017-11-15/bin/gmap -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        -t $threads \
        -f 2 \
        -x 20 \
        -O \
        -n 10000 \
        $infile > $dir/${prefix}.gmap.gff 2>> $dir/${prefix}.gmap.log
	if [[ "$?" -eq 9 ]]
	then
		echo "Genome is too large to run gmap. Running gmapl..." 1>&2
		rm *.gmap.*
		run-gmapl.sh $db/$dbname $threads $dir/$prefix $infile
	fi
else
    /projects/btl/lcoombe/bin/gmap-2017-11-15_hpce/gmap-2017-11-15/bin/gmap -D $db -d $dbname \
        --max-intronlength-ends=1000000 \
        --max-intronlength-middle=1000000 \
        --totallength=20000000 \
        --nosplicing \
        -t $threads \
        -f 2 \
        -x 20 \
        -O \
        -n 10000 \
        $infile > $dir/${prefix}.gmap.gff 2>> $dir/${prefix}.gmap.log
	if [[ "$?" -eq 9 ]]
	then
		echo "Genome is too large to run gmap. Running gmapl..." 1>&2
		rm *.gmap.*
		run-gmapl.sh $db/$dbname $threads $dir/$prefix $infile $splice
	fi
fi

### EOF ###
