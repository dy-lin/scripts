#!/bin/bash
#dir="/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/jackhmmer-transcriptome/defensins/"
if [[ "$#" -ne 3 && "$#" -ne 4 ]]
then
	echo "USAGE: $(basename $0) <transcriptome transcripts> <gene annotation transcripts> <genotype> [GMAP-index path]" 1>&2
	exit 1
fi
query=$1
subject=$2
genotype=$3
if [[ -z "$4" ]]
then
	index=GMAP-index
	build-gmap-index.sh $index $subject
else
	index=$(readlink -f $4)
fi
if [[ "$genotype" == "PG29" ]]
then
	tissues="bark embryo flush_bud mature_needle megagametophyte seed_germination xylem young_buds"
else
	tissues="control wound gallery Cort_Par Dev_SC"
fi


run-gmap.sh $index 48 $genotype $query
blat $subject $query output.psl
blast.sh blastn $query $subject .

