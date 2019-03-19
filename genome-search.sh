#!/bin/bash

if [[ "$(pwd)" != */annotation/amp ]]
then
	echo "Wrong working directory." 1>&2
	echo "Working directory must be: /projects/spruceup/scratch/spruce/genotype/annotation/amp" 1>&2
	exit 1
fi

amps=$1
mRNA=$2
scaffold=$3
prefix=$4

# Make tblastn-genome directory
mkdir -p tblastn-genome
cd tblastn-genome
blast.sh nucl $amps $scaffold
cd ..
# Make gmap-genome directory
mkdir -p gmap-genome
run-gmap.sh GMAP-index 48 $prefix $mRNA

# Run gmap-coords and then jira summary
