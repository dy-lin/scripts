#!/bin/bash
set -euo pipefail
PROGRAM=$(basename $0)
if [[ "$#" -ne 5 ]]
then
	echo "USAGE: $PROGRAM <class of AMPs> <literature AMPs (prot) > <literature AMPs (nucl)> <genome assembly> <output prefix>"
	exit 1
fi
# if [[ "$(pwd)" != /annotation/amp ]]
# then
#	echo "Wrong working directory." 1>&2
#	echo "Working directory must be: /projects/spruceup/scratch/spruce/genotype/annotation/amp" 1>&2
#	exit 1
# fi
class=$1
if [[ "$class" != *s ]]
then
	class="${class}s"
fi
amps=$2
mRNA=$3
scaffold=$4
prefix=$5

# Make tblastn-genome directory
mkdir -p tblastn-genome/${class}
cd tblastn-genome/${class}
blast.sh tblastn $amps $scaffold
blast.sh blastn $mRNA $scaffold
cd ../..


# Make gmap-genome directory
mkdir -p gmap-genome/${class}
cd gmap-genome/${class}
run-gmapl.sh ../GMAP-index 48 $prefix $mRNA
cd ../..


