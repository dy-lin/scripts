#!/bin/bash
set -euo pipefail
if [[ "$#" -ne 3 ]]
then
	echo "USAGE: $(basename $0) <DE genes> <significant GOterms> <GFF file>" 1>&2
	exit 1
fi

DEgenes=$1
sigGO=$2
gff=$3

dir=$(dirname $sigGO)


gene=$(grep -f $sigGO <(grep -f $DEgenes $gff) | awk -F "\t" '/\tgene\t/ {print $9}' | awk -F "ID=" '{print $2}' | awk -F ";" '{print $1}')
name=$(grep -f $sigGO <(grep -f $DEgenes $gff) | awk -F "\t" '/\tgene\t/ {print $9}' | awk -F "Note=" '{print $2}' | awk -F ";" '{print $1}')
echo -e "Gene\tFunction" > $dir/GOterms.sig.genes.tsv
paste <(echo "$gene") <(echo "$name") >> $dir/GOterms.sig.genes.tsv
