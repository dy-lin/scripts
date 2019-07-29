#!/bin/bash

# returns reverse compliment coordinates of the GFF

if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $(basename $0) <GFF file>" 1>&2
	exit 1
fi

gff=$1
outfile=${gff%.*}.rc.gff
# write header in
head -n2 $gff > $outfile
length=$(awk 'NR==2 {print $4}' $gff)
length=$((length+1))
while IFS=$'\n'  read line
do
	echo "$line" | awk -F "\t" -v l="$length" 'BEGIN{OFS="\t"}{if($7=="+") {s="-"} else {s="+"};print $1, $2, $3, l-$5, l-$4, $6, s, $8, $9}' >> $outfile

done < <(tail -n +3 $gff)


