#!/bin/bash

# returns reverse compliment coordinates of the GFF

if [[ "$#" -ne 2 ]]
then
	echo "USAGE: $(basename $0) [OPTIONS] <file>" 1>&2
	echo "DESCRIPTION: Changes the coordinates of a file to correspond with the reverse complement." 1>&2
	echo -e "OPTIONS:\n\t--gff\tIndicates file is a GFF file\n\t--tsv\tIndicates file is a TSV file" 1>&2
	exit 1
fi
file=$1
shift
case $file in
	--gff)
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
		;;
	--tsv)
		tsv=$1
		# where tsv is formatted scaffoldname[dot]tsv
		prefix=${tsv%%.*}
		outfile=${prefix}.rc.tsv
		if [[ -e "$outfile" ]]
		then
			rm "$outfile"
		fi
		if [[ -e ${prefix}.scaffold.fa ]]
		then
			length=$(seqtk comp ${prefix}.scaffold.fa | awk '{print $2}')
			length=$((length+1))
		else
			read -p "Path to scaffold: " scaffold
			length=$(seqtk comp $scaffold | awk '{print $2}')
			length=$((length+1))
		fi
			while IFS=$'\n' read line
			do
				echo "$line" | awk -F "\t" -v l="$length" 'BEGIN{OFS="\t"}{print l-$2, l-$1, $3}' >> $outfile
			done < $tsv
		;;
esac



