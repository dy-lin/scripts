#!/bin/bash
# set -eu -o pipefail
PROGRAM=$(basename $0)
gethelp=false
map2transcript=false
aln=false
while getopts :ahgt opt
do
	case $opt in
		a) aln=true;;
		h) gethelp=true;;
		t) map2transcript=true;;
		\?) echo "$PROGRAM: invalid option $opt" 1>&2; exit 1;;
	esac
done
shift $((OPTIND-1))


if [[ "$#" -eq 0 || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM <GMAP input FASTA file> <GMAP output GFF file> <GMAP index GFF file> [output directory] [output prefix]" 1>&2 
	echo "DESCRIPTION: Checks the GMAP alignments for overlaps." 1>&2
	echo -e "OPTIONS:\n\t-a\tShow alignments to other scaffolds\n\t-h\tShow help menu\n\t-t\tMap to transcripts instead of scaffold" 1>&2
	exit 1
fi
fasta=$1
gmap_gff=$2
gff=$3
if [[ "$#" -eq 4 ]]
then
	# If the fourth argument is a directory. 
	if [[ -d "$4" ]]
	then
		echo "Output directory detected..." 1>&2
		dir=$(echo $4 | sed 's/\/$//')
	else
		echo "Prefix detected. Using default directory..." 1>&2
		prefix=${4}.
		dir=$(dirname $gmap_gff)
	fi
elif [[ "$#" -eq 5 ]]
then
	echo "Output directory and prefix detected..." 1>&2
	prefix=${5}.
	dir=$(echo $4 | sed 's/\/$//')
else
	echo "Using default directory..." 1>&2
	dir=$(dirname $gmap_gff)
fi
if [[ "$map2transcript" = true ]]
then
	feature="\tmRNA\t"
	echo "Aligning to transcripts..." 1>&2
else
	if [[ "$(awk '/\tgene\t/' $gff | wc -l)" -eq 0 ]]
	then
		processGFF.sh $gff
	fi
	feature="\tgene\t"
	echo "Aligning to genes..." 1>&2
fi
echo "Writing significant alignments to ${prefix}transcript-alignments.tsv..." 1>&2
echo -e "ID\tORF\tCDS\tTranscript\tAlignment\tCoverage\tIdentity" > $dir/${prefix}transcript-alignments.tsv
while read transcript
do
	transcript_id=$(echo $transcript | awk '{print $9}' | awk -F ";" '{print $1}')
	transcript_begin=$(echo $transcript | awk '{print $4}')
	transcript_end=$(echo $transcript | awk '{print $5}')
	if [[ "$(grep -c 'lcl' $fasta)" -ne 0 ]]
	then
		transcript_CDS=$(echo $transcript | awk '{print $9}' | awk -F ";" '{print $2}' | sed 's/Name=//')
	else
		transcript_CDS="-"
	fi
	if [[ "$(grep -c 'ORF' $fasta)" -ne 0 ]]
	then
		transcript_ORF=$(grep $transcript_CDS $fasta | awk '{print $2}')
	else
		transcript_ORF="-"
	fi
	if [[ "$(grep -c 'lcl' $fasta)" -ne 0 ]]
	then
		transcript_name=$(echo $transcript | awk '{print $9}' | awk -F ";" '{print $2}' | awk -F ":" '{print $1}' | sed 's/Name=lcl|//')
	else
		transcript_name=$(echo $transcript | awk '{print $9}' | awk -F ";" '{print $2}' | awk -F ":" '{print $1}' | sed 's/Name=//')
	fi
	transcript_cov=$(echo $transcript | awk '{print $9}' | awk -F ";" '{print $4}' | sed 's/coverage=//')
	transcript_pid=$(echo $transcript | awk '{print $9}' | awk -F ";" '{print $5}' | sed 's/identity=//')
	transcript_scaffold=$(echo $transcript | awk '{print $1}')
	mapped=false
	while read gene
	do
		gene_begin=$(echo $gene | awk '{print $4}')
		gene_end=$(echo $gene | awk '{print $5}')
		gene_scaffold=$(echo $gene | awk '{print $1}')
		gene_name=$(echo $gene | awk '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//')

		# If transcript began within the gene
		if [[ "$transcript_begin" -ge "$gene_begin" && "$transcript_begin" -lt "$gene_end" ]]
		then
			overlap=true
		# If transcript ends within the gene
		elif [[ "$transcript_end" -gt "$gene_begin" && "$transcript_end" -le "$gene_end" ]]
		then
			overlap=true
		# If the transcript is within the gene entirely
		elif [[ "$transcript_begin" -ge "$gene_begin" && "$transcript_end" -le "$gene_end"  ]]
		then
			overlap=true
		# If the transcript covers the gene entirely
		elif [[ "$transcript_begin" -le "$gene_begin" && "$transcript_end" -ge "$gene_end" ]]
		then
			overlap=true
		# If the transcript is fully before the gene
		elif [[ "$transcript_begin" -lt "$gene_begin" && "$transcript_end" -le "$gene_begin" ]]
		then
			overlap=false
		# If the transcript is fully after the gene
		elif [[ "$transcript_begin" -ge "$gene_end" && "$transcript_end" -gt "$gene_end" ]]
		then
			overlap=false
		fi

		if [[ "$overlap" = true ]]
		then
			mapped=true
			if [[ "$transcript_scaffold" =~ $gene_scaffold ]]
			then
				echo -e "$transcript_id\t$transcript_ORF\t$transcript_CDS\t$transcript_name\t$gene_name\t$transcript_cov\t$transcript_pid" >> $dir/${prefix}transcript-alignments.tsv
			else
				if [[ "$aln" = true ]]
				then
					echo -e "$transcript_id\t$transcript_ORF\t$transcript_CDS\t$transcript_name\t$gene_scaffold:$transcript_begin:$transcript_end\t$transcript_cov\t$transcript_pid" >> $dir/${prefix}transcript-alignments.tsv
				fi
			fi
		fi	
	done < <(grep -v '^#' $gff | awk -v var=$feature '$0 ~ var' )
	if [[ "$mapped" = false ]]
	then
		echo -e "$transcript_id\t$transcript_ORF\t$transcript_CDS\t$transcript_name\tunmapped\t-\t-" >> $dir/${prefix}.transcript-alignments.tsv
	fi

# done < <(grep -v '^#' $gmap_gff | awk '/\tmRNA\t/' | grep 'identity=100' | grep 'coverage=100')
# done < <(grep -v '^#' $gmap_gff | awk '/\tmRNA\t/' | grep 'path1')
done < <(grep -v '^#' $gmap_gff | awk '/\tmRNA\t/')

# List all genes that were mapped to
mapped_transcripts=$(awk '!/unmapped/{if($6>=95 && $7>=95) print $4}' <(tail -n +2 $dir/${prefix}transcript-alignments.tsv) | sort -u)
# unmapped_transcripts=$(awk '/unmapped/ {print $4}' <(tail -n +2 $dir/transcript-alignments.tsv) | sort -u)
mapped_genes=$(awk '!/unmapped/{if($6>=95 && $7>=95) print $5}' <(tail -n +2 $dir/${prefix}transcript-alignments.tsv) | sort -u)
# unmapped_genes=$(awk '/unmapped/ {print $5}' <(tail -n +2 $dir/transcript-alignments.tsv) | sort -u )
leftover=$(grep -vf <(awk '{print $1}' <(echo "$mapped_genes")) <(awk -F "ID=" '/\tgene\t/ {print $2}' $gff | awk -F ";" '{print $1}'| sort -u))
#echo $leftover
# echo $mapped_genes

if [[ -e "$dir/${prefix}transcript2gene.tsv" ]]
then
	rm $dir/${prefix}transcript2gene.tsv
fi

echo "Writing transcripts and their corresponding genes to ${prefix}transcript2gene.tsv..." 1>&2
for i in $mapped_transcripts
do
	echo -en "$i\t" >> $dir/${prefix}transcript2gene.tsv
	result=$(awk -v var=$i '/.mrna1/{if($4==var) print $5}' <(tail -n +2 $dir/${prefix}transcript-alignments.tsv) | sort -u | tr '\n' ' ')
	if [[ -z "$result" ]]
	then
		echo "-" >> $dir/${prefix}transcript2gene.tsv
	else
		echo "$result" >> $dir/${prefix}transcript2gene.tsv
	fi
done
if [[ -e "$dir/${prefix}gene2transcript.tsv" ]]
then
	rm $dir/${prefix}gene2transcript.tsv
fi

echo "Writing genes and their corresponding transcripts to ${prefix}gene2transcript.tsv..." 1>&2
for i in $mapped_genes
do
	echo -en "$i\t" >> $dir/${prefix}gene2transcript.tsv
	result=$(awk -v var=$i '/.mrna1/ {if ($5==var) print $4}' <(tail -n +2 $dir/${prefix}transcript-alignments.tsv) | sort -u | tr '\n' ' ')
	if [[ -z "$result" ]]
	then
		echo "-" >> $dir/${prefix}gene2transcript.tsv
	else
		echo "$result" >> $dir/${prefix}gene2transcript.tsv
	fi
done
if [[ "$(basename $dir)" == "." ]]
then
	echo -e "Gene\t$(basename $(pwd))" > $dir/${prefix}unique_transcripts.tsv
else
	echo -e "Gene\t$(basename $dir)" > $dir/${prefix}unique_transcripts.tsv
fi
echo "Writing genes and their corresponding CDS to ${prefix}gene2cds.tsv..." 1>&2
echo "Writing genes and their unique number of transcripts to ${prefix}unique_transcripts.tsv..." 1>&2
if [[ -e "$dir/${prefix}gene2cds.tsv" ]]
then
	rm $dir/${prefix}gene2cds.tsv
fi
for i in $mapped_genes
do
	echo -en "$i\t" >> $dir/${prefix}gene2cds.tsv
	awk -v var=$i '/.mrna1/ {if ($5==var) print $3}' <(tail -n +2 $dir/${prefix}transcript-alignments.tsv) | sort -u | tr '\n' ' ' >> $dir/${prefix}gene2cds.tsv
	# GET CDS SEQ FOR EACH MATCH
	num=$(uniqtag -k 252 <(seqtk subseq $fasta <(awk -v var=$i '/.mrna1/ {if ($5==var) print $3}' <(tail -n +2 $dir/${prefix}transcript-alignments.tsv) | sort -u)) | grep -c '\-1$')
	echo -e "$i\t$num" >> $dir/${prefix}unique_transcripts.tsv
	echo >> $dir/${prefix}gene2cds.tsv
done
# echo "hi" 1>&2
# Write missed genes to each file
for gene in $leftover
do
#	echo $gene
	echo -e "$gene\t-" >> $dir/${prefix}gene2transcript.tsv
	echo -e "$gene\t0" >> $dir/${prefix}unique_transcripts.tsv
done
if [[ "$(grep -c 'unmapped' $dir/${prefix}transcript-alignments.tsv)" -ne 0 ]]
then
	echo "Writing unmapped transcripts to ${prefix}transcript-alignments.unmapped.tsv..." 1>&2
	head -n1 $dir/${prefix}transcript-alignments.tsv > $dir/${prefix}transcript-alignments.unmapped.tsv
	awk 'BEGIN{OFS="\t"}/unmapped/' $dir/${prefix}transcript-alignments.tsv >> $dir/${prefix}transcript-alignments.unmapped.tsv
fi
# SORT ALL THE writen files
sort -k1 -o $dir/${prefix}gene2transcript.tsv $dir/${prefix}gene2transcript.tsv
sort -k1 -o $dir/${prefix}transcript2gene.tsv $dir/${prefix}transcript2gene.tsv
sort -k1 -o $dir/${prefix}gene2cds.tsv $dir/${prefix}gene2cds.tsv
cat <(head -n1 $dir/${prefix}unique_transcripts.tsv) <(tail -n +2 $dir/${prefix}unique_transcripts.tsv | sort -k1) > temp
mv temp $dir/${prefix}unique_transcripts.tsv
echo "...Done." 1>&2
