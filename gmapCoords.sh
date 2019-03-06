#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)

if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <GMAP input FASTA file> <GMAP output GFF file> <GMAP index GFF file>" 1>&2 
	echo "DESCRIPTION: Checks the GMAP alignments for overlaps." 1>&2
	exit 1
fi
fasta=$1
gmap_gff=$2
gff=$3
if [[ "$(awk '/\tgene\t/' $gff | wc -l)" -eq 0 ]]
then
	processGFF.sh $gff
fi
echo "Writing significant alignments to transcript-alignments.tsv..." 1>&2
echo -e "ORF\tCDS\tTranscript\tAlignment\tCoverage\tIdentity" > transcript-alignments.tsv
while read transcript
do
	transcript_begin=$(echo $transcript | awk '{print $4}')
	transcript_end=$(echo $transcript | awk '{print $5}')
	transcript_CDS=$(echo $transcript | awk '{print $9}' | awk -F ";" '{print $2}' | sed 's/Name=//')
	transcript_ORF=$(grep $transcript_CDS $fasta | awk '{print $2}')
	transcript_name=$(echo $transcript | awk '{print $9}' | awk -F ";" '{print $2}' | awk -F ":" '{print $1}' | sed 's/Name=lcl|//')
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
		if [[ "$transcript_begin" -ge "$gene_begin" && "$transcript_begin" -le "$gene_end" ]]
		then
			overlap=true
		# If transcript ends within the gene
		elif [[ "$transcript_end" -ge "$gene_begin" && "$transcript_end" -le "$gene_end" ]]
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
			if [[ "$gene_scaffold" == "$transcript_scaffold" ]]
			then
				echo -e "$transcript_ORF\t$transcript_CDS\t$transcript_name\t$gene_name\t$transcript_cov\t$transcript_pid" >> transcript-alignments.tsv
			else
				echo -e "$transcript_ORF\t$transcript_CDS\t$transcript_name\t$gene_scaffold:$transcript_begin:$transcript_end\t$transcript_cov\t$transcript_pid" >> transcript-alignments.tsv
			fi
		fi	
	done < <(grep -v '^#' $gff | awk '/\tgene\t/')
	if [[ "$mapped" = false ]]
	then
		echo -e "$transcript_ORF\t$transcript_CDS\t$transcript_name\tunmapped\t-\t-" >> transcript-alignments.tsv
	fi

done < <(grep -v '^#' $gmap_gff | awk '/\tmRNA\t/' | grep 'path1')

# List all genes that were mapped to
mapped_transcripts=$(awk '!/unmapped/{print $3}' <(tail -n +2 transcript-alignments.tsv) | sort -u)
mapped_genes=$(awk '!/unmapped/{print $4}' <(tail -n +2 transcript-alignments.tsv) | sort -u)

if [[ -e "transcript2gene.tsv" ]]
then
	rm transcript2gene.tsv
fi
echo "Writing transcripts and their corresponding genes to transcript2gene.tsv..." 1>&2
for i in $mapped_transcripts
do
	echo -en "$i\t" >> transcript2gene.tsv
	awk -v var=$i '{if($3==var) print $4}' <(tail -n +2 transcript-alignments.tsv) | sort -u | tr '\n' ' ' >> transcript2gene.tsv
	echo >> transcript2gene.tsv
done

if [[ -e "gene2transcript.tsv" ]]
then
	rm gene2transcript.tsv
fi
echo "Writing genes and their corresponding transcripts to gene2transcript.tsv..." 1>&2
for i in $mapped_genes
do
	echo -en "$i\t" >> gene2transcript.tsv
	awk -v var=$i '{if ($4==var) print $3}' <(tail -n +2 transcript-alignments.tsv) | sort -u | tr '\n' ' ' >> gene2transcript.tsv
	echo >> gene2transcript.tsv
done

if [[ "$(grep -c 'unmapped' transcript-alignments.tsv)" -ne 0 ]]
then
	echo "Writing unmapped transcripts to transcript-alignments.unmapped.tsv..." 1>&2
	head -n1 transcript-alignments.tsv > transcript-alignments.unmapped.tsv
	awk 'BEGIN{OFS="\t"}/unmapped/' transcript-alignments.tsv >> transcript-alignments.unmapped.tsv
fi

# what to do with unmapped transcripts?

echo "...Done." 1>&2
