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
	echo "USAGE: $PROGRAM <GMAP input FASTA file> <GMAP output GFF file> <GMAP index GFF file>" 1>&2 
	echo "DESCRIPTION: Checks the GMAP alignments for overlaps." 1>&2
	echo -e "OPTIONS:\n\t-a\tShow alignments to other scaffolds\n\t-h\tShow help menu\n\t-t\tMap to transcripts instead of scaffold" 1>&2
	exit 1
fi
fasta=$1
gmap_gff=$2
gff=$3
dir=$(dirname $gmap_gff)
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
echo "Writing significant alignments to transcript-alignments.tsv..." 1>&2
echo -e "ID\tORF\tCDS\tTranscript\tAlignment\tCoverage\tIdentity" > $dir/transcript-alignments.tsv
while read transcript
do
	transcript_id=$(echo $transcript | awk '{print $9}' | awk -F ";" '{print $1}')
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
				echo -e "$transcript_id\t$transcript_ORF\t$transcript_CDS\t$transcript_name\t$gene_name\t$transcript_cov\t$transcript_pid" >> $dir/transcript-alignments.tsv
			else
				if [[ "$aln" = true ]]
				then
					echo -e "$transcript_id\t$transcript_ORF\t$transcript_CDS\t$transcript_name\t$gene_scaffold:$transcript_begin:$transcript_end\t$transcript_cov\t$transcript_pid" >> $dir/transcript-alignments.tsv
				fi
			fi
		fi	
	done < <(grep -v '^#' $gff | awk -v var=$feature '$0 ~ var' )
	if [[ "$mapped" = false ]]
	then
		echo -e "$transcript_id\t$transcript_ORF\t$transcript_CDS\t$transcript_name\tunmapped\t-\t-" >> $dir/transcript-alignments.tsv
	fi

# done < <(grep -v '^#' $gmap_gff | awk '/\tmRNA\t/' | grep 'identity=100' | grep 'coverage=100')
# done < <(grep -v '^#' $gmap_gff | awk '/\tmRNA\t/' | grep 'path1')
done < <(grep -v '^#' $gmap_gff | awk '/\tmRNA\t/')

# List all genes that were mapped to
mapped_transcripts=$(awk '!/unmapped/{if($6>=95 && $7>=95) print $4}' <(tail -n +2 $dir/transcript-alignments.tsv) | sort -u)
# unmapped_transcripts=$(awk '/unmapped/ {print $4}' <(tail -n +2 $dir/transcript-alignments.tsv) | sort -u)
mapped_genes=$(awk '!/unmapped/{if($6>=95 && $7>=95) print $5}' <(tail -n +2 $dir/transcript-alignments.tsv) | sort -u)
# unmapped_genes=$(awk '/unmapped/ {print $5}' <(tail -n +2 $dir/transcript-alignments.tsv) | sort -u )
leftover=$(grep -vf <(awk '{print $1}' <(echo "$mapped_genes")) <(awk -F "ID=" '/\tgene\t/ {print $2}' $gff | awk -F ";" '{print $1}'| sort -u))
# echo $leftover
# echo $mapped_genes

if [[ -e "$dir/transcript2gene.tsv" ]]
then
	rm $dir/transcript2gene.tsv
fi

echo "Writing transcripts and their corresponding genes to transcript2gene.tsv..." 1>&2
for i in $mapped_transcripts
do
	echo -en "$i\t" >> $dir/transcript2gene.tsv
	result=$(awk -v var=$i '/.mrna1/{if($4==var) print $5}' <(tail -n +2 $dir/transcript-alignments.tsv) | sort -u | tr '\n' ' ')
	if [[ -z "$result" ]]
	then
		echo "-" >> $dir/transcript2gene.tsv
	else
		echo "$result" >> $dir/transcript2gene.tsv
	fi
done

if [[ -e "$dir/gene2transcript.tsv" ]]
then
	rm $dir/gene2transcript.tsv
fi

echo "Writing genes and their corresponding transcripts to gene2transcript.tsv..." 1>&2
for i in $mapped_genes
do
	echo -en "$i\t" >> $dir/gene2transcript.tsv
	result=$(awk -v var=$i '/.mrna1/ {if ($5==var) print $4}' <(tail -n +2 $dir/transcript-alignments.tsv) | sort -u | tr '\n' ' ')
	if [[ -z "$result" ]]
	then
		echo "-" >> $dir/gene2transcript.tsv
	else
		echo "$result" >> $dir/gene2transcript.tsv
	fi
done

echo -e "Gene\t$(basename $dir)" > $dir/unique_transcripts.tsv
echo "Writing genes and their corresponding CDS to gene2cds.tsv..." 1>&2
echo "Writing genes and their unique number of transcripts to unique_transcripts.tsv..." 1>&2
if [[ -e "$dir/gene2cds.tsv" ]]
then
	rm $dir/gene2cds.tsv
fi
for i in $mapped_genes
do
	echo -en "$i\t" >> $dir/gene2cds.tsv
	awk -v var=$i '/.mrna1/ {if ($5==var) print $3}' <(tail -n +2 $dir/transcript-alignments.tsv) | sort -u | tr '\n' ' ' >> $dir/gene2cds.tsv
	# GET CDS SEQ FOR EACH MATCH
	num=$(uniqtag -k 252 <(seqtk subseq $fasta <(awk -v var=$i '/.mrna1/ {if ($5==var) print $3}' <(tail -n +2 $dir/transcript-alignments.tsv) | sort -u)) | grep -c '\-1$')
	echo -e "$i\t$num" >> $dir/unique_transcripts.tsv
	echo >> $dir/gene2cds.tsv
done
# echo "hi" 1>&2
# Write missed genes to each file
for gene in $leftover
do
#	echo $gene
	echo -e "$gene\t-" >> $dir/gene2transcript.tsv
	echo -e "$gene\t0" >> $dir/unique_transcripts.tsv
done
if [[ "$(grep -c 'unmapped' $dir/transcript-alignments.tsv)" -ne 0 ]]
then
	echo "Writing unmapped transcripts to transcript-alignments.unmapped.tsv..." 1>&2
	head -n1 $dir/transcript-alignments.tsv > $dir/transcript-alignments.unmapped.tsv
	awk 'BEGIN{OFS="\t"}/unmapped/' $dir/transcript-alignments.tsv >> $dir/transcript-alignments.unmapped.tsv
fi

# SORT ALL THE writen files
sort -k1 -o $dir/gene2transcript.tsv $dir/gene2transcript.tsv
sort -k1 -o $dir/transcript2gene.tsv $dir/transcript2gene.tsv
sort -k1 -o $dir/gene2cds.tsv $dir/gene2cds.tsv
cat <(head -n1 $dir/unique_transcripts.tsv) <(tail -n +2 $dir/unique_transcripts.tsv | sort -k1) > temp
mv temp $dir/unique_transcripts.tsv
echo "...Done." 1>&2
