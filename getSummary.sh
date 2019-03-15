#!/bin/bash
#set -eu -o pipefail

if [[ "$#" -ne 2 ]]
then
	echo "USAGE: $(basename $0) <bit score threshold> <gene to transcript tsv>" 1>&2
	exit 1
fi

threshold=$1
dict=$2

# Get all AMPs of a class
# amp=$(basename $(pwd))
amp="defensins"
class=$(echo $amp | sed 's/s$//')
outfile=${amp}_summary.jira
tissues=$(ls -l | awk '/^d/ {print $9}' | sort -u)

max=$(($(echo $tissues | wc -w) +1 ))
echo -n "||$amp||" > $outfile
count=1

# header
# On the last tissue, write with newline.
for tissue in $tissues
do
	count=$((count+1))
	if [[ "$count" -lt "$max" ]]
	then
		echo -n "$tissue||" >> $outfile
	else
		echo "$tissue||" >> $outfile
	fi
done

# iterate through genes and tally numbers

for gene in $(awk -F "\t" '!/unmapped/ {print $1}' $dict | sort -u)
do
	echo -ne "|$gene|" >> $outfile
	count=1
	for tissue in $tissues
	do
		count=$((count+1))
		num=$(grep -cf <(for file in ${tissue}/bs${threshold}/transcripts/${amp}-only/*.fa; do basename $file ".transcript.fa" ; done | sort -u) <(awk -v var=$gene -F "\t" '{if($1==var) print $2}' $dict | tr ' ' '\n' | sort -u ))
		if [[ "$count" -lt "$max" ]]
		then
			echo -n "$num|" >> $outfile
		else
			echo "$num|" >> $outfile
		fi
	done
done

		
# check if transcripts exist for that tissue, if it does (/) if not, (x)

# for prot in $(for file in $(ls */bs${threshold}/transcripts/${amp}-only/*); do basename $file; done | sort -u)
# do
# 	echo -n "|$(basename $prot ".transcript.fa")|" >> $outfile
# 	count=1
# 	for tissue in $tissues
# 	do
# 		count=$((count+1))
# 		if [[ -s "${tissue}/bs${threshold}/transcripts/${amp}-only/${prot}" ]]
# 		then
# 			if [[ "$count" -lt "$max" ]]
# 			then
# 				echo -n "(/)|" >> $outfile
# 			else
# 				echo "(/)|" >> $outfile
# 			fi
# 		else
# 			if [[ "$count" -lt "$max" ]]
# 			then
# 				echo -n "(x)|" >> $outfile
# 			else
# 				echo "(x)|" >> $outfile
# 			fi
# 		fi
# 	done
# done
# echo -n "|Total|" >> $outfile
# count=1
# for tissue in $tissues
# do
# 	num_transcripts=$(ls ${tissue}/bs${threshold}/transcripts/${amp}-only/* | wc -l)
# 	if [[ "$count" -lt "$max" ]]
# 	then
# 		echo -n "$num_transcripts|" >> $outfile
# 	else
# 		echo "$num_transcripts|" >> $outfile
# 	fi
# done

echo >> $outfile

echo "||tissue||\# transcripts||\# ORFs||jackhmmer hits||blast hits ($amp/$class-like)|| corresponding \# transcripts || $amp (total) || $amp (complete) || $amp (partial) || ${class}-like || threshold ||" >> $outfile
for tissue in $tissues
do
#	echo $tissue >> $outfile
	assembly="../assemblies/${tissue}"
	while [[ ! -e "$assembly" ]]
	do
		assembly="../$assembly"
	done

	if [[ "$(readlink -f $assembly)" == *.gz ]]
	then
		num_transcriptome=$(zgrep -c '^>' $assembly)
	else
		num_transcriptome=$(grep -c '^>' $assembly)
	fi
	ORFfile="${tissue}/${tissue}.ORFs.fa"

	if [[ ! -e "$ORFfile" ]]
	then
		ORFfile="${ORFfile}a"
	fi
	num_ORFs=$(grep -c '^>' $ORFfile)

	hitsfile="${tissue}/bs${threshold}/jackhmmer-hits.fa"
	if [[ ! -e "$hitsfile" ]]
	then
		hitsfile="${hitsfile}a"
	fi
	num_hits=$(grep -c '^>' $hitsfile)

	blasthitsfile="${tissue}/bs${threshold}/jackhmmer-blast-hits.fa"
	if [[ ! -e "$blasthitsfile" ]]
	then
		blasthitsfile="${blasthitsfile}a"
	fi
	num_blast=$(grep -c '^>' $blasthitsfile) 
	num_transcripts=$(ls ${tissue}/bs${threshold}/transcripts/*.fa | wc -l)
	if [[ -e ${tissue}/bs${threshold}/jackhmmer-blast-hits.trimmed.faa ]]
	then
		num_amplike=$(grep -c "$class-like" ${tissue}/bs${threshold}/jackhmmer-blast-hits.trimmed.faa)
		num_amp_partial=$(grep '^>' ${tissue}/bs${threshold}/jackhmmer-blast-hits.trimmed.faa | grep -v "$class-like" | grep -c 'partial')
	else
		num_amplike=0
		num_amp_partial=$(grep '^>' $blasthitsfile | grep -v "$class-like" | grep -c 'partial')
	fi
	num_amp=$((num_blast-num_amplike))
	num_amp_complete=$((num_amp-num_amp_partial))
	echo "|$tissue|$num_transcriptome|$num_ORFs|$num_hits|$num_blast|$num_transcripts|$num_amp|$num_amp_complete|$num_amp_partial|$num_amplike|$threshold|" >> $outfile
done
