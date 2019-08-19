#!/bin/bash
PROGRAM=$(basename $0)
set -euo pipefail
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <reference GFF file> <tabular BLAST output(s)>" 1>&2
	echo "DESCRIPTION: Takes a reference GFF file and tabular BLAST output(s) to check if blast alignments are annotated." 1>&2
	exit 1
fi
gff=$1
if [[ "$(grep -c '	gene	' $gff)" -eq 0 ]]
then
	processGFF.sh -i $gff
fi

shift

for blast in $@
do
	if [[ ! -e "$blast" ]]
	then
		echo "$(basename $blast) does not exist." 1>&2
		exit 1
	fi
	filename=${blast%.*}
	if [[ -e "${filename}.novelhits.gff" ]]
	then
		rm ${filename}.novelhits.gff
	fi
		touch ${filename}.novelhits.gff
	# Read the blast alignment
	while IFS=$'\t' read qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs
	do
		# get the corresponding record in the GFF
		if [[ "$sstart" -lt "$send" ]]
		then
			strand="+"
		else
			strand="-"
		fi
		if [[ "$(grep -c "$saccver" ${filename}.novelhits.gff)" -eq 0 ]]
		then
			count=1
		fi
		
		if [[ "$(grep -c "$saccver" $gff)" -eq 0 ]]
		then
			echo -e "$saccver\tmanual\tgene\t$sstart\t$send\t.\t$strand\t.\tID=manual-${saccver}-blast-gene-0.${count};idenity=$pident;mismatches=$mismatch;coverage=$qcovs" >> ${filename}.novelhits.gff
			echo -e "$saccver\tmanual\tmRNA\t$sstart\t$send\t.\t$strand\t.\tID=manual-${saccver}-blast-gene-0.${count}-mRNA-1;Parent=manual-${saccver}-blast-gene-0.${count};idenity=$pident;mismatches=$mismatch;coverage=$qcovs" >> ${filename}.novelhits.gff
			count=$((count+1))
		else
			IFS=$'\n'
			for line in $(awk -v var="$saccver" '{if($1==var && $3=="gene") print}' $gff)
			do
				begin=$(echo $line | awk '{print $4}')
				end=$(echo $line | awk '{print $5}')
				# If transcript began within the gene
				if [[ "$sstart" -ge "$begin" && "$sstart" -lt "$end" ]]
				then
					overlap=true
				# If transcript ends within the gene
				elif [[ "$send" -gt "$begin" && "$send" -le "$end" ]]
				then
					overlap=true
				# If the transcript is within the gene entirely
				elif [[ "$sstart" -ge "$begin" && "$send" -le "$end"  ]]
				then
					overlap=true
				# If the transcript covers the gene entirely
				elif [[ "$sstart" -le "$begin" && "$send" -ge "$end" ]]
				then
					overlap=true
				# If the transcript is fully before the gene
				elif [[ "$sstart" -lt "$begin" && "$send" -le "$begin" ]]
				then
					overlap=false
				# If the transcript is fully after the gene
				elif [[ "$sstart" -ge "$end" && "$send" -gt "$end" ]]
				then
					overlap=false
				fi

				if [[ "$overlap" = false ]]
				then
					if [[ "$(grep -c "$saccver	manual	gene	$sstart	$send" ${filename}.novelhits.gff)" -eq 0 ]]
					then
						echo -e "$saccver\tmanual\tgene\t$sstart\t$send\t.\t$strand\t.\tID=manual-${saccver}-blast-gene-0.${count};idenity=$pident;mismatches=$mismatch;coverage=$qcovs" >> ${filename}.novelhits.gff
						echo -e "$saccver\tmanual\tmRNA\t$sstart\t$send\t.\t$strand\t.\tID=manual-${saccver}-blast-gene-0.${count}-mRNA-1;Parent=manual-${saccver}-blast-gene-0.${count};idenity=$pident;mismatches=$mismatch;coverage=$qcovs" >> ${filename}.novelhits.gff
						count=$((count+1))
					fi
				fi
			done
		fi

#	done < $blast
	done < <(awk 'BEGIN{OFS="\t"}{if($3>=90 && $13>=70) print}' $blast)
done
