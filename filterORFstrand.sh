#!/bin/bash

for fasta in $@
do
	for header in $(grep '^>' $fasta)
	do
		begin=$(echo $header | awk -F "_" '{print $2}' | awk -F ":" '{print $2}')
		end=$(echo $header | awk -F "_" '{print $2}' | awk -F ":" '{print $3}')

		if [[ "$begin" -gt "$end" ]]
		then
			seqtk subseq $fasta <(echo ${header:1}) >> ${fasta::-2}.neg.fa
		else
			seqtk subseq $fasta <(echo ${header:1}) >> ${fasta::-2}.pos.fa
		fi
	done
done


