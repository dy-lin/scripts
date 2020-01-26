#!/bin/bash

set -uo pipefail

dir=/projects/amp/peptaid/sra_reads

anurans="rtemporaria xlaevis calboguttata rpipiens xborealis xallofraseri xlargeni pnigromaculatus xtropicalis omargaretae"
hymenoptera="aechinatior amellifera acerana cobscurior omonticola nvitripennis nvitripennis_venom pbarbatus trugatulus ccastaneus"

THREADS=8

for i in $anurans
do
	cd $dir/anurans/$i
	if [[ -e sra.log ]]
	then
		rm sra.log
	else
		touch sra.log
	fi

	num_done=0
	num_failed=0
	total=$(cat sra.txt | wc -l)
	while read line
	do
		if [[ "$(ls ${line}*.fastq.gz 2> /dev/null | wc -l)" -ne 0 ]]
		then
			echo "$line: Already existed." >> sra.log
			num_done=$((num_done+1))
		else
			fasterq-dump --split-files --threads $THREADS --progress --force $line &> ${line}.log 
			if [[ "$?" -eq 0 ]]
			then
				echo "$line: Done." >> sra.log
				num_done=$((num_done+1))
				rm /home/dlin/ncbi/public/sra/${line}.sra.cache 2> /dev/null
				for j in ${line}*.fastq
				do
					pigz -p $THREADS $j 
				done
			else
				echo "$line: Failed." >> sra.log
				num_failed=$((num_failed+1))
				rm /home/dlin/ncbi/public/sra/${line}.sra.cache 2> /dev/null
			fi
		fi
	done < <(awk -F "\t" '{print $1}' sra.txt)
	echo "Done: $num_done/$total" >> sra.log
	echo "Failed: $num_failed/$total" >> sra.log
done

for i in $hymenoptera
do
	cd $dir/hymenoptera/$i
	if [[ -e sra.log ]]
	then
		rm sra.log
	else
		touch sra.log
	fi
	num_done=0
	num_failed=0
	total=$(cat sra.txt | wc -l)
	while read line
	do
		if [[ "$(ls ${line}*.fastq.gz 2> /dev/null | wc -l)" -ne 0  ]]
		then
			echo "$line: Already existed." >> sra.log
			num_done=$((num_done+1))
		else
			fasterq-dump --split-files --threads $THREADS --progress --force $line &> ${line}.log
			if [[ "$?" -eq 0 ]]
			then
				echo "$line: Done." >> sra.log
				rm /home/dlin/ncbi/public/sra/${line}.sra.cache 2> /dev/null
				
				num_done=$((num_done+1))

				for j in ${line}*.fastq
				do
					pigz -p $THREADS $j
				done
			else
				echo "$line: Failed." >> sra.log
				num_failed=$((num_failed+1))
				rm /home/dlin/ncbi/public/sra/${line}.sra.cache 2> /dev/null
				rm -rf fast.tmp 2> /dev/null
			fi
		fi
	done < <(awk -F "\t" '{print $1}' sra.txt)
	echo "Done: $num_done/$total" >> sra.log
	echo "Failed: $num_failed/$total" >> sra.log

done
