#!/bin/bash

if [[ "$#" -ne 2 && "$#" -ne 3  ]]
then
	echo "USAGE: $(basename $0) <AMP class> <peptide length>" 1>&2
	echo "DESCRIPTION: Reads jackhmmer-blast-hits.faa and labels each sequence according to user." 1>&2
	exit 1
fi
amp=$1
length=$2

seqtk subseq jackhmmer-blast-hits.faa <(awk -v var=$length '{if($2==var) print $1}' <(seqtk comp jackhmmer-blast-hits.faa)) > jackhmmer-blast-hits.trimmed.faa

while read -u 3 id
do
	read -u 3 seq
	echo $id
	echo $seq
	while true
	do
		echo -e "\ttrb\ttrim before\n\ttra\ttrim after\n\ttrba\ttrim before and after\n\tp\tpartial\n\tl\t$amp-like\n\ttrbp\ttrim before and partial\n\ttrap\ttrim after and partial\n\tnp\tnon-partial\n\tap\talready partial\n\tnptra\tnon-partial, trimmed after\n\tnptrb\tnon-partial, trimmed before" 1>&2
		echo -n "Answer: " 1>&2
		read answer
	
		if [[ "$answer" == "trb" ]]
		then
			echo -n "Trim BEFORE this pattern: " 1>&2
			read pattern
			echo "$id, trimmed" >> jackhmmer-blast-hits.trimmed.faa
			echo $seq | sed "s/^.\+$pattern/$pattern/" >> jackhmmer-blast-hits.trimmed.faa
			break
		elif [[ "$answer" == "tra" ]]
		then
			echo -n "Trim AFTER this pattern: " 1>&2
			read pattern
			echo "$id, trimmed end" >> jackhmmer-blast-hits.trimmed.faa
			echo $seq | sed "s/${pattern}.\+$/$pattern/" >> jackhmmer-blast-hits.trimmed.faa
			break
		elif [[ "$answer" == "trba" ]]
		then
			echo -n "Trim BEFORE this pattern: " 1>&2
			read pattern
			echo -n "Trim AFTER this pattern: " 1>&2
			read pattern2

			echo "$id, trimmed, trimmed end" >> jackhmmer-blast-hits.trimmed.faa
			echo $seq | sed "s/^.\+$pattern/$pattern/" | sed "s/${pattern2}.\+$/${pattern2}/" >> jackhmmer-blast-hits.trimmed.faa
			break
		elif [[ "$answer" == "p" ]]
		then
			echo "$id, partial" >> jackhmmer-blast-hits.trimmed.faa
			echo $seq >> jackhmmer-blast-hits.trimmed.faa
			break
		elif [[ "$answer" == "l" ]]
		then
			echo "$id, $amp-like" >> jackhmmer-blast-hits.trimmed.faa
			echo "$seq" >> jackhmmer-blast-hits.trimmed.faa
			break
		elif [[ "$answer" == "trbp" ]]
		then
			echo -n "Trim BEFORE this pattern: " 1>&2
			read pattern
			echo "$id, trimmed, partial" >> jackhmmer-blast-hits.trimmed.faa
			echo $seq | sed "s/^.\+$pattern/$pattern/" >> jackhmmer-blast-hits.trimmed.faa
			break
		elif [[ "$answer" == "trap" ]]
		then
			echo -n "Trim AFTER this pattern: " 1>&2
			read pattern
			echo "$id, trimmed end, partial" >> jackhmmer-blast-hits.trimmed.faa
			echo $seq | sed "s/${pattern}.\+$/$pattern/" >> jackhmmer-blast-hits.trimmed.faa
			break
		elif [[ "$answer" == "np" ]]
		then
			echo $id | sed 's/, partial$//' >> jackhmmer-blast-hits.trimmed.faa
			echo $seq >> jackhmmer-blast-hits.trimmed.faa
			break
		elif [[ "$answer" == "ap" ]]
		then
			echo $id >> jackhmmer-blast-hits.trimmed.faa
			echo $seq >> jackhmmer-blast-hits.trimmed.faa
			break
		elif [[ "$answer" == "nptrb" ]]
		then
			echo -n "Trim BEFORE this pattern: " 1>&2
			read pattern
			echo $id | sed 's/, partial//' | sed 's/$/, trimmed/' >> jackhmmer-blast-hits.trimmed.faa
			echo $seq | sed "s/^.\+$pattern/$pattern/" >> jackhmmer-blast-hits.trimmed.faa
			break
		elif [[ "$answer" == "nptra" ]]
		then
			echo -n "Trim AFTER this pattern: " 1>&2
			read pattern
			echo $id | sed 's/, partial//' | sed 's/$/, trimmed end/'  >> jackhmmer-blast-hits.trimmed.faa
			echo $seq | sed "s/${pattern}.\+$/$pattern/" >> jackhmmer-blast-hits.trimmed.faa
			break	
		else
			echo "Invalid answer!" 1>&2
			continue
		fi
	done
done 3< <(seqtk subseq jackhmmer-blast-hits.faa <(awk -v var=$length '{if($2 != var) print $1}' <(seqtk comp jackhmmer-blast-hits.faa)))

echo $(wc -l jackhmmer-blast-hits.faa)
echo $(wc -l jackhmmer-blast-hits.trimmed.faa)

