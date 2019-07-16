#!/bin/bash

if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $(basename $0) <SV file(s)>" 1>&2
	echo "DESCRIPTION: Turns the header of a SV file into the first column and vice versa." 1>&2
	exit 1
fi

for file in $@
do
	ext=${file##*.}
	output=${file%.*}
	outfile=${output}.flip.${ext}

	comma_count=$(grep -c ',' $file)
	tab_count=$(grep -cP '\t' $file)
	jira_count=$(grep -c '|' $file)

	most=$(sort -g -r -k2,2 <(echo -e "commas $comma_count\ntabs $tab_count\npipes $jira_count") | head -n1 | awk '{print $1}')
	# echo $most
	flip.py $file $most > $outfile
	if [[ "$?" -ne 0 ]]
	then
		if [[ -e "$outfile" ]]
		then
			rm $outfile
		fi
	fi
done


