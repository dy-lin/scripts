#!/bin/bash
set -euo pipefail

if [[ "$#" -eq 1 ]]
then
	count=$(grep -n $1 /projects/amp/peptaid/sra_reads/*/*/sra.txt | awk -F ":" '{print $2}')
	dir=$(grep -n $1 /projects/amp/peptaid/sra_reads/*/*/sra.txt | awk -F ":" '{print $1}')
	total=$(cat $dir | wc -l)
	echo "$(dirname $dir | awk -F "/" '{print $NF}'): $count/$total"
else
	pid=$(top -u dlin -n1 | grep -E 'fasterq-dump | pigz' | awk '{print $1}')

	wd=$(pwdx $pid | awk '{print $2}')

	acc=$(ps $pid | tail -n1 | awk '{print $NF}' | awk -F "_" '{print $1}')

	total=$(cat $wd/sra.txt | wc -l)

	done=$(grep -n $acc $wd/sra.txt | awk -F ":" '{print $1}')

	echo "$(basename $wd): $done/$total"
fi
