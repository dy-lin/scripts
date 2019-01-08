#!/bin/bash

if [ "$#" -ne 1 ]
then
	echo "USAGE: $(basename $0) <list of reads>"
	exit 1
fi

file=$1
total=$(sort -u $file | wc -l)

if [[ $(( $total % 2 )) -ne 0 ]]
then
	echo "Some reads are missing their pair."
	exit 1
else
	count=$((total/2))
fi
for line in $(head -n $count <(sort -u $file))
do
	if [[ $(echo $line | grep -c '_[0-9]_1_001[_.]') -gt 0 ]]
	then
		pair=$(echo $line | sed 's/_1_001/_2_001/')
		echo "$line $pair" >> reads.in
	elif [[ $(echo $line | grep -c '_R1[_.]') -gt 0 ]]
	then
		pair=$(echo $line | sed 's/_R1/_R2/')
		echo "$line $pair" >> reads.in
	elif [[ $(echo $line | grep -c '[_-]first') -gt 0 ]]
	then
		pair=$(echo $line | sed 's/first/second/')
		echo "$line $pair" >> reads.in
	else
		echo "Incompatible read names."
		exit 1
	fi
done
