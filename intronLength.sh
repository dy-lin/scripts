#!/bin/bash

gff=$1

IFS=$'\n'

for line in $(cat $gff)
do
	feature=$(echo $line | awk -F "\t" '{print $3}')
	if [ "$feature" == "mRNA" ]
	then
		name=$(echo $line | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//')
	fi
	if [ "$feature" == "intron" ]
	then
		begin=$(echo $line | awk -F "\t" '{print $4}')
		end=$(echo $line | awk -F "\t" '{print $5}')
		length=$((end-begin+1))
		echo -e "$name\t$length"
	fi
done
	
