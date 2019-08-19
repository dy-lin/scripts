#!/bin/bash

gff=$1

while IFS=$'\n' read line
do
	gene=$(echo "$line" | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//')
	scaffold=$(echo "$line" | awk -F "\t" '{print $1}')
	echo -e "$gene\t$scaffold" 
done < <(awk '/\tgene\t/' $gff) | sort 
