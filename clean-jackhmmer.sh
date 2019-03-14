#!/bin/bash

for i in guide-*
do
	if [[ -e "$i" ]]
	then
		rm $i
	fi
done

for i in jackhmmer_bs*
do
	if [[ -e "$i" ]]
	then
		rm $i
	fi
done
