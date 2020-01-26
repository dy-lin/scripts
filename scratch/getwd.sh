#!/bin/bash

dir=$(readlink -f $1)

for i in $(ls -dtr $dir/*)
do
	echo "Working directory:"
	echo
	echo "{noformat}"
	echo $i
	echo "{noformat}"
	echo
done
