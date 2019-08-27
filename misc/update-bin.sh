#!/bin/bash

PROGRAM=$(basename $0)

for dir in $(ls -d /projects/btl/dlin/scripts/*)
do
	ln -sf $(grep -v "README" <(ls $dir/*) | tr "\n" " ") /projects/btl/dlin/bin/
done


