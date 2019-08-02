#!/bin/bash

if [[ "$#" -ne 2 ]]
then
	echo "USAGE: $(basename $0) <FASTA file> <FASTA file>" 1>&2
	exit 1
fi

# longer scaffold goes on top

length1=$(seqtk comp $1 | awk '{print $2}')
length2=$(seqtk comp $2 | awk '{print $2}')

if [[ "$length1" -ge "$length2" ]]
then
	top_scaffold=$1
	bottom_scaffold=$2
	top_length="$length1"
else
	top_scaffold=$2
	bottom_scaffold=$1
	top_length="$length2"
fi
if [[ "$(echo $top_scaffold | grep -c 'rc')" -eq 0 ]]
then
	top=$(echo $top_scaffold | sed 's/.scaffold.fa//')
else
	top=$(echo $top_scaffold | sed 's/.scaffold.rc.fa//').rc
fi

if [[ "$(echo $bottom_scaffold | grep -c 'rc')" -eq 0 ]]
then
	bottom=$(echo $bottom_scaffold | sed 's/.scaffold.fa//')
else
	bottom=$(echo $bottom_scaffold | sed 's/.scaffold.rc.fa//').rc
fi

rep=${bottom}_vs_${top}.rep
if [[ ! -e "$rep" ]]
then
	/home/pubseq/BioSw/phrap/current/cross_match $top_scaffold $bottom_scaffold -minmatch 10 -minscore 10 -masklevel 101 > $rep 2> cross_match.log
else
	read -p "File already exists. Run cross_match again? " response
	if [[ "$response" == [Yy]* ]]
	then

		/home/pubseq/BioSw/phrap/current/cross_match $top_scaffold $bottom_scaffold -minmatch 10 -minscore 10 -masklevel 101 > $rep 2> cross_match.log
	fi
fi

mismatch=10
block=69

alpha=200
rleap=10
filetype=png

scale=$((top_length/2070))


if [[ ! -e "xmv-${rep}_m${mismatch}_b${block}_r${rleap}_c${scale}.${filetype}" ]]
then
	xmatchview.py -x $rep -s $top_scaffold -q $bottom_scaffold -m $mismatch -b $block -c $scale -a $alpha -r $rleap -f $filetype -p /projects/btl/lcoombe/git/xmatchview/tarballs/fonts -e ${top}.tsv -y ${bottom}.tsv &> xmv.log
	echo "Saved: xmv-${rep}_m${mismatch}_b${block}_r${rleap}_c${scale}.${filetype}" 
else
	read -p "File already exists. Run XMatchView again? " response
	if [[ "$response" == [Yy]* ]]
	then

		xmatchview.py -x $rep -s $top_scaffold -q $bottom_scaffold -m $mismatch -b $block -c $scale -a $alpha -r $rleap -f $filetype -p /projects/btl/lcoombe/git/xmatchview/tarballs/fonts -e ${top}.tsv -y ${bottom}.tsv &> xmv.log

	echo "Saved: xmv-${rep}_m${mismatch}_b${block}_r${rleap}_c${scale}.${filetype}" 
	fi
fi




