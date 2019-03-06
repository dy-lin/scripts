#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)
gethelp=false
scaffolds=false
transcripts=false
gff=false
if [[ "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM -s <SCAFFOLD FASTA> -t <TRANSCRIPT FASTA> -g <GFF file>" 1>&2 
	echo "DESCRIPTION: Fetches necessary files for IGV." 1>&2
	exit 1
fi


while getopts :hs:t:g: opt
do
	case $opt in
		h) gethelp=true;;
		s) scaffolds=$OPTARG; if [ ! -s "$scaffolds" ]; then echo "$(basename $scaffolds) is an invalid file."1>&2;exit 1;fi;;
		t) transcripts=$OPTARG; if [ ! -s "$transcripts" ]; then echo "$(basename $transcripts) is an invalid file."1>&2;exit 1; fi;;
		g) gff=$OPTARG; if [ ! -s "$gff" ]; then echo "$(basename $gff) is an invalid file." 1>&2;exit 1; fi;;
		\?) echo "$PROGRAM: Invalid option $OPTARG" 1>&2; exit 1;;
	esac
done
shift $((OPTIND-1))
# Write help/error messages here
if [[ "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM -s <SCAFFOLD FASTA> -t <TRANSCRIPT FASTA> -g <GFF file>" 1>&2 
	echo "DESCRIPTION: Fetches necessary files for IGV." 1>&2
	exit 1
fi
# if scaffolds, gffs and transcripts are given-- typical case or proteome jackhmmer spruce
if [[ "$scaffolds" != false && "$transcripts" != false && "$gff" != false  ]]
then
	echo "Scaffolds, transcripts and GFF files detected!" 1>&2
	# GET scaffolds, gffs and transcripts
	echo "Fetching scaffolds, transcripts and GFF files..." 1>&2
	if [ ! -e "scaffolds" ]
	then
		mkdir scaffolds
	fi

	if [ ! -e "transcripts" ]
	then
		mkdir transcripts
	fi

	if [ ! -e "gffs" ]
	then
		mkdir gffs
	fi

	if [ ! -e "igv" ]
	then
		mkdir igv
	fi

	for i in $(awk '/^>/' jackhmmer-blast-hits.faa | awk -F ">" ' {print $2}' | awk '{print $1}')
	do
		# Get scaffold name
		temp=${i#*-}
		scaffold=${temp%-*-gene-*-mRNA-?}
		if [ ! -e "scaffolds/${scaffold}.scaffold.fa" ]
		then
			seqtk subseq $scaffolds <(echo $scaffold) > scaffolds/${scaffold}.scaffold.fa
			length=$(tail -n 1 "scaffolds/${scaffold}.scaffold.fa" | head -c -1 | wc -m)
			echo -e "##gff-version 3\n##sequence-region $scaffold 1 $length" > gffs/${scaffold}.gff
			cd igv
			# Create IGV softlinks
			ln -sf $(readlink -f ../scaffolds/${scaffold}.scaffold.fa)
			ln -sf $(readlink -f ../gffs/${scaffold}.gff)
			cd ..
		fi
		seqtk subseq $transcripts <(echo $i) >> transcripts/${scaffold}.transcripts.fa
		grep ${i::-7} $gff >> gffs/${scaffold}.gff
	done

	date=$(date | awk '{if($3<10) {print "0" $3 $2 $6} else {print $3 $2 $6}}')
	cd scaffolds
	cat *.scaffold.fa > ../all.scaffolds.raw.${date}.fa
	for file in $(ls ../all.scaffolds.*.fa)
	do
		if [[ "../all.scaffolds.raw.${date}.fa" -nt "$file" ]]
		then
			rm "$file"
		fi
	done
	cd ../transcripts
	cat *.transcripts.fa > ../all.transcripts.raw.${date}.fa
	for file in $(ls ../all.transcripts.*.fa)
	do
		if [[ "../all.transcripts.raw.${date}.fa" -nt "$file" ]]
		then
			rm "$file"
		fi
	done
	cd ../gffs
	cat *.gff > ../all.raw.${date}.gff
	for file in $(ls ../all.*.gff)
	do
		if [[ "../all.raw.${date}.gff" -nt "$file" ]]
		then
			rm "$file"
		fi
	done
	cd ..
	echo "...Done." 1>&2

# if for spruce transcriptome, where transcripts
elif [[ "$transcripts" != false && "$scaffolds" = false && "$gff" = false ]]
then
	echo "Transcripts detected!" 1>&2
	# in this case, must MAKE a gff.
	echo "Fetching transcripts..." 1>&2
	echo "Making GFFs..." 1>&2 
	mkdir -p transcripts
	mkdir -p gffs
	mkdir -p igv
	for i in $(awk '/^>/' jackhmmer-blast-hits.faa | awk -F ">" '{print $2}' | awk '{print $1}')
	do
		# Get transcript name
		transcript_name=$(echo $i | awk '{print $1}' | awk -F "_" '{print $2}' | awk -F ":" '{print $1}')

		# Get transcript sequence and deposit into the folder
		if [[ ! -e "transcripts/${transcript_name}.transcript.fa" ]]
		then
			seqtk subseq $transcripts <(echo $transcript_name) > transcripts/${transcript_name}.transcript.fa
		fi
		length=$(seqtk comp transcripts/${transcript_name}.transcript.fa | awk '{print $2}')


		cds_begin=$(echo $i | awk '{print $1}' | awk -F "_" '{print $2}' | awk -F ":" '{print $2}')
		cds_end=$(echo $i | awk '{print $1}' | awk -F "_" '{print $2}' | awk -F ":" '{print $3}')
		if [[ "$cds_begin" -gt "$cds_end" ]]
		then
			strand="-"
		else
			strand="+"
		fi
		
		# Start making GFF
		# First line of GFF
		if [[ ! -e "gffs/${transcript_name}.gff" ]]
		then
			echo -e "##gff-version 3\n##sequence-region $transcript_name 1 $length" > gffs/${transcript_name}.gff
		fi
		# remove gene annotation as it's not a gene -- intron has been removed
	#	echo -e "$transcript_name\tORFfinder\tgene\t1\t$length\t.\t$strand\t.\tID=${transcript_name}" >> gffs/${transcript_name}.gff
		echo -e "$transcript_name\tORFfinder\tmRNA\t1\t$length\t.\t$strand\t.\tID=${transcript_name}-mRNA-1;Parent=${transcript_name}" >> gffs/${transcript_name}.gff
		echo -e "$transcript_name\tORFfinder\texon\t$cds_begin\t$cds_end\t.\t$strand\t0\tID=${transcript_name}:exon;Parent=${transcript_name}-mRNA-1" >> gffs/${transcript_name}.gff
		echo -e "$transcript_name\tORFfinder\tCDS\t$cds_begin\t$cds_end\t.\t$strand\t0\tID=${transcript_name}:cds;Parent=${transcript_name}-mRNA-1" >> gffs/${transcript_name}.gff 
		cd igv
		ln -sf $(readlink -f ../transcripts/${transcript_name}.transcript.fa)
		ln -sf $(readlink -f ../gffs/${transcript_name}.gff)
		cd ..
	done
	if [[ -s "jackhmmer-blast-hits.trimmed.faa" ]]
	then
		amp=$(pwd | awk -F "/" ' {print $(NF-2)}')
		mkdir -p transcripts/${amp}-only/cds
		# GET AMP only
		for i in $(grep '^>' jackhmmer-blast-hits.trimmed.faa | grep -v "\-like" | awk -F "_" '{print $2}' | awk -F ":" '{print $1}' | sort -u )
		do
			cp transcripts/${i}.transcript.fa transcripts/${amp}-only/
		done
		
		# GET CDS
		for i in $(grep '^>' jackhmmer-blast-hits.trimmed.faa | grep -v '\-like' | awk '{print $1}' | awk -F "|" '{print $2}')
		do
			j=$(echo $i | awk -F "_" '{print $2}' | awk -F ":" '{print $1}')
			grep -A1 --no-group-separator $i ../*.cds.fa > transcripts/${amp}-only/cds/${j}.cds.fa
		done

		cat transcripts/${amp}-only/cds/*.cds.fa > transcripts/${amp}-only/all.cds.fa
	fi
	
	date=$(date | awk '{if($3<10) {print "0" $3 $2 $6} else {print $3 $2 $6}}')
	cd transcripts
	cat *.transcript.fa > ../all.transcripts.${date}.fa
	# if older files exist, remove them
	for file in $(ls ../all.transcripts.*.fa)
	do
		if [[ "../all.transcripts.${date}.fa" -nt "$file" ]]
		then
			rm "$file"
		fi
	done
	cd ../gffs
	cat *.gff > ../all.${date}.gff
	for file in $(ls ../all.*.gff)
	do
		if [[ "../all.${date}.gff" -nt "$file" ]]
		then
			rm "$file"
		fi
	done
	cd ..
	echo "...Done." 1>&2
elif [[ "$transcripts" != false && "$scaffolds" = false && "$gff" != false ]]
then
	echo "Transcripts and GFF detected!" 1>&2
	# necessary GFF headers are in line 1 and 7 of NCBI GFFs
	awk 'NR==1 || NR==7' $gff > amp.gff
	# NCBI data is named differently
	for i in $(awk '/^>/' jackhmmer-blast-hits.faa | awk -F ">" '{print $2}' | awk '{print $1}')
	do
		# No need to fetch scaffolds when genome is assembled into one scaffold
		# Fetch transcripts
		seqtk subseq $transcripts <(echo $i) >> amp.transcripts.fa
		# Fetch GFF
		grep $i $gff >> amp.gff
	done
	echo "....Done." 1>&2
else
	echo "Invalid combination of options." 1>&2
	echo -e "OPTIONS:\n\t-s\tscaffold FASTA file\n\t-t\ttranscript FASTA file\n\t-g\tGFF file" 1>&2
	exit 1
fi
