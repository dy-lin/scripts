#!/bin/bash
PROGRAM=$(basename $0)
gethelp=false
scaffolds=false
transcripts=false
gff=false
while getopts :hs:t:g: opt
do
	case $opt in
		h) gethelp=true;;
		s) scaffolds=$OPTARG; if [ ! -s "$scaffolds" ]; then echo "$(basename $scaffolds) is an invalid file."1>&2;fi;exit 1;;
		t) transcripts=$OPTARG; if [ ! -s "$transcripts" ]; then echo "$(basename $transcripts) is an invalid file."1>&2;fi;exit 1;;
		g) gff=$OPTARG; if [ ! -s "$gff" ]; then echo "$(basename $gff) is an invalid file." 1>&2;fi;exit 1;;
		\?) echo "$PROGRAM: Invalid option $OPTARG" 1>&2; exit 1;;
	esac
done
shift $((OPTIND-1))

# Write help/error messages here
if [[ "$#" -ne 0 || "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM -s <SCAFFOLD FASTA> -t <TRANSCRIPT FASTA> -g <GFF file>" 1>&2 
	echo "DESCRIPTION: Fetches necessary files for IGV." 1>&2
	exit 1
fi

# if scaffolds, gffs and transcripts are given-- typical case or proteome jackhmmer spruce
if [[ "$scaffolds" != false && "$transcripts" != false && "$gff" != false  ]]
then
	# GET scaffolds, gffs and transcripts
	echo "Fetching scaffolds, transcripts and GFF files..."
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
			ln -sf ../scaffolds/${scaffold}.scaffold.fa
			ln -sf ../gffs/${scaffold}.gff
			cd ..
		fi
		seqtk subseq $transcripts <(echo $i) >> transcripts/${scaffold}.transcripts.fa
		grep $i $gff >> gffs/${scaffold}.gff
	done

	date=$(date | awk '{if($3<10) {print "0" $3 $2 $6} else {print $3 $2 $6}}')
	cd scaffolds
	cat *.scaffold.fa > ../all.scaffolds.${date}.fa
	cd ../transcripts
	cat *.transcripts.fa > ../all.transcripts.${date}.fa
	cd ../gffs
	cat *.gff > ../all.${date}.gff
	cd ..


# if for spruce transcriptome, where transcripts and tbl are given
elif [[ "$transcripts" != false && "$scaffolds" = false && "$gff" = false ]]
then
	# in this case, must MAKE a gff.
	echo "Fetching transcripts..." 1>&2
	echo "Making GFF..." 1>&2 
	mkdir -p transcripts
	mkdir -p gffs
	mkdir -p igv
	for i in $(awk '/^>/' jackhmmer-blast-hits.faa | awk -F ">" '{print $2}' | awk '{print $1}')
	do
		# Get transcript name
		transcript_name=$(echo $i | awk '{print $1}' | awk -F "_" '{print $2}' | awk -F ":" '{print $1}')

		# Get transcript sequence and deposit into the folder
		seqtk subseq $transcripts <(echo $transcript_name) > transcripts/${transcript_name}.transcript.fa
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
		echo -e "##gff-version 3\n##sequence-region $transcript_name 1 $length" > gffs/${transcript_name}.gff
		echo -e "$transcript_name\tORFfinder\tgene\t1\t$length\t.\t$strand\t.\tID=${transcript_name}" >> gffs/${transcript_name}.gff
		echo -e "$transcript_name\tORFfinder\tmRNA\t1\t$length\t.\t$strand\t.\tID=${transcript_name}-mRNA-1;Parent=${transcript_name}" >> gffs/${transcript_name}.gff
		echo -e "$transcript_name\tORFfinder\texon\t$cds_begin\t$cds_end\t.\t$strand\t0\tID=${transcript_name}:exon;Parent=${transcript_name}-mRNA-1" >> gffs/${transcript_name}.gff
		echo -e "$transcript_name\tORFfinder\tCDS\t$cds_begin\t$cds_end\t.\t$strand\t0\tID=${transcript_name}:cds;Parent=${transcript_name}-mRNA-1" >> gffs/${transcrip_name}.gff 
		cd igv
		ln -sf ../transcripts/${transcript_name}.transcript.fa
		ln -sf ../gffs/${transcript_name}.gff
		cd ..
	done
	date=$(date | awk '{if($3<10) {print "0" $3 $2 $6} else {print $3 $2 $6}}')
	cd transcripts
	cat *.transcript.fa > ../all.transcripts.${date}.fa
	cd ../gffs
	cat *.gff > ../all.${date}.gff
	cd ..
elif [[ "$transcripts" != false && "$scaffolds" = false && "$gff" != false ]]
then
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
else
	echo "Invalid combination of options." 1>&2
	exit 1
fi
