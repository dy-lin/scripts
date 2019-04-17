#!/bin/bash

if [[ "$#" -ne 2 ]]
then
	echo "$(basename $0) <spruce genotype> <file>" 1>&2
	exit 1
fi

spruce=$1
file=$2

if [[ "$spruce" == "all" && "$file" == "all" ]]
then
	trees="WS77111 PG29 Q903"
	types="Scaffolds Transcripts Proteins GFF Map"
	for i in $trees
	do
		echo "$i:"
		for j in $types
		do
			echo -ne "> $j:\t"
			spruce.sh $i $j
		done | column -t -s$'\t'
		echo
	done
	exit 0
fi

case $file in 
	[Ss]caf*)
		case $spruce in 
			[Ww][Ss]77111)
				echo "/projects/spruceup/pglauca/WS77111/assemblies/releases/version2/WS77111v2_release/WS77111-v2_1000plus_LGs.fa"
				;;
			[Pp][Gg]29)
				echo "/projects/spruceup/interior_spruce/PG29/assemblies/releases/PG29-v5/PG29-v5_1000plus.fa"
				;;
			[Qq]903)
				echo "/projects/spruceup/psitchensis/Q903/assembly/releases/version1/Q903_v1_1000plus.fa"
				;;
			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903:\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29:\t\t"
				spruce.sh PG29 $file
				;;
			\?)
				echo "Invalid genotype." 1>&2
				exit 1
				;;
		esac
		;;
	[Tt]ranscript*) 
		case $spruce in 
			[Ww][Ss]77111) 
				echo "/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.all.maker.transcriptsLoweAED.rename.function.fasta"
				;;
			[Pp][Gg]29)
				echo "/projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29v5.all.maker.transcriptsLoweAED.renamed.function.fasta"
				;;
			[Qq]903)
				echo "/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.all.maker.transcriptsLoweAED.renamed.function.fasta"
				;;
			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903:\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29:\t\t"
				spruce.sh PG29 $file
				;;	
			\?)
				echo "Invalid genotype." 1>&2
				exit 1
				;;
		esac
		;;
	[Gg][Ff][Ff]*) 
		case $spruce in
			[Ww][Ss]77111)
				echo "/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.all.polished.genesLoweAED.function_domain.gff"
				;;
			[Pp][Gg]29)
				echo "/projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29.all.polished.genesLoweAED.function_domain.gff"
				;;
			[Qq]903)
				echo "/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.all.polished.genesLoweAED.function_domain.gff"
				;;
			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903:\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29:\t\t"
				spruce.sh PG29 $file
				;;	
			\?)
				echo "Invalid genotype" 1>&2
				exit 1
				;;
		esac
		;;
	[Mm]ap) 
		case $spruce in 
			[Ww][Ss]77111)
				echo "/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.CompleteFinal.map"
				;;
			[Pp][Gg]29)
				echo "/projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29.CompleteFinal.map"
				;;
			[Qq]903)
				echo "/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.CompleteFinal.map"
				;;
			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29\t\t"
				spruce.sh PG29 $file
				;;	
			\?)
				echo "Invalid genotype." 1>&2
				exit 1
				;;
		esac
		;;
	[Pp]rot*|[Aa]nnotation) 
		case $spruce in
			[Ww][Ss]77111)
				echo "/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.all.maker.proteinsLoweAED.rename.function.fasta"
				;;
			[Pp][Gg]29)
				echo "/projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29v5.all.maker.proteinsLoweAED.renamed.function.fasta"
				;;
			[Qq]903)
				echo "/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.all.maker.proteinsLoweAED.renamed.function.fasta"
				;;
			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903:\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29:\t\t"
				spruce.sh PG29 $file
				;;
		esac
		;;
	[A][Ll][Ll])
		echo -ne "Scaffolds:\t"
		spruce.sh $spruce scaffolds
		echo -ne "Transcripts:\t"
		spruce.sh $spruce transcripts
		echo -ne "Proteins:\t"
		spruce.sh $spruce proteins
		echo -ne "GFF:\t\t"
		spruce.sh $spruce gff
		echo -ne "Map:\t\t"
		spruce.sh $spruce map
		;;
	\?) 
		echo "Invalid file type." 1>&2; exit 1
		;;
esac
