#!/bin/bash
set -eu -o pipefail
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
			[Ww][Ss]77111 | [Ww][Ss] | 77111)
				echo "/projects/spruceup/pglauca/WS77111/assemblies/releases/version2/WS77111v2_release/WS77111-v2_1000plus_LGs.fa"
				;;
			[Pp][Gg]29 | [Pp][Gg] | 29)
				echo "/projects/spruceup/interior_spruce/PG29/assemblies/releases/PG29-v5/PG29-v5_1000plus.fa"
				;;
			[Qq]903 | [Qq] | 903)
				echo "/projects/spruceup/psitchensis/Q903/assembly/releases/version1/Q903_v1_1000plus.fa"
				;;
			[Ss][Ee]404-851|[Ss][Ee]|[Ss][Ee]404|404-851|851)
				# Preparing for Engelmann files.
				;;
			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903:\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29:\t\t"
				spruce.sh PG29 $file
				# echo -ne "Se404-851:\t\t"
				# spruce.sh Se404-851 $file
				;;
			\?)
				echo "Invalid genotype." 1>&2
				exit 1
				;;
		esac
		;;
	[Tt]ranscript*) 
		case $spruce in 
			[Ww][Ss]77111 | [Ww][Ss] | 77111) 
				echo "/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.all.maker.transcriptsLoweAED.rename.function.fasta"
				;;
			[Pp][Gg]29 | [Pp][Gg] | 29)
				echo "/projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29v5.all.maker.transcriptsLoweAED.renamed.function.fasta"
				;;
			[Qq]903 | [Qq] | 903)
				echo "/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.all.maker.transcriptsLoweAED.renamed.function.fasta"
				;;
			[Ss][Ee]404-851|[Ss][Ee]|[Ss][Ee]404|404-851|851)
				# Preparing for Engelmann files.
				;;

			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903:\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29:\t\t"
				spruce.sh PG29 $file
				# echo -ne "Se404-851:\t\t"
				# spruce.sh Se404-851 $file
				;;	
			\?)
				echo "Invalid genotype." 1>&2
				exit 1
				;;
		esac
		;;
	[Gg][Ff][Ff]*) 
		case $spruce in
			[Ww][Ss]77111 | [Ww][Ss] | 77111)
				echo "/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.all.polished.genesLoweAED.function_domain.gff"
				;;
			[Pp][Gg]29 | [Pp][Gg] | 29)
				echo "/projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29.all.polished.genesLoweAED.function_domain.gff"
				;;
			[Qq]903 | [Qq] | 903)
				echo "/projects/btl/dlin/scratch/spruce/defensins/Q903.all.polished.genesLoweAED.function_domain.scaffoldrenamed.gff"
				# echo "/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.all.polished.genesLoweAED.function_domain.gff"
				;;
			[Ss][Ee]404-851|[Ss][Ee]|[Ss][Ee]404|404-851|851)
				# Preparing for Engelmann files.
				;;

			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903:\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29:\t\t"
				spruce.sh PG29 $file
				# echo -ne "Se404-851:\t\t"
				# spruce.sh Se404-851 $file
				;;	
			\?)
				echo "Invalid genotype" 1>&2
				exit 1
				;;
		esac
		;;
	[Mm]ap) 
		case $spruce in 
			[Ww][Ss]77111 | [Ww][Ss] | 77111)
				echo "/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.CompleteFinal.map"
				;;
			[Pp][Gg]29 | [Pp][Gg] | 29)
				echo "/projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29.CompleteFinal.map"
				;;
			[Qq]903 | [Qq] | 903)
				echo "/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.CompleteFinal.map"
				;;
			[Ss][Ee]404-851|[Ss][Ee]|[Ss][Ee]404|404-851|851)
					# Preparing for Engelmann files.
				;;

			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29\t\t"
				spruce.sh PG29 $file
				# echo -ne "Se404-851:\t\t"
				# spruce.sh Se404-851 $file
				;;	
			\?)
				echo "Invalid genotype." 1>&2
				exit 1
				;;
		esac
		;;
	[Pp]rot*|[Aa]nnotation) 
		case $spruce in
			[Ww][Ss]77111 | [Ww][Ss] | 77111)
				echo "/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.all.maker.proteinsLoweAED.rename.function.fasta"
				;;
			[Pp][Gg]29 | [Pp][Gg] | 29)
				echo "/projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29v5.all.maker.proteinsLoweAED.renamed.function.fasta"
				;;
			[Qq]903 | [Qq] | 903)
				echo "/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.all.maker.proteinsLoweAED.renamed.function.fasta"
				;;
			[Ss][Ee]404-851|[Ss][Ee]|[Ss][Ee]404|404-851|851)
				# Preparing for Engelmann files.
				;;

			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903:\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29:\t\t"
				spruce.sh PG29 $file
				# echo -ne "Se404-851:\t\t"
				# spruce.sh Se404-851 $file
				;;
		esac
		;;
	[Dd]efensin* | [Dd]ef)
		case $spruce in
			[Ww][Ss]77111 | [Ww][Ss] | 77111)
				echo "/projects/spruceup/scratch/pglauca/WS77111/annotation/amp/jackhmmer-proteome/defensins/final"
				# since this is a directory, ask to cd to it
				read -p "Go to this directory? " response
# 				case $response in
# 					[Yy]*) cd /projects/spruceup/scratch/pglauca/WS77111/annotation/amp/jackhmmer-proteome/defensins/final
# 						;;
# 					[Nn]*) # do nothing
# 						;;
# 				esac
 				;;
			[Pp][Gg]29 | [Pp][Gg] | 29)
				echo "/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/jackhmmer-proteome/defensins/final"

				read -p "Go to this directory? " response
# 				case $response in
# 					[Yy]*) cd /projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/jackhmmer-proteome/defensins/final
# 						;;
# 					[Nn]*) # do nothing
# 						;;
# 				esac
 				;;
			[Qq]903 | [Qq] | 903)
				echo "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/jackhmmer-proteome/defensins/final"
# 				read -p "Go to this directory? " response
# 				case $response in
# 					[Yy]*) cd /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/jackhmmer-proteome/defensins/final
# 						;;
# 					[Nn]*) # do nothing
# 						;;
# 				esac
			;;
			[Ss][Ee]404-851|[Ss][Ee]|[Ss][Ee]404|404-851|851)
				# Preparing for Engelmann files.
				;;

			[Aa][Ll][Ll])
				echo -ne "WS77111:\t"
				spruce.sh WS77111 $file
				echo -ne "Q903:\t\t"
				spruce.sh Q903 $file
				echo -ne "PG29:\t\t"
				spruce.sh PG29 $file
				# echo -ne "Se404-851:\t\t"
				# spruce.sh Se404-851 $file
				;;
		esac
		;;
	[Kk]allisto|[Kk]all|[Kk]al|[Kk])
		case $spruce in
			[Ww][Ss]77111| [Ww][Ss] | 77111)
				echo "/projects/spruceup/scratch/pglauca/WS77111/annotation/amp/kallisto"
# 				read -p "Go to this directory? " response
# 				case $response in
# 					[Yy]*) cd /projects/spruceup/scratch/pglauca/WS77111/annotation/amp/kallisto
# 						;;
# 					[Nn]*) # do nothing
# 						;;
# 				esac
 				;;
			[Pp][Gg]29 | [Pp][Gg] | 29)
				echo "/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
# 				read -p "Go to this directory? " response
# 				case $response in
# 					[Yy]*) cd /projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto
# 						;;
# 					[Nn]*) # do nothing
# 						;;
# 				esac
 				;;
			[Qq]903|[Qq]|903)
				echo "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
# 				read -p "Go to this directory? " response
# 				case $response in
# 					[Yy]*) cd /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto
# 						;;
# 					[Nn]* ) # do nothing
# 						;;
# 				esac
 				;;
			[Ss][Ee]404-851|[Ss][Ee]|[Ss][Ee]404|404-851|851)
				# Preparing for Engelmann files
				;;
		esac
		;;
	[Ii][Pp][Ss]|[Ii]nter[Pp]ro[Ss]can)
		case $spruce in
			[Ww][Ss]77111| [Ww][Ss] | 77111)
				echo "/projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/InterproScanFuncAnnotation/FourthStepMakerCumulatNonntEditRenamedLongIso/WS77111.all.maker.proteinsLoweAED.renamed.functionAnnotatedFinal.fasta.tsv"
				;;
			[Pp][Gg]29 | [Pp][Gg] | 29)
				echo "/projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/InterproScan/ThirdIterationCumulatNonntEditRenamedLongIso/PG29.all.maker.proteinsLoweAED.renamed.functionAnnotatedFinal.fasta.tsv"
				;;
			[Qq]903|[Qq]|903)
				echo "/projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/InterProScan/ThirdSTepCumulatRenamedLongIso/Q903.all.maker.proteinsLoweAED.renamed.functionAnnotatedFinal.fasta.tsv"
				;;
			[Ss][Ee]404-851|[Ss][Ee]|[Ss][Ee]404|404-851|851)
				;;
		esac
		;;
	[Gg][Mm][Aa][Pp]\ [Ii][Nn][Dd][Ee][Xx] | [Gg][Mm][Aa][Pp] | [Ii][Nn][Dd][Ee][Xx] | [Ii][Dd][Xx])
		case $spruce in
			[Ww][Ss]77111| [Ww][Ss] | 77111)
				echo "/projects/spruceup/pglauca/WS77111/assemblies/releases/version2/WS77111v2_release/GMAPindex/"
				;;
			 [Pp][Gg]29 | [Pp][Gg] | 29)
				echo "/projects/spruceup/interior_spruce/PG29/assemblies/releases/PG29-v5/PG29-v5_1000plus_GMAP-index/" 
				 ;;
			 [Qq]903|[Qq]|903)
				 echo "/projects/spruceup/psitchensis/Q903/assembly/releases/v0.1beta/gnavigator/Q903_v0.1beta_1000plus-gmap-index-dir/Q903_v0.1beta_1000plus-gmap-index/"
				 ;;
			 [Ss][Ee]404-851|[Ss][Ee]|[Ss][Ee]404|404-851|851)
				 ;;
		 esac
		;;
	[Aa][Ll][Ll])
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
