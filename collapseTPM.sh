#!/bin/bash
PROGRAM=$(basename $0)
header=false
fieldgiven=false
filetype=tsv
separator="\t"
while getopts :hf:s opt
do
	case $opt in
		h) header=true
			;;
		F) fieldgiven=true; field=$OPTARG; if [[ "$field" != [0-9] ]];then echo "Field number must be an integer." 1>&2; exit 1 ;fi
			;;
		f) filetype=$OPTARG;
			case $filetype in
				tsv) separator="\t"
					;;
				csv) separator=","
					;;
				jira) separator="|"
					;;
				?) echo "Unsupported filetype."; exit 1
					;;
			esac
			;;
		?) echo "$PROGRAM: invalid option: $opt" 1>&2; exit 1
			;;
	esac
done

shift $((OPTIND-1))
if [[ "$#" -eq 0 ]]
then
	echo "$(basename $0) <(TSV) file(s)>" 1>&2
	echo "DESCRIPTON: Takes a Field # to collapse, and (TSV) files." 1>&2
	echo -e "OPTIONS:\n\t-h\t\t\theader = true (Default: false)\n\t-f <filetype>\t\tavailable filetypes: tsv, csv jira\n\t-F <field number>\tDefault: NF"
	exit 1
fi

for file in $@
do
	if [[ -e "$file" ]]
	then
		newfile=${file%.*}.collapsed.tsv
		# Grab first column and get unique entries to iterate over.
		if [[  "$header" = true ]]
		then
			# grab headers
			# add header detection (look for numbers in the header)??
			if [[ "$fieldgiven" = true ]]
			then
		#		echo "Header extrapolated"
				head -n1 $file | awk -v var="$separator" -v bar="$field" -F $separator 'BEGIN{OFS=var}{print $1, $bar}' > $newfile

			else
				head -n1 $file | awk -v var="$separator" -F $separator 'BEGIN{OFS=var}{print $1, $NF}' > $newfile
			fi

			for gene in $(awk -F "\t" '{print $1}' <(tail -n +2 $file) | sed 's/-R.\?//' |sort -u)
			do
				if [[ "$fieldgiven" = false ]]
				then
					tpm=$(awk -v var="$gene" -F "\t" '$0 ~ var {sum+=$NF} END {if (NR !=0) {print sum} else { print "There is nothing to sum."; exit 1 }}' $file )
				else
					tpm=$(awk -v var="$gene" -v bar="$field" -F "\t" '$0 ~ var {sum+=$bar} END {if (NR !=0) {print sum} else { print "There is nothing to sum."; exit 1 }}' $file )
				fi
		#		echo "TPM: $tpm"
				echo -e "$gene$separator$tpm" >> $newfile
			done
		else
			echo -n > $newfile	
			for gene in $(awk -F "\t" '{print $1}' $file| sed 's/-R.\?//' |sort -u)
			do
				if [[ "$fieldgiven" = false ]]
				then
					tpm=$(awk -v var="$gene" -F $separator '$0 ~ var {sum+=$NF} END {if (NR !=0) {print sum} else { print "There is nothing to sum."; exit 1 }}' $file )
				else
					tpm=$(awk -v var="$gene" -v bar="$field" -F $separator '$0 ~ var {sum+=$bar} END {if (NR !=0) {print sum} else { print "There is nothing to sum."; exit 1 }}' $file )
				fi
			#	echo "TPM: $tpm"
				echo -e "$gene$separator$tpm" >> $newfile
			done
		fi
	else
		echo "$(basename $file) does not exist." 1>&2
	fi
done

