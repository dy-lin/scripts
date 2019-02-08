#!/bin/bash
set -e -o pipefail
PROGRAM=$(basename $0)
# Experimental -- only FASTQ function tested.
bam=false
fq=false
gethelp=false
delete=false
if [[ "$#" -lt 1 ]]
then
	echo "USAGE: $PROGRAM [OPTIONS] <DIRECTORY>" 1>&2
	echo "DESCRIPTION: Takes a directory, finds all specified files and compresses them." 1>&2
	echo -e "OPTIONS: At least ONE of -b or -f  has to be selected.\n\t-b\tCompresses BAM files\n\t-d\tDeletes original file\n\t-f\tCompresses FASTQ files\n\t-h\tShow help menu" 1>&2

	exit 1
fi


while getopts :bfh opt
do
	case $opt in
		b) bam=true;;
		d) delete=true;;
		f) fq=true;;
		h) gethelp=true;;
		\?) echo "$PROGRAM: Error: Invalid option $opt" 1>&2; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM [OPTIONS] <DIRECTORY>" 1>&2
	echo "DESCRIPTION: Takes a directory, finds all specified files and compresses them." 1>&2
	echo -e "OPTIONS: At least ONE of -b or -f  has to be selected.\n\t-b\tCompresses BAM files\n\t-d\tDeletes original file\n\t-f\tCompresses FASTQ files\n\t-h\tShow help menu" 1>&2

	exit 1
fi

if [[ "$bam" = false && "$fq" = false ]]
then
	echo "ERROR: One of file types must be selected." 1>&2
	echo "USAGE: $PROGRAM <DIRECTORY>" 1>&2
	echo "DESCRIPTION: Takes a directory, finds all specified files and compresses them." 1>&2
	echo -e "OPTIONS: At least ONE of -b or -f has to be selected.\n\t-b\tCompresses BAM files\n\t-d\tDeletes original file\n\t-f\tCompresses FASTQ files\n\t-h\tShow help menu" 1>&2
	exit 1
fi
if [ ! -z "$1" ]
then
	if [[ "$fq" = true ]]
	then
		for file in $(find $1 -name "*.fastq")
		do
			if [[ ! -z "$file" ]]
			then
				dsrc c -t48 $file ${file}.dsrc
				if [[ "$delete" = true ]]
				then
					rm $file
				fi
			else
				break
			fi
		done
		
		for file in $(find $1 -name "*.fq")
		do
			if [[ ! -z "$file" ]]
			then
				dsrc c -t48 $file ${file}.dsrc
				if [[ "$delete" = true ]]
				then
					rm $file
				fi

			else
				break
			fi
		done

		for file in $(find $1 -name "*.fq.gz")
		do
			if [[ ! -z "$file" ]]
			then
				zcat $file | dsrc c -t48 -s ${file}.dsrc
				if [[ "$delete" = true ]]
				then
					rm $file
				fi

			else
				break
			fi
		done

		for file in $(find $1 -name "*.fastq.gz")
		do
			if [[ ! -z "$file" ]]
			then
				zcat $file | dsrc c -t48 -s ${file}.dsrc
				if [[ "$delete" = true ]]
				then
					rm $file
				fi

			else
				break
			fi
		done
	fi
	if [[ "$bam" = true ]]
	then
		for file in $(find $1 -name "*.bam")
		do
			if [[ ! -z "$file" ]]
			then
				filename=$(basename $file ".bam")
				dirname=$(dirname $file)
				fasta=${dirname}/${filename}.fa
				while [ ! -e "$fasta" ]
				do
					echo "Reference FASTA file not found." 1>&2
					echo "FASTA files in the current directory:" 1>&2
					cat <(ls ${dirname}/*.fa 2> /dev/null) <(ls ${dirname}/*.fasta 2> /dev/null) <(ls ${dirname}/*.fna 2> /dev/null) 1>&2
					echo "Please enter the full absolute path of your reference FASTA file:" 1>&2
					read fasta
				done

				fasta2ref --fasta-file $fasta --ref-file ${fasta}.fa.ref
				bam2spec --in $file --out ${file}.spec --ref ${fasta}.fa.ref

			else
				break
			fi
		done
	fi
else
	if [[ "$fq" = true ]]
	then
		for file in $(find . -name "*.fastq")
		do
			if [[ ! -z "$file" ]]
			then
				dsrc c -t48 $file ${file}.dsrc
				if [[ "$delete" = true ]]
				then
					rm $file
				fi

			else
				break
			fi
		done

		for file in $(find . -name "*.fq")
		do
			if [[ ! -z "$file" ]]
			then
				dsrc c -t48 $file ${file}.dsrc
				if [[ "$delete" = true ]]
				then
					rm $file
				fi

			else
				break
			fi
		done

		for file in $(find . -name "*.fq.gz")
		do
			if [[ ! -z "$file" ]]
			then
				zcat $file | dsrc c -t48 -s ${file}.dsrc
				if [[ "$delete" = true ]]
				then
					rm $file
				fi

			else
				break
			fi
		done

		for file in $(find . -name "*.fastq.gz")
		do
			if [[ ! -z "$file" ]]
			then
				zcat $file | dsrc c -t48 -s ${file}.dsrc
				if [[ "$delete" = true ]]
				then
					rm $file
				fi

			else
				break
			fi
		done

	fi
	if [[ "$bam" = true ]]
	then
		for file in $(find . -name "*.bam")
		do
			if [[ ! -z "$file" ]]
			then
				filename=$(basename $file ".bam")
				dirname=$(dirname $file)
				fasta=${dirname}/${filename}.fa
				while [ ! -e "$fasta" ]
				do
					echo "Reference FASTA file not found." 1>&2
					echo "FASTA files in the current directory:" 1>&2
					cat <(ls ${dirname}/*.fa 2> /dev/null) <(ls ${dirname}/*.fasta 2> /dev/null) <(ls ${dirname}/*.fna 2> /dev/null) 1>&2
					echo "Please enter the full absolute path of your reference FASTA file:" 1>&2
					read fasta
				done

				fasta2ref --fasta-file $fasta --ref-file ${fasta}.fa.ref
				bam2spec --in $file --out ${file}.spec --ref ${fasta}.fa.ref


			else
				break
			fi
		done
	fi
fi
