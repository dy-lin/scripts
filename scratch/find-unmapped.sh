#!/bin/bash

PROGRAM=$(basename $0)
combined=false

while getopts :c opt
do
	case $opt in
		c) combined=true
			;;
		\?) echo "$PROGRAM: Invalid option: $OPTARG" 1>&2; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $(basename $0) [-c] <genotype>" 1>&2
	echo "DESCRIPTION: Finds the transcripts in the assembly that unmapped reads (in the annotation)  map to." 1>&2
	exit 1
fi

genotype=$1

if [[ "$genotype" == "PG29" || "$genotype" == "all" ]]
then
	cd /projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/alignment
	tissues="bark embryo flush_bud mature_needle megagametophyte seed_germination xylem young_buds"


	for t in $tissues
	do
		# get unmapped reads in annotation
		for i in PG29.annotation.${t}.sorted.bam
		do
			# get the name of the unmapped reads, feed into seqtk to get those reads
			echo "Getting unmapped reads..." 1>&2
			unmapped=$(samtools view -@ 128 -f4 $i | awk -F "\t" '{print $1}' | sort -u)

			# pull those unmapped reads out of the reads
			# cannot use samtools fastq -f4 since those will pull out all unmapped reads (excluding their pairs if those are mapped)
			# get them individually from each read pair file
			echo "Pulling the first unmapped read of each read pair..." 1>&2
			seqtk subseq ../reads/${t}.R1.fq.gz <(echo "$unmapped") > ${filename}.unmapped.R1.fq
			echo "Pulling the second unmapped read of each read pair..." 1>&2
			seqtk subseq ../reads/${t}.R2.fq.gz <(echo "$unmapped") > ${filename}.unampped.R2.fq
			# align these reads to the assembled transcriptome

			echo "Aligning these unmapped read pairs to the transcriptome assembly..." 1>&2
			# should already be indexed in previous script 
			bwa mem -t128 ../transcriptomes/${t}.fa ${filename}.unmapped.R1.fq ${filename}.unmapped.R2.fq | samtools sort | samtools view -b -o ${filename}.all.unmapped.sorted.bam
			echo "Extracting the transcripts that the reads mapped to..." 1>&2
			seqtk subseq ../transcriptomes/${t}.fa <(samtools view ${filename}.all.unmapped.sorted.bam | awk -F "\t" '{print $3}' | sort -u) > PG29.unmapped.to.mapped.${t}.transcripts.fa

		done
	done
fi


if [[ "$genotype" == "Q903" || "$genotype" == "all" ]]
then
	cd /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/alignment

 	tissues="control gallery wound"
 
 
 	for t in $tissues
 	do
 		# get unmapped reads in annotation
		if [[ "$combined" = false ]]
		then
			for i in Q903.annotation.${t}.?.sorted.bam
			do
				filename=$(basename $i ".sorted.bam")
				# samtools fastq does not pull out unmapped reads and their pairs
			#	samtools fastq -@ 128 -f4 -1 ${filename}.unmapped.R1.fq -2 ${filename}.unmapped.R2.fq $i
				num=$(echo "$filename" | awk -F "." '{print $4}')
	 
				# get the name of the unmapped reads, feed into seqtk to get those reads
				echo "Getting unmapped reads..." 1>&2
				unmapped=$(samtools view -@ 128 -f4 $i | awk -F "\t" '{print $1}' | sort -u)
				echo "Pulling the first unmapped read of each read pair..." 1>&2
				seqtk subseq ../reads/${t}-${num}_R1.fq.gz <(echo "$unmapped") > ${filename}.unmapped.R1.fq
	 
				echo "Pulling the second unmapped read of each read pair..." 1>&2
				seqtk subseq ../reads/${t}-${num}_R2.fq.gz <(echo "$unmapped") > ${filename}.unmapped.R2.fq
				# align these reads to the assembly
				# should already be indexed
				
				echo "Aligning these unmapped read pairs to the transcriptome assembly..." 1>&2
				bwa mem -t128 ../transcriptomes/${t}.fa ${filename}.unmapped.R1.fq ${filename}.unmapped.R2.fq | samtools sort | samtools view -b -o ${filename}.all.unmapped.sorted.bam
				echo "Extracting the transcripts that the reads mapped to..." 1>&2
				seqtk subseq ../transcriptomes/${t}.fa <(samtools view ${filename}.all.unmapped.sorted.bam | awk -F "\t" '{print $3}' | sort -u) > Q903.unmapped.to.mapped.${t}.${num}.transcripts.fa
			done
		else
			i=Q903.annotation.${t}.sorted.bam
			filename=$(basename $i ".sorted.bam")
			echo "Getting unmapped reads..." 1>&2
			unmapped=$(samtools view -@ 128 -f4 $i | awk -F "\t" '{print $1}' | sort -u)
			echo "Pulling the first unmapped read of each read pair..." 1>&2 
			seqtk subseq <(cat ../reads/${t}-?_R1.fq.gz) <(echo "$unmapped") > ${filename}.unmapped.R1.fq

			echo "Pulling the second unmapped read of each read pair..." 1>&2
			seqtk subseq <(cat ../reads/${t}-?_R2.fq.gz) <(echo "$unmapped") > ${filename}.unmapped.R2.fq

			echo "Aligning these unmapped read pairs to the transcriptome assembly..." 1>&2
			bwa mem -t128 ../transcriptomes/${t}.fa ${filename}.unmapped.R1.fq ${filename}.unmapped.R2.fq | samtools sort | samtools view -b -o ${filename}.all.unmapped.sorted.bam
			echo "Extracting the transcripts that the reads mapped to..." 1>&2
			seqtk subseq ../transcriptomes/${t}.fa <(samtools view ${filename}.all.unmapped.sorted.bam | awk -F "\t" '{print $3}' | sort -u) > Q903.unmapped.to.mapped.${t}.transcripts.fa
		fi
 
 	done

	tissues="Cort_Par Dev_SC"
	for t in $tissues
	do
		# get unmapped reads in annotation
		if [[ "$combined" = false ]]
		then
			for i in Q903.annotation.${t}.?.sorted.bam
			do
				filename=$(basename $i ".sorted.bam")
				num=$(echo "$filename" | awk -F "." '{print $4}')

				# get the name of the unmapped reads, feed into seqtk to get those reads
				echo "Getting unmapped reads..." 1>&2
				unmapped=$(samtools view -@ 128 -f4 $i | awk -F "\t" '{print $1}' | sort -u)
				echo "Pulling the first unmapped read of each read pair..." 1>&2
				seqtk subseq <(cat ../reads/${t}-${num}-00?_R1.fq.gz) <(echo "$unmapped") > ${filename}.unmapped.R1.fq
				echo "Pulling the second unmapped read of each read pair..." 1>&2
				seqtk subseq <(cat ../reads/${t}-${num}-00?_R2.fq.gz) <(echo "$unmapped") > ${filename}.unmapped.R2.fq
				# align these reads to the assembly
				# should already be indexed
				
				echo "Aligning these unmapped read pairs to the transcriptome assembly..." 1>&2

				bwa mem -t128 ../transcriptomes/${t}.fa ${filename}.unmapped.R1.fq ${filename}.unmapped.R2.fq | samtools sort | samtools view -b -o ${filename}.all.unmapped.sorted.bam
				echo "Extracting the transcripts that the reads mapped to..." 1>&2
				seqtk subseq ../transcriptomes/${t}.fa <(samtools view ${filename}.all.unmapped.sorted.bam | awk -F "\t" '{print $3}' | sort -u) > Q903.unmapped.to.mapped.${t}.${num}.transcripts.fa
			done
		else
			i=Q903.annotation.${t}.sorted.bam
			filename=$(basename $i ".sorted.bam")
			echo "Getting unmapped reads..." 1>&2
			unmapped=$(samtools view -@ 128 -f4 $i | awk -F "\t" '{print $1}' | sort -u)
			echo "Pulling the first unmapped read of each read pair..." 1>&2
			seqtk subseq <(cat ../reads/${t}-${num}-00?_R1.fq.gz) <(echo "$unmapped") > ${filename}.unmapped.R1.fq
			echo "Pulling the second unmapped read of each read pair..." 1>&2
			seqtk subseq <(cat ../reads/${t}-${num}-00?_R2.fq.gz) <(echo "$unmapped") > ${filename}.unmapped.R2.fq

			echo "Aligning these unmapped read pairs to the transcriptome assembly..." 1>&2
			bwa mem -t128 ../transcriptomes/${t}.fa ${filename}.unmapped.R1.fq ${filename}.unmapped.R2.fq | samtools sort | samtools view -b -o ${filename}.all.unmapped.sorted.bam
			echo "Extracting the transcripts that the reads mapped to..." 1>&2
			seqtk subseq ../transcriptomes/${t}.fa <(samtools view ${filename}.all.unmapped.sorted.bam | awk -F "\t" '{print $3}' | sort -u) > Q903.unmapped.to.mapped.${t}.${num}.transcripts.fa
		fi
	done

fi
