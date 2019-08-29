#!/bin/bash
if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $(basename $0) <genotype>" 1>&2
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
			unmapped=$(samtools view -@ 128 -f4 $i | awk -F "\t" '{print $1}' | sort -u)
			seqtk subseq ../reads/${t}.R1.fq.gz <(echo "$unmapped") > ${filename}.unmapped.R1.fq
			seqtk subseq ../reads/${t}.R2.fq.gz <(echo "$unmapped") > ${filename}.unampped.R2.fq
			# align these reads to the assembly
			# should already be indexed

			seqtk subseq ../transcriptomes/${t}.fa <(bwa mem -t128 ../transcriptomes/${t}.fa ${filename}.unmapped.R1.fq ${filename}.unmapped.R2.fq | awk -F "\t" '{print $3}') > Q903.unmapped.to.mapped.${t}.transcripts.fa

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
		for i in Q903.annotation.${t}.?.sorted.bam
		do
			filename=$(basename $i ".sorted.bam")
			samtools fastq -@ 128 -f4 -1 ${filename}.unmapped.R1.fq -2 ${filename}.unmapped.R2.fq $i
			num=$(echo "$filename" | awk -F "." '{print $4}')

			# get the name of the unmapped reads, feed into seqtk to get those reads
			unmapped=$(samtools view -@ 128 -f4 $i | awk -F "\t" '{print $1}' | sort -u)
			seqtk subseq ../reads/${t}-${num}_R1.fq.gz <(echo "$unmapped") > ${filename}.unmapped.R1.fq
			seqtk subseq ../reads/${t}-${num}_R2.fq.gz <(echo "$unmapped") > ${filename}.unmapped.R2.fq
			# align these reads to the assembly
			# should already be indexed

			seqtk subseq .../transcriptomes/${t}.fa <(bwa mem -t128 ../transcriptomes/${t}.fa ${filename}.unmapped.R1.fq ${filename}.unmapped.R2.fq | awk -F "\t" '{print $3}') > Q903.unmapped.to.mapped.${t}.${num}.transcripts.fa
		done

	done

	tissues="Cort_Par Dev_SC"
	for t in $tissues
	do
		# get unmapped reads in annotation
		for i in Q903.annotation.${t}.?.sorted.bam
		do
			filename=$(basename $i ".sorted.bam")
			num=$(echo "$filename" | awk -F "." '{print $4}')

			# get the name of the unmapped reads, feed into seqtk to get those reads
			unmapped=$(samtools view -@ 128 -f4 $i | awk -F "\t" '{print $1}' | sort -u)
			seqtk subseq <(cat ../reads/${t}-${num}-00?_R1.fq.gz) <(echo "$unmapped") > ${filename}.unmapped.R1.fq
			seqtk subseq <(cat ../reads/${t}-${num}-00?-R2.fq.gz) <(echo "$unmapped") > ${filename}.unmapped.R2.fq
			# align these reads to the assembly
			# should already be indexed

			seqtk subseq .../transcriptomes/${t}.fa <(bwa mem -t128 ../transcriptomes/${t}.fa ${filename}.unmapped.R1.fq ${filename}.unmapped.R2.fq | awk -F "\t" '{print $3}') > Q903.unmapped.to.mapped.${t}.${num}.transcripts.fa
		done

	done

fi
