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
	echo "USAGE: $(basename $0) [-c] <Genotype>" 1>&2
	echo "DESCRIPTION: Aligns RNAseq reads to the annotated or assembled transcriptome." 1>&2
	echo -e "OPTIONS:\n\t-c\tcombine all replicates"
	exit 1
fi

if [[ "$1" == "PG29" || "$1" == "all" ]]
then
	# PG29 first
	tissues="bark embryo flush_bud mature_needle megagametophyte seed_germination xylem young_buds"
	outdir="/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/alignment"
	mkdir -p $outdir

	echo -e "Transcriptome\tReads\t\# of Unmapped Reads\tPercentage (%)" > $outdir/PG29-map-rna.tsv

	transcriptome="/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/PG29.maker.transcripts.edited.fa"
	transcriptome_dir="/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/transcriptomes"
	
	if [[ ! -e "${transcriptome}.sa" ]]
	then
		echo "Indexing annotated PG29 transcriptome..." 1>&2
		bwa index $transcriptome 2>> $outdir/bwa-index.log
	fi

	for t in $tissues
	do
		total=$(countReads.sh /projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/PG29_reads/${t}.R?.fq.gz | tail -n1 | awk '{print $2}')
		if [[ ! -e "$outdir/PG29.annotation.${t}.sorted.bam" ]]
		then
			echo "Aligning PG29 $t reads to annotated transcriptome..." 1>&2
			bwa mem -t128 $transcriptome /projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/PG29_reads/${t}.R?.fq.gz 2>> $outdir/bwa-mem-${t}.log | samtools sort -@ 128 -O BAM -o $outdir/PG29.annotation.${t}.sorted.bam
		fi
		count=$(samtools view -f4 -c -@ 128 $outdir/PG29.annotation.${t}.sorted.bam)
		echo "Number of unmapped PG29 $t reads: $count/$total" 1>&2
		percentage=$(echo "scale=4;$count / $total * 100" | bc)
		echo "Percentage of unmapped PG29 $t reads: $percentage" 1>&2
		echo -e "Annotation\t$t\t$count\t$percentage" >> $outdir/PG29-map-rna.tsv

		if [[ ! -e "$transcriptome_dir/${t}.fa.sa" ]]
		then
			echo "Indexing assembled PG29 $t transcriptome..." 1>&2
			bwa index $transcriptome_dir/${t}.fa 2>> $outdir/bwa-index-${t}.log
		fi
		if [[ ! -e "$outdir/PG29.assembly.${t}.sorted.bam" ]]
		then
			echo "Aligning PG29 $t reads to assembled transcriptome..." 1>&2
			bwa mem -t128 $transcriptome_dir/${t}.fa /projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/PG29_reads/${t}.R?.fq.gz 2>> $outdir/bwa-mem-${t}.log | samtools sort -@ 128 -O BAM -o $outdir/PG29.assembly.${t}.sorted.bam
		fi
		count=$(samtools view -@ 128 -f4 -c $outdir/PG29.assembly.${t}.sorted.bam)
		echo "Number of unmapped PG29 $t reads: $count/$total"
		percentage=$(echo "scale=4;$count / $total * 100" | bc)
		echo "Percentage of unmapped PG29 $t reads: $percentage" 1>&2
		echo -e "Assembly\t$t\t$count\t$percentage" >> $outdir/PG29-map-rna.tsv

	done
fi

if [[ "$1" == "Q903" || "$1" == "all" ]]
then
	# Q903
	tissues="Dev_SC Cort_Par"
	outdir="/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/alignment"
	mkdir -p $outdir


	echo -e "Transcriptome\tReads\t\# of Unmapped Reads\tPercentage (%)" > $outdir/Q903-map-rna.tsv

	transcriptome="/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903.maker.transcripts.edited.fa"
	transcriptome_dir="/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/transcriptomes"
	if [[ ! -e "${transcriptome}.sa" ]]
	then
		echo "Indexing annotated Q903 transcriptome..." 1>&2
		bwa index $transcriptome 2>> $outdir/bwa-index.log
	fi
	for t in $tissues
	do
		if [[ "$combined" = false ]]
		then
			for num in $(seq 4 2 8)
			do
				echo "Counting the total number of reads..." 1>&2
				total=$(countReads.sh /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}-00?_R?.fq.gz | tail -n1 | awk '{print $2}')
				if [[ ! -e "$outdir/Q903.annotation.${t}.${num}.sorted.bam" ]]
				then
					echo "Aligning Q903 $t reads to annotated transcriptome..." 1>&2
					bwa mem -t128 $transcriptome <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}-00?_R1.fq.gz) <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}-00?_R2.fq.gz) 2>> $outdir/bwa-mem-${t}.log | samtools sort -@ 128 -O BAM -o $outdir/Q903.annotation.${t}.${num}.sorted.bam
				fi
				count=$(samtools view -@ 128 -f4 -c $outdir/Q903.annotation.${t}.${num}.sorted.bam)
				echo "Number of unmapped Q903 $t reads: $count/$total" 1>&2
				percentage=$(echo "scale=4;$count / $total * 100" | bc)
				echo "Percentage of unmapped Q903 $t reads: $percentage" 1>&2

				echo -e "Annotation\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv
				if [[ ! -e "${transcriptome_dir}/${t}.fa.sa" ]]
				then
					echo "Indexing assembled Q903 $t transcriptome..." 1>&2
					bwa index $transcriptome_dir/${t}.fa 2>> $outdir/bwa-index-${t}.log
				fi
				if [[ ! -e "$outdir/Q903.assembly.${t}.${num}.sorted.bam" ]]
				then

					echo "Aligning Q903 $t reads to assembled transcriptome..." 1>&2
					bwa mem -t128 $transcriptome_dir/${t}.fa <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}-00?_R1.fq.gz) <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}-00?_R2.fq.gz) 2>> $outdir/bwa-mem-${t}.log | samtools sort -@ 128 -O BAM -o $outdir/Q903.assembly.${t}.${num}.sorted.bam
				fi
				count=$(samtools view -@ 128 -f4 -c $outdir/Q903.assembly.${t}.${num}.sorted.bam)
				echo "Number of unmapped Q903 $t reads: $count/$total" 1>&2
				percentage=$(echo "scale=4;$count / $total * 100" | bc)
				echo "Percentage of unmapped Q903 $t reads: $percentage" 1>&2
				echo -e "Assembly\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv

			done
		else
			# Combine them
			echo "Counting the total number of reads..." 1>&2
			total=$(countReads.sh /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-?-00?_R?.fq.gz | tail -n1 | awk '{print $2}')
			if [[ ! -e "$outdir/Q903.annotation.${t}.sorted.bam" ]]
			then
				echo "Aligning Q903 $t reads to annotated transcriptome..." 1>&2
				bwa mem -t128 $transcriptome <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-?-00?_R1.fq.gz) <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-?-00?_R2.fq.gz) 2>> $outdir/bwa-mem-${t}.log | samtools sort -@ 128 -O BAM -o $outdir/Q903.annotation.${t}.sorted.bam
			fi

			count=$(samtools view -@ 128 -f4 -c $outdir/Q903.annotation.${t}.sorted.bam)
			echo "Number of unmapped Q903 $t reads: $count/$total" 1>&2
			percentage=$(echo "scale=4;$count / $total * 100" | bc)
			echo "Percentage of unmapped Q903 $t reads: $percentage" 1>&2
			echo -e "Annotation\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv
			if [[ -e "${transcriptome_dir}/${t}.fa.sa" ]]
			then
				echo "Indexing assembled Q903 $t transcriptome..." 1>&2
				bwa index $transcriptome_dir/${t}.fa 2>> $outdir/bwa-index-${t}.log
			fi

			if [[ ! -e "$outdir/Q903.assembly.${t}.sorted.bam" ]]
			then
				echo "Aligning Q903 $t reads to assembled transcriptome..." 1>&2
				bwa mem -t128 $transcriptome_dir/${t}.fa <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-?-00?_R1.fq.gz) <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-?-00?_R2.fq.gz) 2>> $outdir/bwa-mem-${t}.log | samtools sort -@ 128 -O BAM -o $outdir/Q903.assembly.${t}.sorted.bam
			fi

			count=$(samtools view -@ 128 -f4 -c $outdir/Q903.assembly.${t}.sorted.bam)
			echo "Number of unmapped Q903 $t reads: $count/$total" 1>&2
			percentage=$(echo "scale=4;$count / $total * 100" | bc)
			echo "Percentage of unmapped Q903 $t reads: $percentage" 1>&2
			echo -e "Assembly\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv
		fi
	done

	tissues="control gallery wound"
	for t in $tissues
	do
		if [[ "$combined" = false ]]
		then
			for num in $(seq 4)
			do

				echo "Counting the total number of reads..." 1>&2
				total=$(countReads.sh /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}_R?.fq.gz | tail -n1 | awk '{print $2}')
				if [[ ! -e "$outdir/Q903.annotation.${t}.${num}.sorted.bam" ]]
				then
					echo "Aligning Q903 $t reads to annotated transcriptome..." 1>&2
					bwa mem -t128 $transcriptome <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}_R1.fq.gz) <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}_R2.fq.gz) 2>> $outdir/bwa-mem-${t}.log | samtools sort -@ 128 -O BAM -o $outdir/Q903.annotation.${t}.${num}.sorted.bam
				fi
				count=$(samtools view -f4 -c -@ 128 $outdir/Q903.annotation.${t}.${num}.sorted.bam)
				echo "Number of unmapped Q903 $t reads: $count/$total" 1>&2
				percentage=$(echo "scale=4;$count / $total * 100" | bc)
				echo "Percentage of unmapped Q903 $t reads: $percentage" 1>&2

				echo -e "Annotation\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv
				if [[ ! -e "$transcriptome_dir/${t}.fa.sa" ]]
				then
					echo "Indexing assembled Q903 $t transcriptome..." 1>&2
					bwa index $transcriptome_dir/${t}.fa 2>> $outdir/bwa-index-${t}.log
				fi
				if [[ ! -e "$outdir/Q903.assembly.${t}.${num}.sorted.bam" ]]
				then
					echo "Aligning Q903 $t reads to assembled transcriptome..." 1>&2
					bwa mem -t128 $transcriptome_dir/${t}.fa <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}_R1.fq.gz) <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}_R2.fq.gz) 2>> $outdir/bwa-mem-${t}.log | samtools sort -O BAM -@ 128 -o $outdir/Q903.assembly.${t}.${num}.sorted.bam
				fi
				count=$( samtools view -f4 -c -@ 128 $outdir/Q903.assembly.${t}.${num}.sorted.bam)
				echo "Number of unmapped Q903 $t reads: $count/$total" 1>&2

				percentage=$(echo "scale=4;$count / $total * 100" | bc)
				echo "Percentage of unmapped Q903 $t reads: $percentage" 1>&2

				echo -e "Assembly\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv
			done
		else
			echo "Counting the total number of reads..." 1>&2
			total=$(countReads.sh /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-?_R?.fq.gz | tail -n1 | awk '{print $2}')
			if [[ ! -e "$outdir/Q903.annotation.${t}.sorted.bam" ]]
			then
				echo "Aligning Q903 $t reads to annotated transcriptome..." 1>&2
				bwa mem -t128 $transcriptome <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-?_R1.fq.gz) <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-?_R2.fq.gz) 2>> $outdir/bwa-mem-${t}.log | samtools sort -@ 128 -O BAM -o $outdir/Q903.annotation.${t}.sorted.bam
			fi
			count=$(samtools view -f4 -c -@ 128 $outdir/Q903.annotation.${t}.sorted.bam)
			echo "Number of unmapped Q903 $t reads: $count/$total" 1>&2
			percentage=$(echo "scale=4;$count / $total * 100" | bc)
			echo "Percentage of unmapped Q903 $t reads: $percentage" 1>&2
			echo -e "Annotation\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv

			if [[ ! -e "$transcriptome_dir/${t}.fa.sa" ]]
			then
				echo "Indexing assembled Q903 $t transcriptome..." 1>&2
				bwa index $transcriptome_dir/${t}.fa 2>> $outdir/bwa-index-${t}.log
			fi

			if [[ ! -e "$outdir/Q903.assembly.${t}.sorted.bam" ]]
			then
				echo "Aligning Q903 $t reads to assembled transcriptome..." 1>&2
				bwa mem -t128 $transcriptome_dir/${t}.fa <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-?_R1.fq.gz) <(pigz -d -c /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-?_R2.fq.gz) 2>> $outdir/bwa-mem-${t}.log | samtools sort -O BAM -@ 128 -o $outdir/Q903.assembly.${t}.sorted.bam
			fi
			count=$( samtools view -f4 -c -@ 128 $outdir/Q903.assembly.${t}.sorted.bam)
			echo "Number of unmapped Q903 $t reads: $count/$total" 1>&2
			percentage=$(echo "scale=4;$count / $total * 100" | bc)
			echo "Percentage of unmapped Q903 $t reads: $percentage" 1>&2
			echo -e "Assembly\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv

		fi
	done
fi
