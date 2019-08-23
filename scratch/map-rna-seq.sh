#!/bin/bash


# PG29 first
tissues="bark embryo flush_bud mature_needle megagametophyte seed_germination xylem young_buds"
outdir="/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/alignment"
mkdir -p $outdir

echo -e "Transcriptome\tReads\t\# of Unmapped Reads\tPercentage (%)" > $outdir/PG29-map-rna.tsv

transcriptome="/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/PG29.maker.transcripts.edited.fa"
transcriptome_dir="/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/transcriptomes"

bwa index $transcriptome 2>> $outdir/bwa-index.log
for t in $tissues
do
#	bwa index $transcriptome 2>> $outdir/bwa-index-${t}.log
	count=$(bwa mem -t128 $transcriptome /projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/PG29_reads/${t}.R?.fq.gz 2>> $outdir/bwa-mem-${t}.log | samtools sort | samtools view -f4 -c)
	total=$(countReads.sh /projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/PG29_reads/${t}.R?.fq.gz | tail -n1 | awk '{print $2}')
	percentage=$(echo "scale=4;$count / $total * 100" | bc)

	echo -e "Annotation\t$t\t$count\t$percentage" >> $outdir/PG29-map-rna.tsv

	bwa index $transcriptome_dir/${t}.fa 2>> $outdir/bwa-index-${t}.log
	count=$(bwa mem -t128 $transcriptome_dir/${t}.fa /projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/PG29_reads/${t}.R?.fq.gz 2>> $outdir/bwa-mem-${t}.log | samtools sort | samtools view -f4 -c)
	total=$(countReads.sh /projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto/PG29_reads/${t}.R?.fq.gz | tail -n1 | awk '{print $2}')
	percentage=$(echo "scale=4;$count / $total * 100" | bc)

	echo -e "Assembly\t$t\t$count\t$percentage" >> $outdir/PG29-map-rna.tsv

done

# Q903

tissues="Dev_SC Cort_Par"
outdir="/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/alignment"
mkdir -p $outdir


echo -e "Transcriptome\tReads\t\# of Unmapped Reads\tPercentage (%)" > $outdir/Q903-map-rna.tsv

transcriptome="/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903.maker.transcripts.edited.fa"
transcriptome_dir="/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/transcriptomes"


bwa index $transcriptome 2>> $outdir/bwa-index.log

for t in $tissues
do
	for num in $(seq 4 2 8)
	do
#		bwa index $transcriptome 2>> $outdir/bwa-index-${t}.log
		count=$(bwa mem -t128 $transciptome /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}-00?_R?.fq.gz 2>> $outdir/bwa-mem-${t}.log | samtools sort | samtools view -f4 -c)
		total=$(countReads.sh /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}-00?_R?.fq.gz | tail -n1 | awk '{print $2}')
		percentage=$(echo "scale=4;$count / $total * 100" | bc)

		echo -e "Annotation\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv

		bwa index $transcriptome_dir/${t}.fa 2>> $outdir/bwa-index-${t}.log
		count=$(bwa mem -t128 $transciptome_dir/${t}.fa /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}-00?_R?.fq.gz 2>> $outdir/bwa-mem-${t}.log | samtools sort | samtools view -f4 -c)
		total=$(countReads.sh /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}-00?_R?.fq.gz | tail -n1 | awk '{print $2}')
		percentage=$(echo "scale=4;$count / $total * 100" | bc)
		echo -e "Assembly\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv

	done
done

tissues="control gallery wound"
for t in $tissues
do
	for num in $(seq 4)
	do
#		bwa index $transcriptome 2>> $outdir/bwa-index-${t}.log
		count=$(bwa mem -t128 $transcriptome /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}_R?.fq.gz 2>> $outdir/bwa-mem-${t}.log | samtools sort | samtools view -f4 -c)
		total=$(countReads.sh /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}_R?.fq.gz | tail -n1 | awk '{print $2}')
		percentage=$(echo "scale=4;$count / $total * 100" | bc)

		echo -e "Annotation\t$t\t$count\t$percentage" >> Q903-map-rna.tsv

		bwa index $transcriptome_dir/${t}.fa 2>> $outdir/bwa-index-${t}.log
		count=$(bwa mem -t128 $transcriptome_dir/${t}.fa /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}_R?.fq.gz 2>> $outdir/bwa-mem-${t}.log | samtools sort | samtools view -f4 -c)
		total=$(countReads.sh /projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/Q903_reads/${t}-${num}_R?.fq.gz | tail -n1 | awk '{print $2}')
		percentage=$(echo "scale=4;$count / $total * 100" | bc)

		echo -e "Assembly\t$t\t$count\t$percentage" >> $outdir/Q903-map-rna.tsv
	done
done
