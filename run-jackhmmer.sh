#!/bin/bash
AMP=$1
lit=$2
database=$3
N=$4
begin=$5
end=$6
step=10
scaffolds=$7
transcripts=$8
gff=$9

if [ "$#" -ne 9 ]
then
	echo "USAGE: $(basename $0) <literature AMPs> <NCBI defensins> <maker predicted proteins> <# of iterations> <sweep start> <sweep end> <spruce scaffolds> <spruce transcripts> <spruce GFF>"
	echo "To run jackhmmer without a sweep, set:\n<sweep start> = desired threshold\n<sweep end> = 0."
	exit 1
fi

#Guide Blast - directly blast known/literature defensins against the database to see which proteins we cannot lose - threshold 99% identity

if [ "$end" != "0" ]
then
	echo "Making BLAST database..."
	makeblastdb -dbtype prot -in $database -out blastpdb
	echo -e "\nBLASTing..."
	blastp -db blastpdb -query $lit -out guide-blast.blastp -outfmt '6 std qcovs' -num_threads 48
	#filter blastp results for 99% identity sequences
	awk '{if($3>99) print $2}' guide-blast.blastp | sort -u > guide-proteins.txt

	#see what threshold we lose these proteins
	#check for convergence
	echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
	jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
	sweep=$?
	while [ "$sweep" -gt 0 ]
	do

		if [ "$sweep" -eq 1 ]
		then
			N=$((N+5))
			echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
			jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
			sweep=$?
		elif [ "$sweep" -gt 1 ]
		then
			echo "Conducting jackhmmer sweep from $((sweep-step)) to $((sweep)) in 1-step intervals for $N iterations..."
			jackhmmer-sweep.sh $AMP $lit $database $N $((sweep-step+1)) $((sweep-1)) 1
			sweep=$?
			step=1
			if [ "$step" -eq 1 ]
			then
				break
			fi
		fi		
	done
	echo "$((sweep-1)) is the optimal jackhmmer threshold!"
	outfile="jackhmmer_bs$((sweep-1))_N${N}.out"
	for file in jackhmmer_bs*_N${N}.out
	do
		if [ "$file" == "$outfile" ]
		then
			continue
		fi
		rm $file
	done
else
	outfile="jackhmmer_bs${begin}_N${N}.out"
	while true
	do
		echo "Running jackhmmer with a threshold of $begin for $N iterations..."
		jackhmmer --noali -T $begin -N $N -o $outfile $AMP $database
		converged=$(grep -c 'CONVERGED' $outfile)
		total=$(grep -c 'Query:' $outfile)
		if [ "$converged" -ne "$total" ]
		then
			rm $outfile
			N=$((N+5))
			echo "Not all queries converged. Increasing N to $N."
		else
			break
		fi
	done
fi

echo "Running seqtk..."
seqtk subseq $database <(awk '/^>>/ {print $2}' $outfile | sort -u) > jackhmmer-hits.faa

echo "Making BLAST database..."
makeblastdb -dbtype prot -in jackhmmer-hits.faa -out jackhmmer-db

echo -e "\nBLASTing..."
blastp -db jackhmmer-db -query $lit -out jackhmmer.blastp -outfmt '6 std qcovs' -num_threads 48
echo "Running seqtk..."
seqtk subseq $database <(awk '{if ($3>90) print $2}' jackhmmer.blastp | sort -u) > jackhmmer-blast-hits.faa

#GET scaffolds, gffs and transcripts and then run in WS777111-proteins/test
echo "Fetching scaffolds, transcripts and GFF files..."
for i in $(awk '{if ($3>90) print $2}' jackhmmer.blastp | sort -u)
do
	temp=${i#*-}
	scaffold=${temp%-*-gene-*-mRNA-?}
	seqtk subseq $scaffolds <(echo $scaffold) > ${scaffold}.scaffold.fa
	seqtk subseq $transcripts <(echo $i) > ${scaffold}.transcripts.fa
	grep $i $gff > ${scaffold}.gff
done
echo "DONE!"
