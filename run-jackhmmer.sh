#!/bin/bash
AMP=$1
lit=$2
database=$3
N=$4
begin=$5
end=$6
scaffolds=$7
transcripts=$8
gff=$9

if [ "$#" -ne 9 ]
then
	echo "USAGE: $(basename $0) <literature AMPs> <NCBI defensins> <maker predicted proteins> <# of iterations> <sweep start> <sweep end> <spruce scaffolds> <spruce transcripts> <spruce GFF>"
	echo "To run jackhmmer without a sweep, set:\n<sweep start> = <sweep end>."
	exit 1
fi

function run_jackhmmer() {
	AMP=$1
	lit=$2
	database=$3
	N=$4
	threshold=$5

	echo "Bit score $threshold detected..."
	#in the case that <sweep end> = <sweep start>, meaning no sweep is needed, just a straight-forward jackhmmer with a specified threshold
	outfile="jackhmmer_bs${threshold}_N${N}.out"
	while true
	do
		echo "Running jackhmmer with a threshold of $threshold for $N iterations..."
		jackhmmer --noali -T $threshold -N $N -o $outfile $AMP $database
		converged=$(grep -c 'CONVERGED' $outfile)
		total=$(grep -c 'Query:' $outfile)
		#If not converged, increase iterations and delete the file
		if [ "$converged" -ne "$total" ]
		then
			rm $outfile
			N=$((N+5))
			echo "At bit score threshold $threshold, not all queries converged. Increasing N to $N."
		else
			#Once converged, stop increasing iterations and break
			break
		fi
	done
}

#Guide Blast - directly blast known/literature defensins against the database to see which proteins we cannot lose - threshold 99% identity

if [ "$end" -ne "$begin" ]
then
	difference=$((end-begin))
	if [ "$difference" -gt 100 ]
	then
		step=100
	elif [ "$difference" -gt 10 ]
	then
		step=10
	else
		step=1
	fi
	#if difference is not divisble by the assigned step, ask user to change it.
	while [[ $((difference%step)) -ne 0 && "$difference" -ne 0 ]]
	do
		echo "The difference between <sweep start> and <sweep end> must be a multiple of $step."
		echo "Please enter a new <sweep start> and <sweep end> on the line below (separated by a space). To run jackhmmer at a specific threshold (no sweep), please enter the desired threshold as a single value."
		read newbegin newend
		if [ ! -z "$newbegin" ] && [ -z "$newend" ]
		then
			end=0
			begin=$newbegin
			break
		elif [[ ! -z "$newbegin" && ! -z "$newend" ]]
		then
			difference=$((newend-newbegin))
			if [ "$difference" -gt 100 ]
			then
				step=100
			elif [ "$difference" -gt 10 ]
			then
				step=10
			else
				step=1
			fi
			begin=$newbegin
			end=$newend
		else
			continue
		fi
	done

	if [ "$end" -ne 0 ] && [ "$begin" -ne "$end" ]
	then
		#Testing this script for bugs -- do not build BLASTDB and BLASTP unnecessarily
		if [ ! -e "guide-proteins.txt" ]
		then
			echo "Making BLAST database..."
			makeblastdb -dbtype prot -in $database -out blastpdb
			echo -e "\nBLASTing..."
			blastp -db blastpdb -query $lit -out guide-blast.blastp -outfmt '6 std qcovs' -num_threads 48
			#filter blastp results for 99% identity sequences
			awk '{if($3>99) print $2}' guide-blast.blastp | sort -u > guide-proteins.txt
		fi

		#see what threshold we lose these proteins - do a jackhmmer sweep where grep for guide-proteins at each threshold
		echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
		jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
		sweep=$?
		#if sweep > 0, then script executed with no problems
		while [ "$sweep" -gt 0 ]
		do
			#if sweep = 1, jackhmmer did not converge
			if [ "$sweep" -eq 1 ]
			then
				N=$((N+5))
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
				jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
				sweep=$?
			#if sweep = 2, the sweep start is too low
			elif [ "$sweep" -eq 2 ]
			then
				end=$begin
				begin=$((begin-difference))
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
				jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
				sweep=$?
			#if sweep > 2 (aka a threshold), then reduce the interval and find the threshold
			elif [ "$sweep" -gt 2 ]
			then
				end=$sweep
				begin=$((end-step))
				if [ "$step" -eq 100 ]
				then
					step=10
				elif [ "$step" -eq 10 ]
				then
					step=1
				else
					break
				fi
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
				jackhmmer-sweep.sh $AMP $lit $database $N $((begin+step)) $((end-step)) $step
				sweep=$?
			fi		
		done
		threshold=$((sweep-1))
		echo "$threshold is the optimal jackhmmer threshold!"
		outfile="jackhmmer_bs${threshold}_N${N}.out"
		#Delete all jackhmmer output files that are unnecessary
		for file in jackhmmer_bs*_N*.out
		do
			if [ "$file" == "$outfile" ]
			then
				continue
			fi
			rm $file
		done
	else
		run_jackhmmer $AMP $lit $database $N $begin
		outfile="jackhmmer_bs${begin}_N${N}.out"
	fi
else
	run_jackhmmer $AMP $lit $database $N $begin
	outfile="jackhmmer_bs${begin}_N${N}.out"
fi

# Filter jackhmmer results for all hits (start with >>), and remove duplicates, then get their sequences using their names
echo "Running seqtk..."
seqtk subseq $database <(awk '/^>>/ {print $2}' $outfile | sort -u) > jackhmmer-hits.faa

#For debugging, avoid unnecessary computations
if [ ! -e "jackhmmer-blast-hits.faa" ]
then
	echo "Making BLAST database..."
	makeblastdb -dbtype prot -in jackhmmer-hits.faa -out jackhmmer-db

	echo -e "\nBLASTing..."
	blastp -db jackhmmer-db -query $lit -out jackhmmer.blastp -outfmt '6 std qcovs' -num_threads 48
	echo "Running seqtk..."
	seqtk subseq $database <(awk '{if ($3>90) print $2}' jackhmmer.blastp | sort -u) > jackhmmer-blast-hits.faa
fi

#GET scaffolds, gffs and transcripts and then run in WS777111-proteins/test
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

for i in $(awk '{if ($3>90) print $2}' jackhmmer.blastp | sort -u)
do
	#Get scaffold name
	temp=${i#*-}
	scaffold=${temp%-*-gene-*-mRNA-?}
	if [ ! -e "scaffolds/${scaffold}.scaffold.fa" ]
	then
		seqtk subseq $scaffolds <(echo $scaffold) > scaffolds/${scaffold}.scaffold.fa
		length=$(tail -n 1 "scaffolds/${scaffold}.scaffold.fa" | head -c -1 | wc -m)
		echo -e "##gff-version 3\n##sequence-region $scaffold 1 $length" > gffs/${scaffold}.gff
		cd igv
		#Create IGV softlinks
		ln -sf ../scaffolds/${scaffold}.scaffold.fa
		ln -sf ../gffs/${scaffold}.gff
		cd ..
	fi
	seqtk subseq $transcripts <(echo $i) >> transcripts/${scaffold}.transcripts.fa
	grep $i $gff >> gffs/${scaffold}.gff
done

cd scaffolds
cat *.scaffold.fa > all.scaffolds.fa
cd ../transcripts
cat *.transcripts.fa > all.transcripts.fa
cd ../gffs
cat *.gff > all.gff
cd ..
ln -sf scaffolds/all.scaffolds.fa
ln -sf transcripts/all.transcripts.fa
ln -sf gffs/all.gff
echo "DONE!"
