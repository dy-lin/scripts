#!/bin/bash
# Using a non-zero exit status as a 'return value', so commenting out pipefail
# set -eu -o pipefail
PROGRAM=$(basename $0)
if [[ "$#" -ne 6 && "$#" -ne 5 ]]
then
	if [[ "$#" -ne 6  && "$#" -ne 5 && "$#" -gt 0 ]]
	then
		echo "ERROR: Incorrect number of arguments!" 1>&2
	fi
	echo "USAGE: $PROGRAM <literature AMPs> <NCBI defensins> <protein database> <# of iterations> <sweep start> <sweep end>" 1>&2
	echo -e "\tTo run jackhmmer without a sweep, just input one value instead of two." 1>&2
	echo "DESCRIPTION: Runs jackhmmer of literature AMPs against a protein database in order to find a specific class of AMP using homology-based search. Sweep is run as to find the best threshold to run jackhmmer." 1>&2
	exit 1
fi

AMP=$1
lit=$2
database=$3
N=$4
begin=$5
if [[ "$#" -eq 6 ]]
then
	end=$6
else
	end=$5
fi

function run_jackhmmer() {
	AMP=$1
	lit=$2
	database=$3
	N=$4
	threshold=$5

	echo "Bit score $threshold detected..."
	# In the case that <sweep end> = <sweep start>, meaning no sweep is needed, just a straight-forward jackhmmer with a specified threshold
	outfile="jackhmmer_bs${threshold}_N${N}.out"
	while true
	do
		echo "Running jackhmmer with a threshold of $threshold for $N iterations..."
		jackhmmer --noali -T $threshold -N $N -o $outfile $AMP $database
		converged=$(grep -c 'CONVERGED' $outfile)
		total=$(grep -c 'Query:' $outfile)
		# If not converged, increase iterations and delete the file
		if [ "$converged" -ne "$total" ]
		then
			rm $outfile
			N=$((N+5))
			echo "At bit score threshold $threshold, not all queries converged. Increasing N to $N."
		else
			# Once converged, stop increasing iterations and break
			break
		fi
	done
}

# Guide Blast - directly blast known/literature defensins against the database to see which proteins we cannot lose - threshold 99% identity

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
	# If difference is not divisble by the assigned step, ask user to change it
	while [[ $((difference%step)) -ne 0 && "$difference" -ne 0 ]]
	do
		echo "The difference between <sweep start> and <sweep end> must be a multiple of $step." 1>&2
		echo "Please enter a new <sweep start> and <sweep end> on the line below (separated by a space). To run jackhmmer at a specific threshold (no sweep), please enter the desired threshold as a single value." 1>&2
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
		# Existence test used for debugging purposes
		if [ ! -e "guide-proteins.txt" ]
		then
			echo "Making BLAST database..."
			makeblastdb -dbtype prot -in $database -out guide-blast
			echo -e "\nBLASTing..."
			blastp -db guide-blast -query $lit -out guide-blast.blastp -outfmt '6 std qcovs' -num_threads 48
			# Filter blastp results for 99% identity sequences
			awk '{if($3>99) print $2}' guide-blast.blastp | sort -u > guide-proteins.txt
			seqtk subseq $database guide-proteins.txt > guide-proteins.faa
		fi

		# See what threshold we lose these proteins - do a jackhmmer sweep where grep for guide-proteins at each threshold
		echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
		jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
		sweep=$?
		# If sweep > 0, then script executed with no problems
		while [ "$sweep" -gt 0 ]
		do
			# If sweep = 1, jackhmmer did not converge
			if [ "$sweep" -eq 1 ]
			then
				N=$((N+5))
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
				jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
				sweep=$?
			# If sweep = 2, the sweep start is too low
			elif [ "$sweep" -eq 2 ]
			then
				if [ "$begin" -eq 1 ]
				then
					echo "<sweep start> cannot be lowered anymore. Your guide proteins are nowhere to be found." 1>&2
					exit 1
				fi
				end=$begin
				begin=$((begin-difference))
				if [ "$begin" -eq 0 ]
				then
					begin=1
				fi
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
				jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
				sweep=$?
			# If sweep = 3, then the whole sweep finished with no problems -- need to lower the whole interval
			elif [[ "$sweep" -eq 3 ]]
			then
				begin=$end
				end=$((begin + step))
			# If sweep > 2 (aka a threshold), then reduce the interval and find the threshold

			elif [ "$sweep" -gt 3 ]
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
		echo -e "\nBit score threshold $threshold is the optimal threshold to use when running jackhmmer!\n"
		outfile="jackhmmer_bs${threshold}_N${N}.out"
		# Delete all jackhmmer output files that are unnecessary
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

# Existence Test for debugging purposes
if [ ! -e "jackhmmer-blast-hits.faa" ]
then
	echo "Making BLAST database..."
	makeblastdb -dbtype prot -in jackhmmer-hits.faa -out jackhmmer

	echo -e "\nBLASTing..."
	blastp -db jackhmmer -query $lit -out jackhmmer.blastp -outfmt '6 std qcovs' -num_threads 48
	echo "Running seqtk..."
	seqtk subseq $database <(awk '{if ($3>90) print $2}' jackhmmer.blastp | sort -u) > jackhmmer-blast-hits.faa
fi


