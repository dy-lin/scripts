#!/bin/bash
# Using a non-zero exit status as a 'return value', so commenting out pipefail
# set -eu -o pipefail
PROGRAM=$(basename $0)
if [[ "$#" -ne 6 && "$#" -ne 5 ]]
then
	if [[ "$#" -gt 0 ]]
	then
		echo "ERROR: Incorrect number of arguments!" 1>&2
	fi
	echo "USAGE: $PROGRAM <literature AMPs> <class of AMPs> <protein database> <# of iterations> <sweep start> <sweep end>" 1>&2
	echo -e "\tTo run jackhmmer without a sweep, just input one value instead of two." 1>&2
	echo "DESCRIPTION: Runs jackhmmer of literature AMPs against a protein database in order to find a specific class of AMP using homology-based search. Sweep is run as to find the best threshold to run jackhmmer." 1>&2
	exit 1
fi

AMP=$1
lit=$2
database=$3
N=$4
begin=$5

if [[ ! -e $AMP ]]
then
	echo "ERROR: $(basename $AMP) does not exist."
	exit 1
fi

if [[ ! -e $lit ]]
then
	echo "ERROR: $(basename $lit) does not exist."
	exit 1
fi

if [[ ! -e $database ]]
then
	echo "ERROR: $(basename $lit) does not exist."
	exit 1
fi

if [[ "$N" -le 0 ]]
then
	echo "ERROR: Invalid number of iterations: $N"
	exit 1
fi

if [[ ! -z $6 ]]
then
	if [[ "$begin" -le 0 ]]
	then
		echo "ERROR: Invalid starting threshold: $begin"
		exit 1
	fi
	end=$6
	if [[ "$end" -le 0 ]]
	then
		echo "ERROR: Invalid ending threshold: $end"
		exit 1
	fi
else
	if [[ "$begin" -le 0 ]]
	then
		echo "ERROR: Invalid threshold: $begin"
		exit 1
	fi
	end=$begin
fi

function run_jackhmmer() {
	AMP=$1
	lit=$2
	database=$3
	N=$4
	threshold=$5

	# In the case that <sweep end> = <sweep start>, meaning no sweep is needed, just a straight-forward jackhmmer with a specified threshold
	outfile="jackhmmer_bs${threshold}_N${N}.out"
	while true
	do
		echo "Running jackhmmer with a threshold of $threshold for $N iterations..."
		echo "COMMAND: jackhmmer --noali -T $threshold -N $N -o $outfile $AMP $database"
		jackhmmer -h | head -n 5
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
			return
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

	if [[ "$end" -gt 0 || "$begin" -ne "$end" ]]
	then
		if [[ ! -e guide-proteins.txt ]]
		then
			echo "Making BLAST database..."
			makeblastdb -dbtype prot -in $database -out guide-blast
			echo -e "\nBLASTing..."
			blastp -db guide-blast -query $lit -out guide-blast.blastp -outfmt '6 std qcovs' -num_threads 48
			# Filter blastp results for 99% identity sequences
			cutoff=99
			while true
			do
				awk -v var=$cutoff '{if($3>var) print $2}' guide-blast.blastp | sort -u > guide-proteins.txt
				size=$(ls -l guide-proteins.txt | awk '{print $5}')
				if [[ "$size" -eq 0 ]]
				then
					cutoff=$((cutoff-1))
				else
					if [[ "$cutoff" -le 50 ]]
					then
						echo "There are no proteins that align significantly that can be used as your guide proteins. All proteins align with percent identity 50% or lower."
						rm guide-proteins.txt
						echo "Status: Failure."
						exit 1
					else
						echo "Guide proteins align with percent identity ${cutoff}% or higher."
						break
					fi
				fi
			done
		fi
		seqtk subseq $database guide-proteins.txt > guide-proteins.faa
		# See what threshold we lose these proteins - do a jackhmmer sweep where grep for guide-proteins at each threshold
		echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
		jackhmmer -h | head -n 5
		jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step 1 
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
			# If sweep = 2, the sweep start is too high
			elif [ "$sweep" -eq 2 ]
			then
				end=$begin
				if [ "$begin" -ne 0 ]
				then
					begin=$((begin-difference))
				fi
				if [ "$begin" -le 0 ]
				then
					begin=0
				fi
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..." 1>&2
				jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
				sweep=$?

				if [[ "$begin" -eq 0 && "$step" -eq 1 ]]
				then
					echo "<sweep start> cannot be lowered anymore. Your guide proteins are nowhere to be found." 1>&2
					rm guide-proteins.txt
					rm jackhmmer_bs*_N*.out
					zero=true
					break
				fi
			# If sweep = 3, then the whole sweep finished with no problems -- need to increase the whole interval
			elif [[ "$sweep" -eq 3 ]]
			then
				interval=$((end-begin))
				begin=$end
				end=$((begin + interval))
				if [ "$step" -eq 1 ]
				then
					step=10
				fi
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..." 1>&2
				jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
				sweep=$?
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
				elif [ "$step" -eq 1 ]
				then
					break
				fi
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..."
				jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step
				sweep=$?
			fi
		done
		if [[ "$zero" = true ]]
		then
			run_jackhmmer $AMP $lit $database $N 0
			outfile="jackhmmer_bs0_N${N}.out"
		else
			if [[ "$sweep" -eq 0 ]]
			then
				run_jackhmmer $AMP $lit $database $N 0
				run_jackhmmer $AMP $lit $database $N 1
				run_jackhmmer $AMP $lit $database $N 2
				run_jackhmmer $AMP $lit $database $N 3
				for num in $(seq 0 3)
				do
					outfile="jackhmmer_bs${num}_N${N}.out"
					for prot in $(cat guide-proteins.txt)
					do
						count=$(grep -c $prot $outfile)
						if [[ "$count" -eq 0 ]]
						then
							if [[ "$num" -le 1 ]]
							then
								echo "There is no need to run jackhmmer with a threshold. The lowest threshold possible still excludes all your guide proteins."
								mv jackhmmer_bs0_N${N}.out jackhmmer_N${N}.out
								rm jackhmmer_bs*_N*.out
								mv jackhmmer_N${N}.out jackhmmer_bs0_N${N}.out
								outfile="jackhmmer_bs0_N${N}.out"
							else
								echo "Bit score ${num} loses guide protein ${prot} so the ideal threshold is $((num-1))!"
								mv jackhmmer_bs$((num-1))_N${N}.out jackhmmer_N${N}.out
								rm jackhmmer_bs*_N*.out
								mv jackhmmer_N${N}.out jackhmmer_bs$((num-1))_N${N}.out
								outfile="jackhmmer_bs$((num-1))_N${N}.out"
							fi
							finished=true
							break
						fi
					done
					if [[ "$finished" = true ]]
					then
						break
					fi
				done
			else
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
			fi
		fi
	else
		outfile="jackhmmer_bs${begin}_N${N}.out"
		echo "User input bit score $begin threshold detected!"
		run_jackhmmer $AMP $lit $database $N $begin
	fi
else
	outfile="jackhmmer_bs${begin}_N${N}.out"
	echo "Bit score $begin threshold detected!"
	run_jackhmmer $AMP $lit $database $N $begin
fi

# Filter jackhmmer results for all hits (start with >>), and remove duplicates, then get their sequences using their names
echo "Running seqtk..."
seqtk subseq $database <(awk '/^>>/ {print $2}' $outfile | sort -u) > jackhmmer-hits.faa

if [ -e jackhmmer-hits.faa ]
then
	file_size=$(ls -l jackhmmer-hits.faa | awk '{print $5}')
	if [ "$file_size" -eq 0 ]
	then
		echo "There were no jackhmmer hits."
		echo "Status: Failure."
		rm jackhmmer-hits.faa
		exit 1
	fi
fi
jackhmmer-blast.sh $lit $database $outfile
