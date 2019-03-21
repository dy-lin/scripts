#!/bin/bash
set -eu -o pipefail
PROGRAM=$(basename $0)
gethelp=false
verbose=false

while getopts :hv opt
do
	case $opt in
		h) gethelp=true;;
		v) export verbose=true;;
		\?) echo "ERROR: $PROGRAM: Invalid option: $opt" 1>&2 ; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ "$#" -ne 6 && "$#" -ne 5 ]] || [[ "$gethelp" = true ]]
then
	if [[ "$#" -gt 0 ]]
	then
		echo "ERROR: Incorrect number of arguments!" 1>&2
	fi
	echo "USAGE: $PROGRAM <literature AMPs> <class of AMPs> <protein database> <# of iterations> <sweep start> <sweep end>" 1>&2
	echo -e "\tTo run jackhmmer without a sweep, just input one value instead of two." 1>&2
	echo "DESCRIPTION: Runs jackhmmer of literature AMPs against a protein database in order to find a specific class of AMP using homology-based search. Sweep is run as to find the best threshold to run jackhmmer." 1>&2
	echo -e "OPTIONS:\n\t-h\tShow help menu\n\t-v\tverbose logging" 1>&2
	exit 1
fi
if [[ "$verbose" = true ]]
then
	echo "COMMAND: $PROGRAM $*" 1>&2
fi

AMP=$1
lit=$2
database=$3
N=$4
begin=$5

if [[ ! -e $AMP ]]
then
	echo "ERROR: $(basename $AMP) does not exist." 1>&2
	exit 1
fi

if [[ ! -e $lit ]]
then
	echo "ERROR: $(basename $lit) does not exist." 1>&2
	exit 1
fi

if [[ ! -e $database ]]
then
	echo "ERROR: $(basename $lit) does not exist." 1>&2
	exit 1
fi

if [[ "$N" -le 0 ]]
then
	echo "ERROR: Invalid number of iterations: $N" 1>&2
	exit 1
fi

if [[ ! -z $6 ]]
then
	if [[ "$begin" -le 0 ]]
	then
		echo "ERROR: Invalid starting threshold: $begin" 1>&2
		exit 1
	fi
	end=$6
	if [[ "$end" -le 0 ]]
	then
		echo "ERROR: Invalid ending threshold: $end" 1>&2
		exit 1
	fi
else
	if [[ "$begin" -le 0 ]]
	then
		echo "ERROR: Invalid threshold: $begin" 1>&2
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
		echo "Running jackhmmer with a threshold of $threshold for $N iterations..." 1>&2
		if [[ "$verbose" = true ]]
		then
			echo "COMMAND: jackhmmer --noali -T $threshold -N $N -o $outfile $AMP $database" 1>&2
			jackhmmer -h | head -n 5
		fi
		if [[ ! -e "$outfile" ]]
		then
			jackhmmer --noali -T $threshold -N $N -o $outfile $AMP $database
		fi
		converged=$(grep -c 'CONVERGED' $outfile)
		total=$(grep -c 'Query:' $outfile)
		# If not converged, increase iterations and delete the file
		if [ "$converged" -ne "$total" ]
		then
			if [[ "$threshold" -eq 0 && "$N" -ge 100 ]]
			then
				return
			fi
			rm $outfile
			N=$((N+5))
			echo "At bit score threshold $threshold, not all queries converged. Increasing N to $N." 1>&2
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
			echo "Making BLAST database..." 1>&2
			if [[ "$verbose" = true ]]
			then
				echo "COMMAND: makeblastdb -dbtype prot -in $database -out guide-blast" 1>&2
			fi
			makeblastdb -dbtype prot -in $database -out guide-blast
			echo -e "\nBLASTing..." 1>&2
			if [[ "$verbose" = true ]]
			then
				echo "COMMAND: blastp -db guide-blast -query $lit -out guide-blast.blastp -outfmt '6 std qcovs' -num_threads 48" 1>&2
			fi
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
						echo "There are no proteins that align significantly that can be used as your guide proteins. All proteins align with percent identity 50% or lower." 1>&2
						echo "Using all guide blast results as guide proteins." 1>&2
					else
						echo "Guide proteins align with percent identity ${cutoff}% or higher." 1>&2
					fi
						break
				fi
			done
		fi
		seqtk subseq $database guide-proteins.txt > guide-proteins.faa
		# See what threshold we lose these proteins - do a jackhmmer sweep where grep for guide-proteins at each threshold
		echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..." 1>&2
		if [[ "$verbose" = true ]]
		then
			jackhmmer -h | head -n 5
			echo "COMMAND: jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step" 1>&2
		fi
		sweep=$(jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step)
		# If sweep > 0, then script executed with no problems
		while true
		do
			# If sweep = 1, jackhmmer did not converge
			if [ "$sweep" == "nc" ]
			then
				N=$((N+5))
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..." 1>&2
				if [[ "$verbose" = true ]]
				then
					echo "COMMAND: jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step" 1>&2
				fi
				sweep=$(jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step)
			# If sweep = 2, the sweep start is too high
			elif [ "$sweep" == "high" ]
			then
				end=$begin
				if [[ "$((begin-difference))" -le 0 && "$begin" -ne 0 ]]
				then
					begin=$((difference/2))
					difference=$((end-begin))
				elif [ "$begin" -le 1 ]
				then
					begin=0
				else
					begin=$((end-difference))
				fi
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..." 1>&2
				if [[ "$verbose" = true ]]
				then
					echo "COMMAND: jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step" 1>&2
				fi
				sweep=$(jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step)
				if [[ "$begin" -eq 0 && "$step" -eq 1 ]]
				then
					echo "<sweep start> cannot be lowered anymore. Your guide proteins are nowhere to be found." 1>&2
					rm guide-proteins.txt
					rm jackhmmer_bs*_N*.out
					zero=true
					break
				fi
			elif [[ "$sweep" == "np" ]]
			then
				interval=$((end-begin))
				begin=$end
				end=$((begin + interval))
				if [ "$step" -eq 1 ]
				then
					step=10
				fi
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..." 1>&2
				if [[ "$verbose" = true ]]
				then
					echo "COMMAND: jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step" 1>&2
				fi
				sweep=$(jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step)
			else
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
				echo "Conducting jackhmmer sweep from $begin to $end in $step-step intervals for $N iterations..." 1>&2
				if [[ "$verbose" = true ]]
				then
					echo "COMMAND: jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step" 1>&2
				fi
				sweep=$(jackhmmer-sweep.sh $AMP $lit $database $N $begin $end $step)
			fi
		done
		if [[ "$zero" = true ]]
		then
			run_jackhmmer $AMP $lit $database $N 0
			outfile="jackhmmer_bs0_N${N}.out"
		fi
		threshold=$((sweep-1))
		echo -e "\nBit score threshold $threshold is the optimal threshold to use when running jackhmmer!\n" 1>&2
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
		outfile="jackhmmer_bs${begin}_N${N}.out"
		echo "User input bit score $begin threshold detected!" 1>&2
		run_jackhmmer $AMP $lit $database $N $begin
	fi
else
	outfile="jackhmmer_bs${begin}_N${N}.out"
	echo "Bit score $begin threshold detected!" 1>&2
	run_jackhmmer $AMP $lit $database $N $begin
fi

# Filter jackhmmer results for all hits (start with >>), and remove duplicates, then get their sequences using their names
echo "Running seqtk..." 1>&2
seqtk subseq $database <(awk '/^>>/ {print $2}' $outfile | sort -u) > jackhmmer-hits.faa

if [ -e jackhmmer-hits.faa ]
then
	file_size=$(ls -l jackhmmer-hits.faa | awk '{print $5}')
	if [ "$file_size" -eq 0 ]
	then
		echo "There were no jackhmmer hits." 1>&2
		echo "Status: Failure." 1>&2
		rm jackhmmer-hits.faa
		exit 1
	fi
fi
jackhmmer-blast.sh $lit $database $outfile
