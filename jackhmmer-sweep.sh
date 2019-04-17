#!/bin/bash
# set -eu -o pipefail
AMP=$1
lit=$2
database=$3
N=$4
begin=$5
end=$6
step=$7
if [[ "$begin" -eq 0 ]]
then
	difference=$end
else
	difference=$((end-begin))
fi

for i in `seq $begin $step $end`
do	
	outfile="jackhmmer_bs${i}_N${N}.out"
	# If N-5 exists, it means the jackhmmer for that iteration already converged
	# So there is no need to run jackhmmer at that threshold again at a higher itereration
	# If N-5 does not exist, either jackhmmer has never been run for that threshold yet, or it had not converged previously.
	if [[ "$(ls jackhmmer_bs${i}_N*.out 2> /dev/null | wc -l)" -eq 0 ]]
	then
		if [ ! -e "jackhmmer_bs${i}_N${N}.out" ]
		then
			if [[ "$verbose" = true ]]
			then
				echo -e "\tCOMMAND: jackhmmer --noali -T $i -N $N -o $outfile $AMP $database" 1>&2
			fi
			jackhmmer --noali --notextw -T $i -N $N -o $outfile $AMP $database
		fi
	else
		# Pick the largest converged iteration.
		current=$(ls jackhmmer_bs${i}_N*.out | awk -F "_N" '{print $2}' | sed 's/.out$//' | sort -g -r | head -n1)
		if [[ "$current" -ne "$N" ]]
		then
			echo "At bit score $i, all queries have previously converged in $current iterations. Renaming the file to match current N value." 1>&2
			mv jackhmmer_bs${i}_N${current}.out $outfile
		fi
	fi
	# If unconverged, delete the file -- all files that have not been deleted will have been converged/proteins found
	converged=$(grep -c 'CONVERGED' $outfile)
	total=$(grep -c 'Query:' $outfile)
	if [ "$converged" -ne "$total" ]
	then
		if [[ "$x" -ge 5 ]]
		then
			iter=$((iter*5))
		fi
		echo "At bit score threshold $i, not all queries have converged. Increasing N to $((N+iter))." 1>&2
		rm $outfile
		echo "nc"
		exit 0
	fi
	threshold=$i
#	threshold=$(echo $outfile | awk -F "_" '{print $2}' | awk -F "s" '{print $2}')
	for guide in $(cat guide-proteins.txt)
	do
		count=$(grep -c $guide $outfile)
#		echo "Count: $count" 1>&2
		# If count is 0, guide protein has been lost, exit script
		if [ "$count" -eq 0 ]
		then
			# If count is 0 AND on the first threshold tried, the threshold is too high and needs to be lowered
			if [ "$threshold" -eq "$begin" ]
			then
#				echo "First threshold tried loses protein!" 1>&2
#				difference=$((end-begin))
				if [[ "$((begin-difference))" -le 0 && "$begin" -ne 0 ]]
				then
#					echo "begin-difference is less than or equal to zero and begin is not 0" 1>&2
					if [[ "$difference" -gt 2 ]]
					then
						echo "Your <sweep start> value is too high...Lowering it to $((difference/2))." 1>&2
						echo "high"
						exit 0
					else
						echo "Your <sweep start> value is too high...Lowering it to the lowest possible value: 0." 1>&2
						echo "high"
						exit 0
					fi
				else

					# If first threshold tried loses proteins and the step is one, it means that the ideal threshold is i-1
					if [ "$step" -eq 1 ]
					then
						echo "The first threshold in the sweep loses guide protein ${guide}..." 1>&2
						echo "$threshold"
						exit 0
					fi
					if [ "$begin" -ne 1 ]
					then
						echo "Your <sweep start> value is too high...Lowering it to $((begin-difference))." 1>&2
						echo "high"
						exit 0
					fi
				fi
			else
				echo "Bit score threshold $i loses guide protein ${guide}..." 1>&2
			#	echo "...${guide} is not present when the threshold is ${threshold}." 
				if [ "$step" -ne 1 ]
				then
					echo "Therefore, the ideal threshold will be between $((threshold-step)) and $threshold!" 1>&2
				fi
				echo "$threshold"
				exit 0
			fi
		fi
	done
	echo "Bit score threshold ${threshold} contains all guide proteins..." 1>&2
done
echo "np"
exit 0
