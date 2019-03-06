#!/bin/bash

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
	if [ ! -e "jackhmmer_bs${i}_N$((N-5)).out" ]
	then
		if [ ! -e "jackhmmer_bs${i}_N${N}.out" ]
		then
			echo "COMMAND: jackhmmer --noali -T $i -N $N -o $outfile $AMP $database"
			jackhmmer --noali -T $i -N $N -o $outfile $AMP $database
		fi
	else
		echo "Bit score thresold $i contains all guide proteins..."
		continue
	fi
	# If unconverged, delete the file -- all files that have not been deleted will have been converged/proteins found
	converged=$(grep -c 'CONVERGED' $outfile)
	total=$(grep -c 'Query:' $outfile)
	if [ "$converged" -ne "$total" ]
	then
		echo "At bit score threshold $i, not all queries have converged. Increasing N to $((N+5))."
		rm $outfile
		echo "Exit code: 1"
		exit 1
	fi
	threshold=$i
#	threshold=$(echo $outfile | awk -F "_" '{print $2}' | awk -F "s" '{print $2}')
	for guide in $(cat guide-proteins.txt)
	do
		count=$(grep -c $guide $outfile)
		# If count is 0, guide protein has been lost, exit script
		if [ "$count" -eq 0 ]
		then
			# If count is 0 AND on the first threshold tried, the threshold is too high and needs to be lowered
			if [ "$threshold" -eq "$begin" ]
			then
#				difference=$((end-begin))
				if [[ "$((begin-difference))" -le 0 && "$begin" -ne 0 ]]
				then
					echo "Your <sweep start> value is too high...Lowering it to the lowest possible value: 0. Exit code: 2"
				else

					# If first threshold tried loses proteins and the step is one, it means that the ideal threshold is i-1
					if [ "$step" -eq 1 ]
					then
						if [[ "$threshold" -ge 1 && "$threshold" -le 3 ]]
						then
							echo "Threshold is between 1 and 3."
							echo "Exit code: 0"
							exit 0
						else
							echo "The first threshold in the sweep loses guide protein ${guide}..."
							echo "Exit code: $threshold"
							exit $threshold
						fi
					fi
					if [ "$begin" -ne 0 ]
					then
						echo "Your <sweep start> value is too high...Lowering it to $((begin-difference))."
					fi
				fi
					echo "Exit code: 2"
					exit 2
			else
				echo "Bit score threshold $i loses guide protein ${guide}..."
			#	echo "...${guide} is not present when the threshold is ${threshold}." 
				if [ "$step" -ne 1 ]
				then
					echo "Therefore, the ideal threshold will be between $((threshold-step)) and $threshold!"
				fi
				echo "Exit code: $threshold"
				exit $threshold
			fi
		fi
	done
	echo "Bit score threshold ${threshold} contains all guide proteins..."
done
echo "Exit code: 3"
exit 3
