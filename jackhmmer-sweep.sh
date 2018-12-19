#!/bin/bash

AMP=$1
lit=$2
database=$3
N=$4
begin=$5
end=$6
step=$7
for i in `seq $begin $step $end`
do	
	outfile="jackhmmer_bs${i}_N${N}.out"
	# If N-5 exists, it means the jackhmmer for that iteration already converged
	# So there is no need to run jackhmmer at that threshold again at a higher itereration
	# If N-5 does not exist, either jackhmmer has never been run for that threshold yet, or it had not converged previously.
	if [ ! -e "jackhmmer_bs${i}_N$((N-5)).out" ]
	then
		jackhmmer --noali -T $i -N $N -o $outfile $AMP $database
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
		exit 1
	fi
	threshold=$(echo $outfile | awk -F "_" '{print $2}' | awk -F "s" '{print $2}')
	for guide in $(cat guide-proteins.txt)
	do
		count=$(grep -c $guide $outfile)
		# If count is 0, guide protein has been lost, exit script
		if [ "$count" -eq 0 ]
		then
			# If count is 0 AND on the first threshold tried, the threshold is too high and needs to be lowered
			if [ "$i" -eq "$begin" ]
			then
				echo "Your <sweep start> value is too high...Lowering it by $step."
				exit 2
			else
				echo "Bit score threshold $i loses guide protein ${guide}..."
			#	echo "...${guide} is not present when the threshold is ${threshold}." 
				if [ "$step" -ne 1 ]
				then
					echo "Therefore, the ideal threshold will be bewteen $((threshold-step)) and $threshold!"
				fi
				exit $threshold
			fi
		fi
	done
	echo "Bit score threshold ${threshold} contains all guide proteins..."	
done
exit 0
