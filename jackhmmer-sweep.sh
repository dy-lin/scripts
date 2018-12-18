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
	jackhmmer --noali -T $i -N $N -o $outfile $AMP $database
	converged=$(grep -c 'CONVERGED' $outfile)
	total=$(grep -c 'Query:' $outfile)
	if [ "$converged" -ne "$total" ]
	then
		echo "Not all queries have converged. Increasing N to $((N+5))."
		exit 1
	fi
	threshold=$(echo $outfile | awk -F "_" '{print $2}' | awk -F "s" '{print $2}')
	for guide in $(cat guide-proteins.txt)
	do
		count=$(grep -c $guide $outfile)
		if [ "$count" -eq 0 ]
		then
			echo "...${guide} is not present when the threshold is ${threshold}!"
			exit $threshold
		fi
	done
	echo "Bit score threshold ${threshold} contains all guide proteins..."	
done
exit 0
