#!/bin/bash
set -eu -o pipefail
if [[ "$#" -ne 1 ]]
then
	echo "USAGE: $(basename $0) <genotype>" 1>&2
	exit 1
fi

genotype=$1
if [[ "$genotype" == "PG29" ]]
then
#	tissues="bark"
	tissues="bark embryo flush_bud mature_needle megagametophyte seed_germination xylem young_buds"
	for i in $tissues
	do
		echo -e  "Tissue:\t$i"
		cd "/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/jackhmmer-transcriptome/defensins/${i}/bs103/transcripts/defensins-only/gene_assignment"
		## Get the transcripts queried
		transcripts=$(awk '{print $1}' clustered_transcripts.tsv.blastn | sort -u)
		# READ the BLAST output file
		if [[ ! -e "PG29.gmap.tsv" ]]
		then
			processGMAP.sh PG29.gmap.gff
		fi

		if [[ ! -e "output.scores.psl" ]]
		then
			pslScore output.psl > output.scores.psl
		fi
		for t in $transcripts
		do
			# Get the alignments for each transcript, sorting by percent identity
			blast_highest_pid_line=$(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr | head -n1)
			blast_highest_pid=$(echo "$blast_highest_pid_line" | awk '{print $3}')
			blast_lowest_eval=$(echo "$blast_highest_pid_line" | awk '{print $11}')
			blast_count=$(awk -v var="$blast_highest_pid" -v bar="$blast_lowest_eval" '{if($3==var && $11==bar) print}' <(grep $t clustered_transcripts.tsv.blastn) | wc -l)
			if [[ "$blast_count" -gt 1 ]]
			then
				# do other stuff
				blast="multiple"
			else
				blast=$(echo "$blast_highest_pid_line" | awk '{print $2}' | sed 's/-R.//')
			fi

			# READ from GMAP TSV

			gmap_highest_pid_line=$(grep $t PG29.gmap.tsv | sort -k8,8gr -k3,3gr | head -n1)
			gmap_highest_pid=$(echo "$gmap_highest_pid_line" | awk '{print $3}')
			gmap_highest_cov=$(echo "$gmap_highest_pid_line" | awk '{print $8}')
			gmap_count=$(awk -v var="$gmap_highest_pid" -v bar="$gmap_highest_cov" '{if($3==var && $8==bar) print}' <(grep $t PG29.gmap.tsv) | wc -l)
			if [[ "$gmap_count" -gt 1 ]]
			then
				# do other stuff?
				gmap="multiple"
			else
				gmap=$(echo "$gmap_highest_pid_line" | awk '{print $2}' |sed 's/-R.//')
			fi

		# READ FROM BLAT 
			
			highest_match_line=$(grep $t output.scores.psl | sort -k6,6gr -k5,5gr | head -n1)
			highest_match=$(echo "$highest_match_line" | awk '{print $5}')
			highest_pid=$(echo "$highest_match_line" | awk '{print $6}')
			blat_count=$(awk -v var="$highest_match" -v bar="$highest_pid" '{if($5==var && $6==bar) print}' <(grep $t output.scores.psl) | wc -l )
			if [[ "$blat_count" -gt 1 ]]
			then
				blat="multiple"
			else
				blat=$(echo "$highest_match_line" | awk '{print $1}' | sed 's/-R.//')
			fi

			if [[ "$blast" == "$blat" && "$blast" == "$gmap"  && "$blast" != "multiple" ]]
			then
				echo -e "$t\t$blast"
			elif [[ "$blast" == "multiple" || "$blat" == "multiple" || "$gmap" == "multiple" ]]
			then

				# CHECK IF IT IS TRULY MULTIPLE
#				echo "Multiple check"
				if [[ "$blast" == "multiple" ]]
				then
					unique=$(head -n $blast_count <(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u | wc -l)
					if [[ "$unique" -eq 1 ]]
					then
						blast=$(echo "$blast_highest_pid_line" | awk '{print $2}' | sed 's/-R.//')
					fi
				fi

				if [[ "$gmap" == "multiple" ]]
				then
					unique=$(head -n $gmap_count <(grep $t PG29.gmap.tsv | sort -k8,8gr -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u | wc -l)
					if [[ "$unique" -eq 1 ]]
					then
						gmap=$(echo "$gmap_highest_pid_line" | awk '{print $2}' |sed 's/-R.//')
					fi
				fi

				if [[ "$blat" == "multiple" ]]
				then
					unique=$(head -n $blat_count <(grep $t output.scores.psl | sort -k6,6gr -k5,5gr )| awk '{print $1}' | sed 's/-R.//' | sort -u | wc -l)
					if [[ "$unique" -eq 1 ]]
					then
						blat=$(echo "$highest_match_line" | awk '{print $1}' | sed 's/-R.//')
					fi
				fi

				if [[ "$blast" == "multiple" || "$blat" == "multiple" || "$gmap" == "multiple" ]]
				then
					echo -e "$t\tRESOLVED MULTIPLE - CONFLICT"
					if [[ "$blast" == "multiple" ]]
					then
						blast=$(head -n $blast_count <(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u)
						blast=$(echo "$blast" | tr '\n' ' ')
					fi
						echo -e "\tBLAST:\t$blast"

					if [[ "$gmap" == "multiple" ]]
					then
						gmap=$(head -n $gmap_count <(grep $t PG29.gmap.tsv | sort -k8,8gr -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u)
						gmap=$(echo "$gmap" | tr '\n' ' ')
					fi
						echo -e "\tGMAP:\t$gmap"
					

					if [[ "$blat" == "multiple" ]]
					then
						blat=$(head -n $blat_count <(grep $t output.scores.psl | sort -k6,6gr -k5,5gr )| awk '{print $1}' | sed 's/-R.//' | sort -u)
						blat=$(echo "$blat" | tr '\n' ' ')
					fi
						echo -e "\tBLAT:\t$blat"
					

					final=$(sort <(echo -e "$blast\n$gmap\n$blat") | uniq -c | sort -k1,1 -g -r | head -n 1|  awk '{print $2}')
					echo -e "\tCONSENSUS:\t$final"

				#	echo -e "$t\tUNRESOLVED MULTIPLE"
				#	echo -e "\tBLAST:"
				#	head -n $blast_count <(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr | sed 's/^/\t\t/')
				#	echo -e "\tGMAP:"
				#	head -n $gmap_count <(grep $t PG29.gmap.tsv | sort -k8,8gr -k3,3gr | sed 's/^/\t\t/')
				#	echo -e "\tBLAT:"
				#	head -n $blat_count <(grep $t output.psl | sort -k1,1 -g -r | sed 's/^/\t\t/')
				else
					echo -e "$t\tRESOLVED MULTIPLE - CONFLICT"
					echo -e "\tBLAST:\t$blast"
					echo -e "\tBLAT:\t$blat"
					echo -e "\tGMAP:\t$gmap"
					# Find consensus
					highest_count=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $1}')
					# Since there are only 3 items to be sort/uniq, the combos are 1 1 1 or 2 1.
					if [[ "$highest_count" -eq 1 ]]
					then
						echo -e "\tCONSENSUS:\tunresolved"
					else
						consensus=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $2}')
						echo -e "\tCONSENSUS:\t$consensus"
					fi
				fi
			else
				echo -e "$t\tCONFLICT"
				echo -e "\tBLAST:\t$blast"
				echo -e "\tBLAT:\t$blat"
				echo -e "\tGMAP:\t$gmap"
				# Find consensus
				highest_count=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $1}')
				# Since there are only 3 items to be sort/uniq, the combos are 1 1 1 or 2 1.
				if [[ "$highest_count" -eq 1 ]]
				then
					echo -e "\tCONSENSUS:\tunresolved"
				else
					consensus=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $2}')
					echo -e "\tCONSENSUS:\t$consensus"
				fi
				
			fi
		done
		echo
	done
else
	for i in "Dev_SC" "Cort_Par"
	do
		echo "Tissue: $i"
		cd "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/jackhmmer-transcriptome/defensins/stonecell/${i}/bs13/transcripts/defensins-only/gene_assignment"
		# do stuff
		transcripts=$(awk '{print $1}' clustered_transcripts.tsv.blastn | sort -u)
		# READ the BLAST output file
		if [[ ! -e "Q903.gmap.tsv" ]]
		then
			processGMAP.sh Q903.gmap.gff
		fi

		if [[ ! -e "output.scores.psl" ]]
		then
			pslScore output.psl > output.scores.psl
		fi
		for t in $transcripts
		do
			# Get the alignments for each transcript, sorting by percent identity
			blast_highest_pid_line=$(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr | head -n1)
			blast_highest_pid=$(echo "$blast_highest_pid_line" | awk '{print $3}')
			blast_lowest_eval=$(echo "$blast_highest_pid_line" | awk '{print $11}')
			blast_count=$(awk -v var="$blast_highest_pid" -v bar="$blast_lowest_eval" '{if($3==var && $11==bar) print}' <(grep $t clustered_transcripts.tsv.blastn) | wc -l)
			if [[ "$blast_count" -gt 1 ]]
			then
				# do other stuff
				blast="multiple"
			else
				blast=$(echo "$blast_highest_pid_line" | awk '{print $2}' | sed 's/-R.//')
			fi

			# READ from GMAP TSV

			gmap_highest_pid_line=$(grep $t Q903.gmap.tsv | sort -k8,8gr -k3,3gr | head -n1)
			gmap_highest_pid=$(echo "$gmap_highest_pid_line" | awk '{print $3}')
			gmap_highest_cov=$(echo "$gmap_highest_pid_line" | awk '{print $8}')
			gmap_count=$(awk -v var="$gmap_highest_pid" -v bar="$gmap_highest_cov" '{if($3==var && $8==bar) print}' <(grep $t Q903.gmap.tsv) | wc -l) 
			if [[ "$gmap_count" -gt 1 ]]
			then
				# do other stuff?
				gmap="multiple"
			else
				gmap=$(echo "$gmap_highest_pid_line" | awk '{print $2}' | sed 's/-R.//')
			fi

		# READ FROM BLAT 

			highest_match_line=$(grep $t output.scores.psl | sort -k1,1 -g -r | head -n1)
			highest_match=$(echo "$highest_match_line" | awk '{print $5}')
			highest_pid=$(echo "$highest_match_line" | awk '{print $6}')
			blat_count=$(awk -v var="$highest_match" -v bar="$highest_pid" '{if($5==var && $6==bar) print}' <(grep $t output.scores.psl) | wc -l)
			if [[ "$blat_count" -gt 1 ]]
			then
				blat="multiple"
			else
				blat=$(echo "$highest_match_line" | awk '{print $1}' | sed 's/-R.//')
			fi

			if [[ "$blast" == "$blat" && "$blast" == "$gmap"  && "$blast" != "multiple" ]]
			then
				echo -e "$t\t$blast"
			elif [[ "$blast" == "multiple" || "$blat" == "multiple" || "$gmap" == "multiple" ]]
			then
				# CHECK IF IT IS TRULY MULTIPLE
#				echo "Multiple check"
				if [[ "$blast" == "multiple" ]]
				then
					unique=$(head -n $blast_count <(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u | wc -l)
					if [[ "$unique" -eq 1 ]]
					then
						blast=$(echo "$blast_highest_pid_line" | awk '{print $2}' | sed 's/-R.//')
					fi
				fi

				if [[ "$gmap" == "multiple" ]]
				then
					unique=$(head -n $gmap_count <(grep $t Q903.gmap.tsv | sort -k8,8gr -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u | wc -l)
					if [[ "$unique" -eq 1 ]]
					then
						gmap=$(echo "$gmap_highest_pid_line" | awk '{print $2}' |sed 's/-R.//')
					fi
				fi

				if [[ "$blat" == "multiple" ]]
				then
					unique=$(head -n $blat_count <(grep $t output.scores.psl | sort -k6,6gr -k5,5gr)| awk '{print $1}' | sed 's/-R.//' | sort -u | wc -l)
					if [[ "$unique" -eq 1 ]]
					then
						blat=$(echo "$highest_match_line" | awk '{print $1}' | sed 's/-R.//')
					fi
				fi

				if [[ "$blast" == "multiple" || "$blat" == "multiple" || "$gmap" == "multiple" ]]
				then

					echo -e "$t\tRESOLVED MULTIPLE - CONFLICT"
					if [[ "$blast" == "multiple" ]]
					then
						blast=$(head -n $blast_count <(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u)
						blast=$(echo "$blast" | tr '\n' ' ')
					fi
						echo -e "\tBLAST:\t$blast"

					if [[ "$gmap" == "multiple" ]]
					then
						gmap=$(head -n $gmap_count <(grep $t Q903.gmap.tsv | sort -k8,8gr -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u)
						gmap=$(echo "$gmap" | tr '\n' ' ')
					fi
						echo -e "\tGMAP:\t$gmap"
					

					if [[ "$blat" == "multiple" ]]
					then
						blat=$(head -n $blat_count <(grep $t output.scores.psl | sort -k6,6gr -k5,5gr) | awk '{print $1}' | sed 's/-R.//' | sort -u)
						blat=$(echo "$blat" | tr '\n' ' ')
					fi
						echo -e "\tBLAT:\t$blat"
					

					final=$(sort <(echo -e "$blast\n$gmap\n$blat") | uniq -c | sort -k1,1 -g -r | head -n 1 | awk '{print $2}')
					echo -e "\tCONSENSUS:\t$final"

				#	echo -e "$t\tUNRESOLVED MULTIPLE"
				#	echo -e "\tBLAST:"
				#	head -n $blast_count <(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr | sed 's/^/\t\t/')
				#	echo -e "\tGMAP:"
				#	head -n $gmap_count <(grep $t Q903.gmap.tsv | sort -k8,8gr -k3,3gr | sed 's/^/\t\t/')
				#	echo -e "\tBLAT:"
				#	head -n $blat_count <(grep $t output.psl | sort -k1,1 -g -r | sed 's/^/\t\t/')
				else
					echo -e "$t\tRESOLVED MULTIPLE - CONFLICT"
					echo -e "\tBLAST:\t$blast"
					echo -e "\tBLAT:\t$blat"
					echo -e "\tGMAP:\t$gmap"
					# Find consensus
					highest_count=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $1}')
					# Since there are only 3 items to be sort/uniq, the combos are 1 1 1 or 2 1.
					if [[ "$highest_count" -eq 1 ]]
					then
						echo -e "\tCONSENSUS:\tunresolved"
					else
						consensus=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $2}')
						echo -e "\tCONSENSUS:\t$consensus"
					fi
				fi

			else
				echo -e "$t\tCONFLICT"
				echo -e "\tBLAST:\t$blast"
				echo -e "\tBLAT:\t$blat"
				echo -e "\tGMAP:\t$gmap"
				# Find consensus
				highest_count=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $1}')
				# Since there are only 3 items to be sort/uniq, the combos are 1 1 1 or 2 1.
				if [[ "$highest_count" -eq 1 ]]
				then
					echo -e "\tCONSENSUS:\tunresolved."
				else
					consensus=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $2}')
					echo -e "\tCONSENSUS:\t$consensus"
				fi

			fi
		done
		echo
	done
	for i in "control" "gallery" "wound"
	do
		echo "Tissue: $i"
		cd "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/jackhmmer-transcriptome/defensins/wpw/${i}/bs13/transcripts/defensins-only/gene_assignment"
		# do stuff
		transcripts=$(awk '{print $1}' clustered_transcripts.tsv.blastn | sort -u)
		# READ the BLAST output file
		if [[ ! -e "Q903.gmap.tsv" ]]
		then
			processGMAP.sh Q903.gmap.gff
		fi

		if [[ ! -e "output.scores.psl" ]]
		then
			pslScore output.psl > output.scores.psl
		fi
		for t in $transcripts
		do
			# Get the alignments for each transcript, sorting by percent identity
			blast_highest_pid_line=$(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr | head -n1)
			blast_highest_pid=$(echo "$blast_highest_pid_line" | awk '{print $3}')
			blast_lowest_eval=$(echo "$blast_highest_pid_line" | awk '{print $11}')
			blast_count=$(awk -v var="$blast_highest_pid" -v bar="$blast_lowest_eval" '{if($3==var && $11==bar) print}' <(grep $t clustered_transcripts.tsv.blastn) | wc -l)
			if [[ "$blast_count" -gt 1 ]]
			then
				# do other stuff
				blast="multiple"
			else
				blast=$(echo "$blast_highest_pid_line" | awk '{print $2}' | sed 's/-R.//')
			fi

			# READ from GMAP TSV

			gmap_highest_pid_line=$(grep $t Q903.gmap.tsv | sort -k8,8gr -k3,3gr | head -n1)
			gmap_highest_pid=$(echo "$gmap_highest_pid_line" | awk '{print $3}')
			gmap_highest_cov=$(echo "$gmap_highest_pid_line" | awk '{print $8}')
			gmap_count=$(awk -v var="$gmap_highest_pid" -v bar="$gmap_highest_cov" '{if($3==var && $8==bar) print}'  <(grep $t Q903.gmap.tsv) | wc -l)
			if [[ "$gmap_count" -gt 1 ]]
			then
				# do other stuff?
				gmap="multiple"
			else
				gmap=$(echo "$gmap_highest_pid_line" | awk '{print $2}' | sed 's/-R.//' )
			fi

		# READ FROM BLAT 

			highest_match_line=$(grep $t output.scores.psl | sort -k6,6gr -k5,5gr | head -n1)
			highest_match=$(echo "$highest_match_line" | awk '{print $5}')
			highest_pid=$(echo "$highest_match_line" | awk '{print $6}')
			blat_count=$(awk -v var="$highest_match" -v bar="$highest_pid" '{if($5==var && $6==bar) print}' <(grep $t output.scores.psl)| wc -l)
			if [[ "$blat_count" -gt 1 ]]
			then
				blat="multiple"
			else
				blat=$(echo "$highest_match_line" | awk '{print $1}' | sed 's/-R.//')
			fi

			if [[ "$blast" == "$blat" && "$blast" == "$gmap"  && "$blast" != "multiple" ]]
			then
				echo -e "$t\t$blast"
			elif [[ "$blast" == "multiple" || "$blat" == "multiple" || "$gmap" == "multiple" ]]
			then
				# CHECK IF IT IS TRULY MULTIPLE
#				echo "Multiple check"
				if [[ "$blast" == "multiple" ]]
				then
					unique=$(head -n $blast_count <(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u | wc -l)
					if [[ "$unique" -eq 1 ]]
					then
						blast=$(echo "$blast_highest_pid_line" | awk '{print $2}' | sed 's/-R.//')
					fi
				fi

				if [[ "$gmap" == "multiple" ]]
				then
					unique=$(head -n $gmap_count <(grep $t Q903.gmap.tsv | sort -k8,8gr -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u | wc -l)
					if [[ "$unique" -eq 1 ]]
					then
						gmap=$(echo "$gmap_highest_pid_line" | awk '{print $2}' |sed 's/-R.//')
					fi
				fi

				if [[ "$blat" == "multiple" ]]
				then
					unique=$(head -n $blat_count <(grep $t output.scores.psl | sort -k6,6gr -k5,5gr)| awk '{print $1}' | sed 's/-R.//' | sort -u | wc -l)
					if [[ "$unique" -eq 1 ]]
					then
						blat=$(echo "$highest_match_line" | awk '{print $1}' | sed 's/-R.//')
					fi
				fi

				if [[ "$blast" == "multiple" || "$blat" == "multiple" || "$gmap" == "multiple" ]]
				then
					echo -e "$t\tRESOLVED MULTIPLE - CONFLICT"
					if [[ "$blast" == "multiple" ]]
					then
						blast=$(head -n $blast_count <(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u)
						blast=$(echo "$blast" | tr '\n' ' ')
					fi
						echo -e "\tBLAST:\t$blast"

					if [[ "$gmap" == "multiple" ]]
					then
						gmap=$(head -n $gmap_count <(grep $t Q903.gmap.tsv | sort -k8,8gr -k3,3gr) | awk '{print $2}' | sed 's/-R.//' | sort -u)
						gmap=$(echo "$gmap" | tr '\n' ' ')
					fi
						echo -e "\tGMAP:\t$gmap"
					

					if [[ "$blat" == "multiple" ]]
					then
						blat=$(head -n $blat_count <(grep $t output.scores.psl | sort -k6,6gr -k5,5gr )| awk '{print $1}' | sed 's/-R.//' | sort -u)
						blat=$(echo "$blat" | tr '\n' ' ')
					fi
						echo -e "\tBLAT:\t$blat"
					

					final=$(sort <(echo -e "$blast\n$gmap\n$blat") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $2}')
					echo -e "\tCONSENSUS:\t$final"

#					if [[ "$blast" == "multiple" || "$blat" == "multiple" || "$gmap" == "multiple" ]]
#					then
#						echo -e "$t\tUNRESOLVED MULTIPLE"
#						echo -e "\tBLAST:"
#						head -n $blast_count <(grep $t clustered_transcripts.tsv.blastn | sort -k11,11g -k3,3gr | sed 's/^/\t\t/')
#						echo -e "\tGMAP:"
#						head -n $gmap_count <(grep $t Q903.gmap.tsv | sort -k8,8gr -k3,3gr | sed 's/^/\t\t/')
#						echo -e "\tBLAT:"
#						head -n $blat_count <(grep $t output.scores.psl | sort -k6,6gr -k5,5gr | sed 's/^/\t\t/')
#					fi
				else
					echo -e "$t\tRESOLVED MULTIPLE - CONFLICT"
					echo -e "\tBLAST:\t$blast"
					echo -e "\tBLAT:\t$blat"
					echo -e "\tGMAP:\t$gmap"
					# Find consensus
					highest_count=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $1}')
					# Since there are only 3 items to be sort/uniq, the combos are 1 1 1 or 2 1.
					if [[ "$highest_count" -eq 1 ]]
					then
						echo -e "\tCONSENSUS:\tunresolved"
					else
						consensus=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $2}')
						echo -e "\tCONSENSUS:\t$consensus"
					fi
				fi
			else
				echo -e "$t\tCONFLICT"
				echo -e "\tBLAST:\t$blast"
				echo -e "\tBLAT:\t$blat"
				echo -e "\tGMAP:\t$gmap"
				# Find consensus
				highest_count=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $1}')
				# Since there are only 3 items to be sort/uniq, the combos are 1 1 1 or 2 1.
				if [[ "$highest_count" -eq 1 ]]
				then
					echo -e "\tCONSENSUS:\tunresolved."
				else
					consensus=$(sort <(echo -e "$blast\n$blat\n$gmap") | uniq -c | sort -k1,1 -g -r | head -n1 | awk '{print $2}')
					echo -e "\tCONSENSUS:\t$consensus"
				fi

			fi
		done
		echo
	done
fi

