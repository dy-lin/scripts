#!/bin/bash
set -eo pipefail
PROGRAM=$(basename $0)
threads=48
gethelp=false

if [[  "$#" -eq 0 ]]
then
	echo "USAGE: $PROGRAM <BLAST program> <query> <subject> [output directory]" 1>&2
	echo "DESCRIPTION: Blasts the query file against the subject file." 1>&2
	echo -e "OPTIONS:\n\t-h\tShow help menu\n\t-t\tnumber of threads" 1>&2
	exit 1
fi
while getopts :ht: opt
do
	case $opt in
		h) gethelp=true;;
		t) threads=$OPTARG;;
		\?) echo "ERROR: $PROGRAM: Invalid option $OPTARG" 1>&2; exit 1;;
	esac
done
shift $((OPTIND-1))

if [[ "$#" -ne 3 && "$#" -ne 4 ]] || [[ "$gethelp" = true ]]
then
	echo "USAGE: $PROGRAM <BLAST program> <query> <subject> [output directory]" 1>&2
	echo "DESCRIPTION: Blasts the query file against the subject file." 1>&2
	echo -e "OPTIONS:\n\t-h\tShow help menu\n\t-t\tnumber of threads" 1>&2
	exit 1
fi
blast=$1

query=$2
db=$3
qname=$(basename ${query%.*})

sname=$(basename ${db%.*})
if [[ -z "$4" ]]
then
	dir=$(pwd)
else
	dir=$4
fi

if [[ "$db" == nr* ]]
then
	echo "BLASTing against non-redundant protein sequences..." 1>&2
#	export BLASTDB=/projects/btl/db/blast-20190104
	export BLASTDB=/projects/btl/dlin/datasets/blast
	case $blast in 
		blastp) qtype=prot;stype=prot;dbtpe=$stype; blastp -query $query -db nr_v5 -out $dir/${qname}.nr.aln.blastp -num_threads $threads; blastp -query $query -db nr_v5 -out $dir/${qname}.nr.tsv.blastp -num_threads $threads -outfmt '7 std qcovs' ;;
		blastx) qtype=nucl;stype=prot; dbtype=$stype; blastx -query $query -db nr_v5 -out $dir/${qname}.nr.aln.blastx -num_threads $threads; blastx -query $query -db nr_v5 -out $dir/${qname}.nr.tsv.blastx -num_threads $threads -outfmt '7 std qcovs' ;;
		tblastx) qtype=prot; stype=prot; dbtype=$stype; tblastx -query $query -db nr_v5 -out $dir/${qname}.nr.aln.tblastx -num_threads $threads; tblastx -query $query -db nr_v5 -out $dir/${qname}.nr.tsv.tblastx -num_threads $threads -outfmt '7 std qcovs' ;;
		\?) echo "ERROR: Unrecognized BLAST program." 1>&2; exit 1;;
	esac
elif [[ "$db" == "sprot" ]]
then
	while [[ "$blast" == *n ]]
	do
		echo "Swissprot is a protein database. BLAST program must be one of blastn or tblastn." 1>&2
		read -p "BLAST program: " blast
	done
	echo "BLASTing against Swissprot protein sequences..." 1>&2
	export BLASTDB=/projects/btl/dlin/datasets/sprot
	case $blast in 
		blastp) qtype=prot;stype=prot;dbtpe=$stype; blastp -query $query -db uniprot_sprot -out $dir/${qname}.sprot.aln.blastp -num_threads $threads; blastp -query $query -db uniprot_sprot -out $dir/${qname}.sprot.tsv.blastp -num_threads $threads -outfmt '7 std qcovs' ;;
		blastx) qtype=nucl;stype=prot; dbtype=$stype; blastx -query $query -db uniprot_sprot -out $dir/${qname}.sprot.aln.blastx -num_threads $threads; blastx -query $query -db uniprot_sprot -out $dir/${qname}.sprot.tsv.blastx -num_threads $threads -outfmt '7 std qcovs' ;;
		tblastx) qtype=prot; stype=prot; dbtype=$stype; tblastx -query $query -db uniprot_sprot -out $dir/${qname}.sprot.aln.tblastx -num_threads $threads; tblastx -query $query -db uniprot_sprot -out $dir/${qname}.sprot.tsv.tblastx -num_threads $threads -outfmt '7 std qcovs' ;;
	blastp-short) qtype=prot; stype=prot;dbtype=$stype; blastp -query $query -db uniprot_sprot -out $dir/${qname}.sprot.short.aln.blastp -num_threads $threads -task blastp-short; blastp -query $query -db uniprot_sprot -out $dir/${qname}.sprot.short.tsv.blastp -num_threads $threads -task blastp-short -outfmt '7 std qcovs';;
		?) echo "ERROR: Unrecognized BLAST program." 1>&2; exit 1;;
	esac
elif [[ "$db" == "trembl" ]]
then
	while [[ "$blast" == *n ]]
	do
		echo "TREMBL is a protein database. BLAST program must be one of blastp, blastx, or tblastx." 1>&2
		read -p "BLAST program: " blast
	done
	echo "BLASTing against TREMBL protein sequences..." 1>&2
	export BLASTDB=/projects/btl/dlin/datasets/trembl
	case $blast in 
		blastp) qtype=prot;stype=prot;dbtpe=$stype; blastp -query $query -db uniprot_trembl -out $dir/${qname}.trembl.aln.blastp -num_threads $threads; blastp -query $query -db uniprot_trembl -out $dir/${qname}.trembl.tsv.blastp -num_threads $threads -outfmt '7 std qcovs' ;;
		blastx) qtype=nucl;stype=prot; dbtype=$stype; blastx -query $query -db uniprot_trembl -out $dir/${qname}.trembl.aln.blastx -num_threads $threads; blastx -query $query -db uniprot_trembl -tout $dir/${qname}.trembl.tsv.blastx -num_threads $threads -outfmt '7 std qcovs' ;;
		tblastx) qtype=prot; stype=prot; dbtype=$stype; tblastx -query $query -db uniprot_trembl -out $dir/${qname}.trembl.aln.tblastx -num_threads $threads; tblastx -query $query -db uniprot_trembl -out $dir/${qname}.trembl.tsv.tblastx -num_threads $threads -outfmt '7 std qcovs' ;;
	blastp-short) qtype=prot; stype=prot; dbtype=$stpe; blastp -query $query -db uniprot_trembl -out $dir/${qname}.trembl.aln.blastp -num_threads $threads -task blastp-short; blastp -query $query -db uniprot_trembl -out $dir/${qname}.trembl.tsv.blastp -num_threads $threads -outfmt '7 std qcovs' -task blastp-short;;
		\?) echo "ERROR: Unrecognized BLAST program." 1>&2; exit 1;;
	esac

elif [[ "$db" == "RefSeq_Plants" ]]
then
	echo "BLASTing against RefSeq Plant sequences..." 1>&2
	export BLASTDB=/projects/btl/dlin/datasets/RefSeq_Plants
	case $blast in
		blastp) qtype=prot;stype=prot;dbtype=$stype; blastp -query $query -db plant.all.protein -out $dir/${qname}.RefSeq_Plants.aln.blastp -num_threads $threads; blastp -query $query -db plant.all.protein -out $dir/${qname}.RefSeq_Plants.tsv.blastp -num_threads $threads -outfmt '7 std qcovs' ;;
	blastp-short) qtype=prot;stype=prot;dbtype=$stype; blastp -query $query -db plant.all.protein -out $dir/${qname}.RefSeq_Plants.short.aln.blastp -num_threads $threads -task blastp-short; blastp -query $query -db plant.all.protein -out $dir/${qname}.RefSeq_Plants.short.tsv.blastp -num_threads $threads -task blastp-short -outfmt '7 std qcovs' ;;
		blastn) read -p "DNA or RNA? " answer
			case $answer in
				[Rr][Nn][Aa]) qtype=nucl;stype=prot;dbtype=$stype; blastn -query $query -db plant.all.rna -out $dir/${qname}.RefSeq_Plants.rna.aln.blastn -num_threads $threads; blastn -query $query -db plant.all.rna -out $dir/${qname}.RefSeq_Plants.rna.tsv.blastn -num_threads $threads -outfmt '7 std qcovs' ;;
				[Dd][Nn][Aa]) qtype=nucl;stype=prot;dbtype=$stype; blastn -query $query -db plant.all.dna -out $dir/${qname}.RefSeq_Plants.dna.aln.blastn -num_threads $threads; blastn -query $query -db plant.all.dna -out $dir/${qname}.RefSeq_Plants.dna.tsv.blastn -num_threads $threads -outfmt '7 std qcovs' ;;
				\?) echo "Invalid nucleotide type." 1>&2 ; exit 1;;
			esac
		#	\?) echo "ERROR: Unrecognized BLAST program." 1>&2; exit 1;;
	esac
elif [[ "$db" == nt* ]]
then
	echo "BLASTing against non-redundant nucleotide sequences..." 1>&2
	export BLASTDB=/projects/btl/dlin/datasets/blast
	case $blast in
		blastn) qtype=nucl;stype=nucl;dbtype=$stype; blastn -query $query -db nt_v5 -out $dir/${qname}.nt.aln.blastn -num_threads $threads; blastn -query $query -db nt_v5 -out $dir/${qname}.nt.tsv.blastn -num_threads $threads -outfmt '7 std qcovs' ;;
		tblastn) qtype=prot;stype=nucl;dbtype=$stype; tblastn -query $query -db nt_v5 -out $dir/${qname}.nt.aln.tblastn -num_threads $threads; tblastn -query $query -db nt_v5 -out $dir/${qname}.nt.tsv.tblastn -num_threads $threads -outfmt '7 std qcovs' ;;
		\?) echo "ERROR: Unrecognized BLAST program." 1>&2; exit 1;;
	esac
else
	sname=$(basename ${db%.*})
	# qdetected=$(findType.sh $query)

	# if [[ "$qdetected" != "$qtype" ]]
	# then
	#	echo "Detected query file type is: $qdetected, but given $qtype." 1>&2
	#	read -p "Continue with $dbtype? (Y/N): " answer
	#	if [[ "$answer" =! [Nn] ]] 
	#	then
	#		exit 1
	#	fi
	#	qdetected=$qtype
	# fi

	case $blast in
		tblastn) qtype=prot; stype=nucl;dbtype=$stype;if [[ "$(ls ${dir}/${sname}*.nhr 2> /dev/null | wc -l)" -eq 0 ]];then makeblastdb -in $db -out ${dir}/$sname -dbtype $dbtype; fi;echo "Running tblastn..." 1>&2;tblastn -query $query -db $dir/$sname -out $dir/${qname}.aln.tblastn -num_threads $threads;tblastn -query $query -db $dir/$sname -out $dir/${qname}.tsv.tblastn -num_threads $threads -outfmt '7 std qcovs';;
	blastn) qtype=nucl; stype=nucl;dbtype=$stype;if [[ "$(ls ${dir}/${sname}*.nhr 2> /dev/null | wc -l)" -eq 0 ]];then makeblastdb -in $db -out ${dir}/$sname -dbtype $dbtype; fi;echo "Running blastn..." 1>&2;blastn -query $query -db $dir/$sname -out $dir/${qname}.aln.blastn -num_threads $threads;blastn -query $query -db $dir/$sname -out $dir/${qname}.tsv.blastn -num_threads $threads -outfmt '7 std qcovs';;
		blastp) qtype=prot; stype=prot;dbtype=$stype;if [[ "$(ls ${dir}/${sname}*.phr 2> /dev/null | wc -l)" -eq 0 ]];then makeblastdb -in $db -out ${dir}/$sname -dbtype $dbtype; fi;echo "Running blastp..." 1>&2;blastp -query $query -db $dir/$sname -out $dir/${qname}.aln.blastp -num_threads $threads;blastp -query $query -db $dir/$sname -out $dir/${qname}.tsv.blastp -num_threads $threads -outfmt '7 std qcovs';;
	blastx) qtype=nucl; stype=prot;dbtype=$stype;if [[ "$(ls ${dir}/${sname}*.phr 2> /dev/null | wc -l)" -eq 0 ]];then makeblastdb -in $db -out ${dir}/$sname -dbtype $dbtype; fi;echo "Running blastx..." 1>&2;blastx -query $query -db $dir/$sname -out $dir/${qname}.aln.blastx -num_threads $threads;blastx -query $query -db $dir/$sname -out $dir/${qname}.tsv.blastx -num_threads $threads -outfmt '7 std qcovs';;
		tblastx) qtype=prot; stype=prot;dbtype=$stype;if [[ "$(ls ${dir}/${sname}*.phr 2> /dev/null | wc -l)" -eq 0 ]];then makeblastdb -in $db -out ${dir}/$sname -dbtype $dbtype; fi;echo "Running tblastx..." 1>&2;tblastx -query $query -db $dir/$sname -out $dir/${qname}.aln.tblastx -num_threads $threads;tblastx -query $query -db $dir/$sname -out $dir/${qname}.tsv.tblastx -num_threads $threads -outfmt '7 std qcovs';;
		\?) echo "ERROR: Unrecognized BLAST program." 1>&2 ; exit 1;;
	esac

fi
