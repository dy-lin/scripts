#!/bin/bash
PROGRAM=$(basename $0)
threads=48
gethelp=false
while getopts :ht $opt
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

dbtype=$stype
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
	tblastn) qtype=prot; stype=nucl;dbtype=$stype;if [[ "$(ls ${dir}/${name}*.nhr 2> /dev/null | wc -l)" -eq 0 ]];then makeblastdb -in $db -out ${dir}/$sname -dbtype $dbtype; fi;echo "Running tblastn..." 1>&2;tblastn -query $query -db $dir/$sname -out $dir/${qname}.aln.tblastn -num_threads $threads;tblastn -query $query -db $dir/$sname -out $dir/${qname}.tsv.tblastn -num_threads $threads -outfmt '6 std qcovs';;
blastn) qtype=nucl; stype=nucl;dbtype=$stype;if [[ "$(ls ${dir}/${name}*.nhr 2> /dev/null | wc -l)" -eq 0 ]];then makeblastdb -in $db -out ${dir}/$sname -dbtype $dbtype; fi;echo "Running blastn..." 1>&2;blastn -query $query -db $dir/$sname -out $dir/${qname}.aln.blastn -num_threads $threads;blastn -query $query -db $dir/$sname -out $dir/${qname}.tsv.blastn -num_threads $threads -outfmt '6 std qcovs';;
	blastp) qtype=prot; stype=prot;dbtype=$stype;if [[ "$(ls ${dir}/${name}*.phr 2> /dev/null | wc -l)" -eq 0 ]];then makeblastdb -in $db -out ${dir}/$sname -dbtype $dbtype; fi;echo "Running blastp..." 1>&2;blastp -query $query -db $dir/$sname -out $dir/${qname}.aln.tblastn -num_threads $threads;blastp -query $query -db $dir/$sname -out $dir/${qname}.tsv.tblastn -num_threads $threads -outfmt '6 std qcovs';;
blastx) qtype=nucl; stype=prot;dbtype=$stype;if [[ "$(ls ${dir}/${name}*.phr 2> /dev/null | wc -l)" -eq 0 ]];then makeblastdb -in $db -out ${dir}/$sname -dbtype $dbtype; fi;echo "Running blastx..." 1>&2;blastx -query $query -db $dir/$sname -out $dir/${qname}.aln.blastx -num_threads $threads;blastx -query $query -db $dir/$sname -out $dir/${qname}.tsv.blastx -num_threads $threads -outfmt '6 std qcovs';;
	tblastx) qtype=prot; stype=prot;dbtype=$stype;if [[ "$(ls ${dir}/${name}*.phr 2> /dev/null | wc -l)" -eq 0 ]];then makeblastdb -in $db -out ${dir}/$sname -dbtype $dbtype; fi;echo "Running tblastx..." 1>&2;tblastx -query $query -db $dir/$sname -out $dir/${qname}.aln.tblastx -num_threads $threads;tblastx -query $query -db $dir/$sname -out $dir/${qname}.tsv.tblastx -num_threads $threads -outfmt '6 std qcovs';;
	\?) echo "ERROR: Unrecognized BLAST program." 1>&2 ; exit 1;;
esac


