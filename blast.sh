#!/bin/bash
PROGRAM=$(basename $0)

if [[ "$#" -ne 3 ]]
then
	echo "USAGE: $PROGRAM <nucl or prot> <query> <subject>" 1>&2
	echo "DESCRIPTION: Blasts the query file against the subject file." 1>&2
	exit 1
fi

dbtype=$1
query=$2
db=$3
qname=${query%.*}
sname=${db%.*}
## Make sure that the subject sequences given are actually the type they claim to be
sdetected=$(findType.sh $db)
if [[ "$sdetected" == "$dbtype" ]]
then
	makeblastdb -in $db -out $sname -dbtype $dbtype
else
	echo "Detected subject file type is: $sdetected, but given $dbtype." 1>&2
	exit 1
fi

qdetected=$(findType.sh $query)

echo "Query file type detected: $qdetected" 1>&2
echo "Subject file type detected: $sdetected" 1>&2
if [[ "$qdetected" == "nucl" && "$dbtype" == "nucl" ]]
then
	echo "Running blastn..." 1>&2
	blastn -query $query -db $sname -out ${qname}.aln.blastn -num_threads 48
	blastn -query $query -db $sname -out ${qname}.tsv.blastn -num_threads 48 -outfmt '6 std qcovs'
fi

if [[ "$qdetected" == "prot" && "$dbtype" == "nucl" ]]
then
	echo "Running tblastn..." 1>&2
	tblastn -query $query -db $sname -out ${qname}.aln.tblastn -num_threads 48
	tblastn -query $query -db $sname -out ${qname}.tsv.tblastn -num_threads 48 -outfmt '6 std qcovs'
fi

if [[ "$qdetected" == "prot" && "$dbtype" == "prot" ]]
then
	echo "Running blastp..." 1>&2
	blastp -query $query -db $sname -out ${qname}.aln.tblastn -num_threads 48
	blastp -query $query -db $sname -out ${qname}.tsv.tblastn -num_threads 48 -outfmt '6 std qcovs'
fi

if [[ "$qdetected" == "nucl" && "$dbtype" == "prot" ]]
then
	echo "Running blastx..." 1>&2
	blastx -query $query -db $sname -out ${qname}.aln.blastx -num_threads 48
	blastx -query $query -db $sname -out ${qname}.tsv.blastx -num_threads 48 -outfmt '6 std qcovs'
fi
