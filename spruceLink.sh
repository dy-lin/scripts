#!/bin/bash
spruce=$1
if [[ "$#" -ne 1 && "$#" -ne 2 ]]
then
	echo "USAGE: $(basename $0) <spruce genotype" 1>&2
	echo "DESCRIPTION: Creates softlinks for all spruce files in the specified directories." 1>&2
	exit 1
fi

if [[ -z $2 ]]
then
	dir=$(pwd)
else
	dir=$2
fi

if [[ "$spruce" == "WS77111" ]]
then
	ln -s /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/uniprot-ncbi-filtered-labeled-nonSubsumed-9April2018.fa ${dir}/literature-AMPs.faa
	ln -s /projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/GeneAndModelsPredFirstStepMaker/FourthStepCumulat/OutputMaker/LoweAED/WS77111_ntEdit_selShuf.all.polishedGeneseAED.gff ${dir}/genes.gff
	ln -s /projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/GeneAndModelsPredFirstStepMaker/FourthStepCumulatNonNtEdit/OutputMaker/LoweAED/WS77111.all.maker.proteinsLoweAEDpolished.fasta ${dir}/proteins.faa
	ln -s /projects/spruceup/pglauca/WS77111/assemblies/releases/version2/WS77111v2_release/WS77111-v2_1000plus_LGs.fa ${dir}/scaffolds.fa
	ln -s /projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/GeneAndModelsPredFirstStepMaker/FourthStepCumulatNonNtEdit/OutputMaker/LoweAED/WS77111.all.maker.transcriptsLoweAEDpolished.fasta ${dir}/transcripts.fa
elif [[ "$spruce" == "PG29" ]]
then
	ln -s /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/uniprot-ncbi-filtered-labeled-nonSubsumed-9April2018.fa ${dir}/literature-AMPs.faa
	ln -s /projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/ThirdIterationCumulatNonntEdit/OutputMaker/LoweAED/PG29v5.all.maker.proteins1LoweAEDnoStopPolished.fasta ${dir}/proteins.faa
	ln -s /projects/spruceup/interior_spruce/PG29/assemblies/releases/PG29-v5/PG29-v5_1000plus.fa ${dir}/scaffolds.fa
	ln -s /projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/ThirdIterationCumulatNonntEdit/OutputMaker/LoweAED/PG29v5.all.maker.transcripts1LoweAED.fasta ${dir}/transcripts.fa
	ln -s /projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/ThirdIterationCumulatNonntEdit/OutputMaker/LoweAED/PG29.all.polished.genesLoweAED.gff ${dir}/genes.gff
elif [[ "$spruce" == "Q903" ]]
then
	ln -s /projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/ThirdStepCumulat/OutputMaker/LoweAED/Q903.all.maker.proteinsLoweAED.fasta ${dir}/proteins.faa
	ln -s /projects/spruceup/psitchensis/Q903/assembly/releases/v0.1beta/Q903_v0.1beta_1000plus.fa ${dir}/scaffolds.fa
	ln -s /projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/ThirdStepCumulat/OutputMaker/LoweAED/Q903.all.maker.transcriptsLoweAED.fasta ${dir}/transcripts.fa
	ln -s /projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/ThirdStepCumulat/OutputMaker/LoweAED/Q903.all.polishedGenesLoweAED.gff ${dir}/genes.gff
	ln -s /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/uniprot-ncbi-filtered-labeled-nonSubsumed-9April2018.fa ${dir}/literature-AMPs.faa
else
	echo "Invalid genotype." 1>&2
fi
