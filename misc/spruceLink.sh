#!/bin/bash
spruce=$1
if [[ "$#" -ne 1 && "$#" -ne 2 && "$#" -ne 3 ]]
then
	echo "USAGE: $(basename $0) <spruce genotype> [AMP Class] [specified directory] " 1>&2
	echo "DESCRIPTION: Creates softlinks for all spruce files in the specified directories." 1>&2
	exit 1
fi

## Cases: #args = 1, no AMP class specified, no directory specified
## #args=2, genotype specified, AMP class specified
## #args=2, genotype specified, directory specified
if [[ "$#" -eq 2 ]]
then
	if [[ ! -d $2 ]]
	then
		dir=$(pwd)
		case $2 in
			thionin*) ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/ncbi-spermatophyta-thionins-02Apr2019.faa ${dir}/thionins.faa;;
			hairpinin*) ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/ncbi-all-hairpinins-02Apr2019.faa ${dir}/hairpinins.faa;;
			hevein*) ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/ncbi-picea-hevein-18Mar2019.faa ${dir}/heveins.faa;;
			lipid*) ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/ncbi-picea-lipidtransfer-18Mar2019.faa ${dir}/lipid_transfer.faa;;
			thaumatin*) ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/ncbi-picea-thaumatins-18Mar2019.faa ${dir}/thaumatins.faa;;
			knottin*) ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/ncbi-spermatophyta-knottins-02Apr2019.faa ${dir}/knottins.faa;;
			snakin*) ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/ncbi-spermatophyta-snakins-02Apr2019.faa ${dir}/snakins.faa;;
			defensin*) ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/psitchensis-defensins.fa ${dir}/defensins.faa;;
			\?) echo "Invalid AMP." 1>&2; exit 1;;
		esac
	else
		dir=$(echo $2 | sed 's/\/$//')
	fi
elif [[ "$#" -eq 3 ]]
then
	dir=$(echo $3 | sed 's/\/$//')
fi

if [[ "$spruce" == "WS77111" ]]
then
	ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/uniprot-ncbi-filtered-labeled-nonSubsumed-9April2018.fa ${dir}/literature-AMPs.faa
	ln -fs /projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.all.polished.genesLoweAED.function_domain.gff ${dir}/genes.gff
	ln -fs /projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.all.maker.proteinsLoweAED.rename.function.fasta ${dir}/proteins.faa
	ln -fs /projects/spruceup_scratch/pglauca/WS77111/annotation/genome-annotation/Maker/Polish/FourthIterationCumulatNonntEdit/NcbiSubmission/WS77111.all.maker.transcriptsLoweAED.rename.function.fasta ${dir}/transcripts.fa
	ln -fs /projects/spruceup/pglauca/WS77111/assemblies/releases/version2/WS77111v2_release/WS77111-v2_1000plus_LGs.fa ${dir}/scaffolds.fa
elif [[ "$spruce" == "PG29" ]]
then
	ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/uniprot-ncbi-filtered-labeled-nonSubsumed-9April2018.fa ${dir}/literature-AMPs.faa
	ln -fs /projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29v5.all.maker.transcriptsLoweAED.renamed.function.fasta ${dir}/transcripts.fa
	ln -fs /projects/spruceup/interior_spruce/PG29/assemblies/releases/PG29-v5/PG29-v5_1000plus.fa ${dir}/scaffolds.fa
	ln -fs  /projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29v5.all.maker.proteinsLoweAED.renamed.function.fasta ${dir}/proteins.faa
	ln -fs /projects/spruceup_scratch/interior_spruce/PG29/annotation/genome-annotation/PG29v5/Maker/Polish/ThirdIterationCumulatNonntEdit/NcbiSubmission/PG29.all.polished.genesLoweAED.function_domain.gff  ${dir}/genes.gff
elif [[ "$spruce" == "Q903" ]]
then
	ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.all.maker.proteinsLoweAED.renamed.function.fasta ${dir}/proteins.faa
	ln -fs /projects/spruceup/psitchensis/Q903/assembly/releases/version1/Q903_v1_1000plus.fa ${dir}/scaffolds.fa
	ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.all.maker.transcriptsLoweAED.renamed.function.fasta ${dir}/transcripts.fa
	ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/genome-annotation/Maker/Polish/ThirdIteration/NcbiSubmission/Q903.all.polished.genesLoweAED.function_domain.gff ${dir}/genes.gff
	ln -fs /projects/spruceup_scratch/psitchensis/Q903/annotation/amp/sequences/uniprot-ncbi-filtered-labeled-nonSubsumed-9April2018.fa ${dir}/literature-AMPs.faa
else
	echo "Invalid genotype." 1>&2
fi

