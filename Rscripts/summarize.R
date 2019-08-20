#!/usr/bin/env Rscript

# summarize transcript TPMs to genes

library(dplyr)

args = commandArgs(trailingOnly=TRUE)
infile <- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/control/annotated_transcripts/replicates/control.tsv"
# infile <- args[1]
# treatment <- args[2]
treatment <- "control"

df <- read.csv(infile, header = TRUE, sep="\t")

defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")

grouped <- group_by(df,target_id)
summarised <- summarise(grouped, sum1=sum(tpm), sum2=sum(tpm.1), sum3=sum(tpm.2), sum4=sum(tpm.3))

metadata<-"/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/samples_wpw.csv"
kallisto_dir<-"/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
specific_dir<-"annotated_transcripts/replicates"
tx2gene_dict<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
outfile<-"test.png"

# load METADATA file
# At least two columns, 'sample' and 'treatment'
samples <- read.csv(metadata, header = TRUE, sep=",")
num_treatments=length(unique(samples$treatment))
num_reps=nrow(samples)/num_treatments

files <- file.path(kallisto_dir,samples$treatment, specific_dir, samples$sample, "abundance.h5")
names(files) <- samples$sample
tx2gene <- read.csv(tx2gene_dict)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)


# load METADATA file
# At least two columns, 'sample' and 'treatment'
samples <- read.csv(metadata, header = TRUE, sep=",")


files <- file.path(kallisto_dir,samples$treatment, specific_dir, samples$sample, "abundance.h5")
names(files) <- samples$sample
tx2gene <- read.csv(tx2gene_dict)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
df <- as.data.frame(txi.kallisto$abundance)

if (num_reps > 1 && num_treatments == 2) {
    df_avg <- df %>% transmute(avg1=rowMeans(df[,1:num_reps]),
                            avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]))
}

if (num_reps > 1 && num_treatments == 3) {
    df_avg <- df %>% transmute(avg1=rowMeans(df[,1:num_reps]),
                            avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]), 
                            avg3=rowMeans(df[,((num_reps*2)+1):(num_reps*3)]))
}

if (num_reps ==1 ) {
    df_avg <- df
}
colnames(df_avg) <- unique(samples$treatment)
df_avg$gene <- rownames(df)
rownames(df_avg) <- rownames(df)
df_avg <- df_avg[defensins,c(4,1,2,3)]
write.table(df_avg, file=file.path(kallisto_dir,"summarized.csv"), sep=",",row.names = FALSE, quote = FALSE)







