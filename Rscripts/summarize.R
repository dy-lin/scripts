#!/usr/bin/env Rscript

# summarize transcript TPMs to genes

library(dplyr)
library(tximport)
library(rhdf5)

# load METADATA file
# At least two columns, 'sample' and 'treatment'

## PG29
# metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples.csv"
# samples <- read.csv(metadata, header = TRUE, sep=",")
# kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
# specific_dir<-"annotated_transcripts"
# files <- file.path(kallisto_dir,samples$treatment,specific_dir,"abundance.h5")
# tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
# defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
# files <- file.path(kallisto_dir,samples$treatment,specific_dir,"abundance.h5")

## Q903
# metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_wpw.csv"
# metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_stonecell.csv"
# samples <- read.csv(metadata, header = TRUE, sep=",")
# kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
# specific_dir <- "annotated_transcripts/replicates"
# files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
# tx2gene_dict<- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
# defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
# files <- file.path(kallisto_dir,samples$treatment, specific_dir, samples$sample, "abundance.h5")

## WS77111
metadata <- "/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/samples.csv"
samples <- read.csv(metadata, header = TRUE, sep=",")
kallisto_dir<-"/projects/spruceup/scratch/pglauca/WS77111/annotation/amp/kallisto"
specific_dir<-"annotated_transcripts"
files <- file.path(kallisto_dir,samples$treatment,specific_dir,"abundance.h5")
tx2gene_dict <- "/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/deseq/tx2gene.csv"
defensins <- c("DB47_00018419",
               "DB47_00018420",
               "DB47_00018421",
               "DB47_00018422",
               "DB47_00018423",
               "DB47_00018424",
               "DB47_00028544",
               "DB47_00044066",
               "DB47_00073581",
               "DB47_00073614",
               "DB47_00080438")

num_treatments=length(unique(samples$treatment))
num_reps=nrow(samples)/num_treatments

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

df_avg <- df_avg[defensins,c(num_treatments+1,1:num_treatments)]

if(num_treatments != 8) {
    metadata<-"/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/samples_stonecell.csv"
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
        df_avg2 <- df %>% transmute(avg1=rowMeans(df[,1:num_reps]),
                                   avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]))
    }
    
    if (num_reps > 1 && num_treatments == 3) {
        df_avg2 <- df %>% transmute(avg1=rowMeans(df[,1:num_reps]),
                                   avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]), 
                                   avg3=rowMeans(df[,((num_reps*2)+1):(num_reps*3)]))
    }
    
    if (num_reps ==1 ) {
        df_avg2 <- df
    }
    
    colnames(df_avg2) <- unique(samples$treatment)
    df_avg2$gene <- rownames(df)
    rownames(df_avg2) <- rownames(df)
    df_avg2 <- df_avg2[defensins,c(3,1,2)]
    
    combined <-merge(df_avg,df_avg2)
    
} else {
    combined <- df_avg
}

write.table(combined, file=file.path(kallisto_dir,"summarized.csv"), sep=",",row.names = FALSE, quote = FALSE)