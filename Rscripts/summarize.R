#!/usr/bin/env Rscript

# summarize transcript TPMs to genes

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
    genotype <- args[1]
    reads <- genotype
} else if (length(args) == 2) {
    genotype <- args[1]
    reads <- args[2]
} else {
    stop("USAGE: summarize.R <Transcriptome Genotype> <Reads Genotype>")
}

library(dplyr)
library(tximport)
library(rhdf5)

# load METADATA file
# At least two columns, 'sample' and 'treatment'

if (transcriptome == "PG29") {
    if (reads == "PG29") {
        metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples_tissue.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
        specific_dir<-"annotated_transcripts"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,"abundance.h5")
        tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
        
        num_treatments=length(unique(samples$treatment))
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        
        df_avg <- df
        
        colnames(df_avg) <- unique(samples$treatment)
        df_avg$gene <- rownames(df)
        rownames(df_avg) <- rownames(df)
        
        df_avg <- df_avg[defensins,c(num_treatments+1,1:num_treatments)]
        
        combined <- df_avg
        
    } else if (reads == "Q903") {
        metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples_wpw.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
        specific_dir<-"annotated_transcripts/replicates"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
        tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
        
        num_treatments=length(unique(samples$treatment))
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        h5closeAll()
        df_avg <- df %>% transmute(avg1=rowMeans(df[,1:num_reps]),
                                   avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]), 
                                   avg3=rowMeans(df[,((num_reps*2)+1):(num_reps*3)]))
        
        colnames(df_avg) <- unique(samples$treatment)
        df_avg$gene <- rownames(df)
        rownames(df_avg) <- rownames(df)
        
        df_avg <- df_avg[defensins,c(num_treatments+1,1:num_treatments)]
        
        metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples_stonecell.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
        specific_dir<-"annotated_transcripts/replicates"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
        tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
        
        num_treatments=length(unique(samples$treatment))
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
     
        df_avg2 <- df %>% transmute(avg1=rowMeans(df[,1:num_reps]),
                                   avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]))
        
        colnames(df_avg2) <- unique(samples$treatment)
        df_avg2$gene <- rownames(df)
        rownames(df_avg2) <- rownames(df)
        
        df_avg2 <- df_avg2[defensins,c(num_treatments+1,1:num_treatments)]
        
        combined <- merge(df_avg,df_avg2)
    } else {
        stop("Invalid reads.")
    }
} else if (transcriptome == "Q903") {
    if (reads == "Q903") {
        metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_wpw.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir<-"/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
        specific_dir<-"annotated_transcripts/replicates"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
        tx2gene_dict<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
        
        num_treatments=length(unique(samples$treatment))
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        h5closeAll()
        df_avg <- df %>% transmute(avg1=rowMeans(df[,1:num_reps]),
                                   avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]), 
                                   avg3=rowMeans(df[,((num_reps*2)+1):(num_reps*3)]))
        
        colnames(df_avg) <- unique(samples$treatment)
        df_avg$gene <- rownames(df)
        rownames(df_avg) <- rownames(df)
        
        df_avg <- df_avg[defensins,c(num_treatments+1,1:num_treatments)]
        
        metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples_stonecell.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
        specific_dir<-"annotated_transcripts/replicates"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
        tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
        
        num_treatments=length(unique(samples$treatment))
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        
        df_avg2 <- df %>% transmute(avg1=rowMeans(df[,1:num_reps]),
                                    avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]))
        
        colnames(df_avg2) <- unique(samples$treatment)
        df_avg2$gene <- rownames(df)
        rownames(df_avg2) <- rownames(df)
        
        df_avg2 <- df_avg2[defensins,c(num_treatments+1,1:num_treatments)]
        
        combined <- merge(df_avg,df_avg2)
        
    } else if (reads == "PG29") {
        metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_tissue.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
        specific_dir <- "annotated_transcripts"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,"abundance.h5")
        tx2gene_dict<- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
        
        num_treatments=length(unique(samples$treatment))
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        
        df_avg <- df
        
        colnames(df_avg) <- unique(samples$treatment)
        df_avg$gene <- rownames(df)
        rownames(df_avg) <- rownames(df)
        
        df_avg <- df_avg[defensins,c(num_treatments+1,1:num_treatments)]
        
        combined <- df_avg
        
    } else {
        
    }
} else if (transcriptome == "WS77111") {
    if (reads == "PG29") {
        metadata <- "/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/samples_tissue.csv"
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
        
        df_avg <- df
        
        colnames(df_avg) <- unique(samples$treatment)
        df_avg$gene <- rownames(df)
        rownames(df_avg) <- rownames(df)
        
        df_avg <- df_avg[defensins,c(num_treatments+1,1:num_treatments)]
        
        combined <- df_avg
        
    } else if (reads == "Q903") {
        metadata<-"/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/samples_wpw.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir<-"/projects/spruceup/scratch/pglauca/WS77111/annotation/amp/kallisto"
        specific_dir<-"annotated_transcripts/replicates"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
        tx2gene_dict<-"/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/deseq/tx2gene.csv"
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
        h5closeAll()
        df_avg <- df %>% transmute(avg1=rowMeans(df[,1:num_reps]),
                                   avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]), 
                                   avg3=rowMeans(df[,((num_reps*2)+1):(num_reps*3)]))
        
        colnames(df_avg) <- unique(samples$treatment)
        df_avg$gene <- rownames(df)
        rownames(df_avg) <- rownames(df)
        
        df_avg <- df_avg[defensins,c(num_treatments+1,1:num_treatments)]
        
        metadata<-"/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/samples_stonecell.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir<-"/projects/spruceup/scratch/pglauca/WS77111/annotation/amp/kallisto"
        specific_dir<-"annotated_transcripts/replicates"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
        tx2gene_dict<-"/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/deseq/tx2gene.csv"
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
        
        df_avg2 <- df %>% transmute(avg1=rowMeans(df[,1:num_reps]),
                                    avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]))
        
        colnames(df_avg2) <- unique(samples$treatment)
        df_avg2$gene <- rownames(df)
        rownames(df_avg2) <- rownames(df)
        
        df_avg2 <- df_avg2[defensins,c(num_treatments+1,1:num_treatments)]
        
        combined <- merge(df_avg,df_avg2)
    } else {
        stop("Invalid reads.")
    }
} else {
    stop("Invalid transcriptome.")
}

write.table(combined, file=file.path(kallisto_dir,"summarized.csv"), sep=",",row.names = FALSE, quote = FALSE)
