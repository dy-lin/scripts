#!/usr/bin/env Rscript

# summarize transcript TPMs to genes

library(dplyr)
library(UpSetR)
library(tximport)
library(rhdf5)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
    transcriptome <- args[1]
    ref <- args[2]
} else {
    stop("USAGE: upset.R <Transcriptome Genotype> <Reads Treatments>")
}
# transcriptome <- "WS77111"
# reads <- "tissue"
# load METADATA file
# At least two columns, 'sample' and 'treatment'

if (transcriptome == "PG29") {
    if (ref == "PG29") {
        metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples_tissue.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
        specific_dir<-"annotated_transcripts"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,"abundance.h5")
        tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
        reads <- "(using PG29 RNAseq reads)"
        
        outfile<-"upset.png"
        
        num_treatments=length(unique(samples$treatment))
        treatments <- unique(samples$treatment)
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        df_avg <- df
        if (num_reps > 1) {
            df_avg <- df
            for (i in 1:num_treatments) {
                varname <- as.character(treatments[i])
                df_avg <- mutate(df_avg, !!varname := rowMeans(df_avg[,((num_reps*i) - (num_reps-1)):(num_reps*i)]) )
            }
        } else {
            df_avg <- df
        }
        df_avg$gene <- rownames(df)
        df_avg <- df_avg[,c(ncol(df_avg),(ncol(df_avg)-num_treatments):(ncol(df_avg)-1))]
        
    } else if (ref == "Q903") {
        metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples_wpw.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
        specific_dir<-"annotated_transcripts/Q903_reads/replicates"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
        tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
        reads <- "(using Q903 RNAseq reads)"

        outfile<-"upset.png"
        
        num_treatments=length(unique(samples$treatment))
        treatments <- unique(samples$treatment)
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        df_avg <- df
        if (num_reps > 1) {
            df_avg <- df
            for (i in 1:num_treatments) {
                varname <- as.character(treatments[i])
                df_avg <- mutate(df_avg, !!varname := rowMeans(df_avg[,((num_reps*i) - (num_reps-1)):(num_reps*i)]) )
            }
        } else {
            df_avg <- df
        }
        df_avg$gene <- rownames(df)
        df_avg <- df_avg[,c(ncol(df_avg),(ncol(df_avg)-num_treatments):(ncol(df_avg)-1))]
        
        metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples_stonecell.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
        specific_dir<-"annotated_transcripts/Q903_reads/replicates"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
        tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
        reads <- "(using Q903 RNAseq reads)"
        
        outfile<-"upset.png"
        
        num_treatments=length(unique(samples$treatment))
        treatments <- unique(samples$treatment)
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        df_avg2 <- df
        if (num_reps > 1) {
            df_avg2 <- df
            for (i in 1:num_treatments) {
                varname <- as.character(treatments[i])
                df_avg2 <- mutate(df_avg2, !!varname := rowMeans(df_avg2[,((num_reps*i) - (num_reps-1)):(num_reps*i)]) )
            }
        } else {
            df_avg2 <- df
        }
        df_avg2$gene <- rownames(df)
        df_avg2 <- df_avg2[,c(ncol(df_avg2),(ncol(df_avg2)-num_treatments):(ncol(df_avg2)-1))]
        
        df_avg <- merge(df_avg,df_avg2)
        
    } else {
        stop("Invalid reads.")
    }
} else if (transcriptome == "Q903") {
    if (ref == "Q903") {
        metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_wpw.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
        specific_dir <- "annotated_transcripts/replicates"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
        tx2gene_dict<- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
        reads <- "(using Q903 RNAseq reads)"
        
        outfile<-"upset.png"
        
        num_treatments=length(unique(samples$treatment))
        treatments <- unique(samples$treatment)
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        df_avg <- df
        if (num_reps > 1) {
            df_avg <- df
            for (i in 1:num_treatments) {
                varname <- as.character(treatments[i])
                df_avg <- mutate(df_avg, !!varname := rowMeans(df_avg[,((num_reps*i) - (num_reps-1)):(num_reps*i)]) )
            }
        } else {
            df_avg <- df
        }
        df_avg$gene <- rownames(df)
        df_avg <- df_avg[,c(ncol(df_avg),(ncol(df_avg)-num_treatments):(ncol(df_avg)-1))]
        
        metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_stonecell.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
        specific_dir <- "annotated_transcripts/replicates"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
        tx2gene_dict<- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
        reads <- "(using Q903 RNAseq reads)"
        
        outfile<-"upset.png"
        
        num_treatments=length(unique(samples$treatment))
        treatments <- unique(samples$treatment)
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        df_avg2 <- df
        if (num_reps > 1) {
            df_avg2 <- df
            for (i in 1:num_treatments) {
                varname <- as.character(treatments[i])
                df_avg2 <- mutate(df_avg2, !!varname := rowMeans(df_avg2[,((num_reps*i) - (num_reps-1)):(num_reps*i)]) )
            }
        } else {
            df_avg2 <- df
        }
        df_avg2$gene <- rownames(df)
        df_avg2 <- df_avg2[,c(ncol(df_avg2),(ncol(df_avg2)-num_treatments):(ncol(df_avg2)-1))]
        
        df_avg <- merge(df_avg,df_avg2)
        
    } else if (ref == "PG29") {
        metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_tissue.csv"
        samples <- read.csv(metadata, header = TRUE, sep=",")
        kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
        specific_dir <- "annotated_transcripts"
        files <- file.path(kallisto_dir,samples$treatment,specific_dir,"abundance.h5")
        tx2gene_dict<- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
        defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
        reads <- "(using PG29 RNAseq reads)"
        
        outfile<-"upset.png"
        
        num_treatments=length(unique(samples$treatment))
        treatments <- unique(samples$treatment)
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        df_avg <- df
        if (num_reps > 1) {
            df_avg <- df
            for (i in 1:num_treatments) {
                varname <- as.character(treatments[i])
                df_avg <- mutate(df_avg, !!varname := rowMeans(df_avg[,((num_reps*i) - (num_reps-1)):(num_reps*i)]) )
            }
        } else {
            df_avg <- df
        }
        df_avg$gene <- rownames(df)
        df_avg <- df_avg[,c(ncol(df_avg),(ncol(df_avg)-num_treatments):(ncol(df_avg)-1))]
        
    } else {
        stop("Invalid reads.")
    }
} else if (transcriptome == "WS77111") {
    if (ref == "PG29") {
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
        reads <- "(using PG29 RNAseq reads)"
        
        outfile<-"upset.png"
        
        num_treatments=length(unique(samples$treatment))
        treatments <- unique(samples$treatment)
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        df_avg <- df
        if (num_reps > 1) {
            df_avg <- df
            for (i in 1:num_treatments) {
                varname <- as.character(treatments[i])
                df_avg <- mutate(df_avg, !!varname := rowMeans(df_avg[,((num_reps*i) - (num_reps-1)):(num_reps*i)]) )
            }
        } else {
            df_avg <- df
        }
        df_avg$gene <- rownames(df)
        df_avg <- df_avg[,c(ncol(df_avg),(ncol(df_avg)-num_treatments):(ncol(df_avg)-1))]
        
    } else if (ref == "Q903") {
        metadata <- "/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/samples_wpw.csv"
        samples <- read.csv(metadata,header = TRUE, sep = ",")
        kallisto_dir <- "/projects/spruceup/scratch/pglauca/WS77111/annotation/amp/kallisto"
        specific_dir <- "annotated_transcripts/replicates"
        files <- file.path(kallisto_dir, samples$treatment,specific_dir, samples$sample, "abundance.h5")
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
        reads <- "(using Q903 RNAseq reads)"
        
        outfile<-"upset.png"
        
        num_treatments=length(unique(samples$treatment))
        treatments <- unique(samples$treatment)
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        df_avg <- df
        if (num_reps > 1) {
            df_avg <- df
            for (i in 1:num_treatments) {
                varname <- as.character(treatments[i])
                df_avg <- mutate(df_avg, !!varname := rowMeans(df_avg[,((num_reps*i) - (num_reps-1)):(num_reps*i)]) )
            }
        } else {
            df_avg <- df
        }
        df_avg$gene <- rownames(df)
        df_avg <- df_avg[,c(ncol(df_avg),(ncol(df_avg)-num_treatments):(ncol(df_avg)-1))]
        
        
        metadata <- "/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/samples_stonecell.csv"
        samples <- read.csv(metadata,header = TRUE, sep = ",")
        kallisto_dir <- "/projects/spruceup/scratch/pglauca/WS77111/annotation/amp/kallisto"
        specific_dir <- "annotated_transcripts/replicates"
        files <- file.path(kallisto_dir, samples$treatment,specific_dir, samples$sample, "abundance.h5")
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
        reads <- "(using Q903 RNAseq reads)"
        
        outfile<-"upset.png"
        
        num_treatments=length(unique(samples$treatment))
        treatments <- unique(samples$treatment)
        num_reps=nrow(samples)/num_treatments
        
        names(files) <- samples$sample
        tx2gene <- read.csv(tx2gene_dict)
        txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
        
        df <- as.data.frame(txi.kallisto$abundance)
        df_avg2 <- df
        if (num_reps > 1) {
            df_avg2 <- df
            for (i in 1:num_treatments) {
                varname <- as.character(treatments[i])
                df_avg2 <- mutate(df_avg2, !!varname := rowMeans(df_avg2[,((num_reps*i) - (num_reps-1)):(num_reps*i)]) )
            }
        } else {
            df_avg2 <- df
        }
        df_avg2$gene <- rownames(df)
        df_avg2 <- df_avg2[,c(ncol(df_avg2),(ncol(df_avg2)-num_treatments):(ncol(df_avg2)-1))]
        
        df_avg <- merge(df_avg,df_avg2)
        
    } else {
        stop("Invalid reads.")
    }
} else {
    stop("Invalid transcriptome.")
}

# flipped <- as.data.frame(df_avg[defensins,] %>% t)[-1,]
# flipped_upset <- flipped
flipped_upset <- filter(df_avg, is.element(gene,defensins))

flipped_upset$gene <- as.factor(flipped_upset$gene)
flipped_upset[,-1][flipped_upset[,-1] >= 1] <-1
flipped_upset[,-1][flipped_upset[,-1] < 1] <-0

png(filename = outfile, width=2560, height=1440, pointsize=18, units="px")
upset(flipped_upset, nsets=ncol(flipped_upset)-1, text.scale=4, point.size = 8, line.size=2)
dev.off()
