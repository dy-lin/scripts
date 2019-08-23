#!/usr/bin/env Rscript

# plot the two p-value histograms

# load packages
library(tximport)
library(dplyr)
library(rhdf5)
library(DESeq2)
library(stringr)

# load METADATA file
# At least two columns, 'sample' and 'treatment'

## PG29
# metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples.csv"
# samples <- read.csv(metadata, header = TRUE, sep=",")
# kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
# specific_dir<-"annotated_transcripts"
# files <- file.path(kallisto_dir,samples$treatment, specific_dir, "abundance.h5")
# tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
# defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
# reads <- ""

## Q903
# metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_wpw.csv"
# metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_stonecell.csv"
# samples <- read.csv(metadata, header = TRUE, sep=",")
# kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
# specific_dir <- "annotated_transcripts/replicates"
# files <- file.path(kallisto_dir,samples$treatment, specific_dir, samples$sample, "abundance.h5")
# tx2gene_dict<- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
# defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
# reads <- ""

## WS77111
metadata <- "/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/samples.csv"
samples <- read.csv(metadata, header = TRUE, sep=",")
kallisto_dir<-"/projects/spruceup/scratch/pglauca/WS77111/annotation/amp/kallisto"
specific_dir<-"annotated_transcripts"
files <- file.path(kallisto_dir,samples$treatment, specific_dir, "abundance.h5")
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

outfile<-"p-val-hist.png"

# Default font sizes
axisfont <- 14
labelfont <- 10
allfont <- 18

names(files) <- samples$sample
tx2gene <- read.csv(tx2gene_dict)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

rownames(samples) <- colnames(txi.kallisto$counts)
dds <- DESeqDataSetFromTximport(txi.kallisto,samples, ~ treatment)

dds <- DESeq(dds)

comparisons <- resultsNames(dds)

for (i in comparisons[-1]) {
    condition <- paste(toupper(substr(strsplit(str_remove(i,"treatment_"),"_vs_")[[1]][1],1,1)),
                       substr(strsplit(str_remove(i,"treatment_"),"_vs_")[[1]][1],2,nchar(strsplit(i,"_")[[1]][1])), 
                       sep = "")
    control <- paste(toupper(substr(strsplit(str_remove(i,"treatment_"),"_vs_")[[1]][2],1,1)),
                     substr(strsplit(str_remove(i,"treatment_"),"_vs_")[[1]][2],2,nchar(strsplit(i,"_")[[1]][1])), 
                     sep = "")
    
    vs <- "vs"
    title=paste(condition,vs,control)
    res <- results(dds, name=i)
    png(filename = paste(condition,vs,control, outfile,sep="_"), width=2560, height=1440, pointsize=18, units="px")
    par(mfrow = c(2, 1))
    pval <- hist(res$pvalue, 
         breaks=20, 
         col="grey", 
         main = paste("P-values of",title), 
         xlab = "P-values",
         sub = reads
    )
    padj <- hist(res$padj, 
                 breaks=20, 
                 col="grey", 
                 main = paste("Adjusted P-values of", title), 
                 xlab = "Adjusted P-values")
    dev.off()
}

