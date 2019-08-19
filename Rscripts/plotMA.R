#!/usr/bin/env Rscript

# plot the MA plots

# Variables to be modified based on conditions
metadata<-"/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/samples_wpw.csv"
kallisto_dir<-"/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
specific_dir<-"annotated_transcripts/replicates"
tx2gene_dict<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
outfile<-"test.png"

# Modify based on genotype defensins
defensins<-c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")

# Default font sizes
axisfont<-10
labelfont<-3
allfont<-12

# number of replicates default, can be modified
num_reps<-4

padj<-0.05
# If any commented out variables are used, remember to add them to the ggplot() line

# fill=""
# color=""
# shape=""

# Default values for vjust and hjust (middle), can be modified 
vjust_label<-0.5
hjust_label<-0


# load packages
library(tximport)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(rhdf5)
library(DESeq2)

# load METADATA file
# At least two columns, 'sample' and 'treatment'
samples <- read.csv(metadata, header = TRUE, sep=",")


files <- file.path(kallisto_dir,samples$treatment, specific_dir, samples$sample, "abundance.h5")
names(files) <- samples$sample
tx2gene <- read.csv(tx2gene_dict)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

conditions <- unique(samples$treatment)

rownames(samples) <- colnames(txi.kallisto$counts)
dds <- DESeqDataSetFromTximport(txi.kallisto,samples, ~ treatment)

dds <- DESeq(dds)

comparisons <- resultsNames(dds)

for (i in comparisons[-1]) {
    res <- results(dds, name=i)
    condition <- paste(toupper(substr(strsplit(i,"_")[[1]][2],1,1)),
                       substr(strsplit(i,"_")[[1]][2],2,nchar(strsplit(i,"_")[[1]][1])), 
                       sep = "")
    control <- paste(toupper(substr(strsplit(i,"_")[[1]][4],1,1)),
                     substr(strsplit(i,"_")[[1]][4],2,nchar(strsplit(i,"_")[[1]][1])), 
                     sep = "")
    
    vs <- strsplit(i,"_")[[1]][3]
    title=paste("Treatment:",condition,vs,control)
    # subtitle=paste("\nred: padj <", padj)
    png(filename = paste(condition,vs,control, outfile,sep="_"), width=2560, height=1440, pointsize=18, units="px")
    par(mar=c(5, 5, 4, 2) + 0.1)
    plotMA(res, main=title, alpha=padj, cex.lab=2, cex.axis=2, cex.main=2, cex=1)
    legend("topright", c(paste("padj <",padj), paste("padj >=",padj)), col= c("red", "black"), pch=20, cex=1.5)
    dev.off()
                                                 
}
