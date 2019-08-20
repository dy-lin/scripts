#!/usr/bin/env Rscript

# plot the two p-value histograms

# Variables to be modified based on conditions
metadata<-"/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/samples_wpw.csv"
# metadata <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/samples_stonecell.csv"
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


# load packages
library(tximport)
library(dplyr)
library(rhdf5)
library(DESeq2)
library(stringr)

# load METADATA file
# At least two columns, 'sample' and 'treatment'
samples <- read.csv(metadata, header = TRUE, sep=",")


files <- file.path(kallisto_dir,samples$treatment, specific_dir, samples$sample, "abundance.h5")
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
         xlab = "P-values")
    padj <- hist(res$padj, 
                 breaks=20, 
                 col="grey", 
                 main = paste("Adjusted P-values of", title), 
                 xlab = "Adjusted P-values")
    dev.off()
}

