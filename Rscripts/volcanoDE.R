#!/usr/bin/env Rscript

# plot the volcano plots

# Variables to be modified based on conditions

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  transcriptome <- args[1]
  ref <- args[2]
} else {
  stop("USAGE: volcanoDE.R <Transcriptome Genotype> <Reads Treatment>")
}

# load METADATA file
# At least two columns, 'sample' and 'treatment'

# test

if (transcriptome == "PG29") {
  if (ref == "tissue") {
    metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples_tissue.csv"
    samples <- read.csv(metadata, header = TRUE, sep=",")
    kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
    specific_dir<-"annotated_transcripts"
    files <- file.path(kallisto_dir,samples$treatment,specific_dir,"abundance.h5")
    tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
    defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
    reads <- "(using PG29 RNAseq reads)"
  } else if (ref == "wpw") {
    metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples_wpw.csv"
    samples <- read.csv(metadata, header = TRUE, sep=",")
    kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
    specific_dir<-"annotated_transcripts/replicates"
    files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
    tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
    defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
    reads <- "(using Q903 RNAseq reads)"
  } else if (ref == "stonecell") {
    metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples_stonecell.csv"
    samples <- read.csv(metadata, header = TRUE, sep=",")
    kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
    specific_dir<-"annotated_transcripts/replicates"
    files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
    tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
    defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
    reads <- "(using Q903 RNAseq reads)"
  } else {
    stop("Invalid reads.")
  }
} else if (transcriptome == "Q903") {
  if (ref == "wpw") {
    metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_wpw.csv"
    samples <- read.csv(metadata, header = TRUE, sep=",")
    kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
    specific_dir <- "annotated_transcripts/replicates"
    files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
    tx2gene_dict<- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
    defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
    reads <- "(using Q903 RNAseq reads)"
    
  } else if (ref == "stonecell") {
    metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_stonecell.csv"
    samples <- read.csv(metadata, header = TRUE, sep=",")
    kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
    specific_dir <- "annotated_transcripts/replicates"
    files <- file.path(kallisto_dir,samples$treatment,specific_dir,samples$sample,"abundance.h5")
    tx2gene_dict<- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
    defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
    reads <- "(using Q903 RNAseq reads)"
  } else if (ref == "tissue") {
    metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_tissue.csv"
    samples <- read.csv(metadata, header = TRUE, sep=",")
    kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
    specific_dir <- "annotated_transcripts"
    files <- file.path(kallisto_dir,samples$treatment,specific_dir,"abundance.h5")
    tx2gene_dict<- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
    defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
    reads <- "(using PG29 RNAseq reads)"
  } else {
    stop("Invalid reads.")
  }
} else if (transcriptome == "WS77111") {
  if (ref == "tissue") {
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
  } else if (ref == "wpw") {
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
  } else if (ref == "stonecell") {
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
  } else {
    stop("Invalid reads.")
  }
} else {
  stop("Invalid transcriptome.")
}

outfile<-"volcano.png"

# Default font sizes
axisfont <- 14
labelfont <- 5
allfont <- 18

# padj threshold default, can be modified
padj_cutoff<-0.001

# lfc threshold default, can be modified
lfc_val <-2
lfc <- c(-lfc_val, lfc_val)

# Default for this type of plot
legendtitle<-"adjusted\np-value"

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
library(stringr)


names(files) <- samples$sample
tx2gene <- read.csv(tx2gene_dict)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

rownames(samples) <- colnames(txi.kallisto$counts)
dds <- DESeqDataSetFromTximport(txi.kallisto,samples, ~ treatment)

dds <- DESeq(dds)

comparisons <- resultsNames(dds)

for (i in comparisons[-1]) {
  res <- results(dds, name=i)
  # for each condition vs control, plot a volcano plot
  # convert to dataframe
  
  condition <- paste(toupper(substr(unlist(strsplit(str_remove(i,"treatment_"),"_vs_"))[1],1,1)),
                     substr(unlist(strsplit(str_remove(i,"treatment_"),"_vs_"))[1],2,nchar(unlist(strsplit(i,"_vs_"))[1])), 
                     sep = "")
  control <- paste(toupper(substr(unlist(strsplit(str_remove(i,"treatment_"),"_vs_"))[2],1,1)),
                   substr(unlist(strsplit(str_remove(i,"treatment_"),"_vs_"))[2],2,nchar(unlist(strsplit(i,"_vs_"))[1])), 
                   sep = "")
  vs <- "vs"
  
  full <- unlist(strsplit(condition,"_"))
  if (length(na.omit(full)) != 1) {
    word1 <- full[1]
    word2 <- full[2]
    word2 <- paste(toupper(substr(word2,1,1)),substr(word2,2,nchar(word2)),sep="")
    condition <- paste(word1, word2)
    
  } else {
    condition <- gsub("_", " ", condition)
  }
  
  full <- unlist(strsplit(control,"_"))
  if (length(na.omit(full)) != 1) {
    word1 <- full[1]
    word2 <- full[2]
    word2 <- paste(toupper(substr(word2,1,1)),substr(word2,2,nchar(word2)),sep="")
    control <- paste(word1, word2)
    
  } else {
    control <- gsub("_", " ", control)
  }
  
  title=paste(condition,vs,control,reads)
  
  condition <- unlist(strsplit(str_remove(i,"treatment_"),"_vs_"))[1]
  control <- unlist(strsplit(str_remove(i,"treatment_"),"_vs_"))[2]
  
  df <- as.data.frame(res)
  
  # add -log(pvalue) column, and padj category
  df_log <- mutate(df, logp=-log10(pvalue), padjcat=ifelse(padj < padj_cutoff, paste("padj <",padj_cutoff), paste("padj >=", padj_cutoff)))
  
  # reassign rownames
  rownames(df_log) <- rownames(df)
  
  # convert padjcat to factor
  df_log$padjcat <- as.factor(df_log$padjcat)
  defensin_subset <- subset(df_log,is.element(rownames(df_log),defensins))
  
  # find minimum and maximum values for each condition to decide ylim/xlim
  # this will zoom into the defensins in the plot (but may leave out important information)
  ymax <-  max(defensin_subset$logp) + 0.25
  ymin <- min(defensin_subset$logp) - 0.25
  
  xmax <- max(defensin_subset$log2FoldChange) + 0.25
  xmin <- min(defensin_subset$log2FoldChange) - 0.25
  
  ylimit <- c(0,ymax)
#  xlimit <- c(xmin, xmax)
#  xlimit <- c(-5,5)
  xlimit <- c(-1*xmax, xmax)
  
  # Plot discrete
  ggplot(df_log, aes(x=log2FoldChange, y=logp, color=padjcat)) + 
      geom_point(na.rm = TRUE,alpha=0.1) +
      geom_point(data=defensin_subset, 
                 alpha=1, 
                 color="green4")  +
      geom_text_repel(label=ifelse(is.element(rownames(df_log),defensins), rownames(df_log), ""),
                      size=labelfont, 
#                      vjust=0.5, 
 #                     hjust=0.5, 
                      color="green4") + 
      geom_vline(xintercept=lfc,linetype="dotted") + 
      theme(text = element_text(size=allfont),
            axis.text = element_text(size=axisfont)) + 
      labs(color=legendtitle) + 
      ggtitle(title) +
      xlim(xlimit) + 
      ylim(ylimit) +
      ylab("-log(p-value)") + 
      xlab("log2FoldChange")
  ggsave(paste("discrete",condition,control,outfile, sep="_"), 
         dpi=300, 
         device="png",
         width=16, 
         height=9, 
 #        units="cm"
  )
  
  # Plot continuous
  ggplot(df_log, 
         aes(x=log2FoldChange,
             y=logp, 
             color=padj)) + 
      geom_point(na.rm = TRUE, 
                 alpha=0.1) + 
      geom_point(data=defensin_subset,
                 alpha=1, 
                 color="red") +
      geom_text_repel(label=ifelse(is.element(rownames(df_log),defensins),rownames(df_log),""), 
                      size=labelfont, 
 #                     vjust=vjust_label, 
 #                     hjust=hjust_label, 
                      color="red") + 
      geom_vline(xintercept=lfc,
                 linetype="dotted") + 
      theme(text = element_text(size=allfont),
            axis.text = element_text(size=axisfont)) + 
      labs(color=legendtitle) + 
      ggtitle(title) +
      ylab("-log(p-value)") + 
      xlab("log2FoldChange") + 
      xlim(xlimit) + 
      ylim(ylimit)
  
  ggsave(paste("continuous",condition,control,outfile, sep="_"), 
         dpi=300, 
         device="png", 
         width=16, 
         height=9, 
        # units="cm"
         )
  # make one without a period
  padj_mod <- gsub("[.]","",as.character(padj_cutoff)) 
  df_log$gene <- rownames(df_log)
  write.table(filter(df_log,abs(log2FoldChange) > lfc_val & df_log$padj < padj_cutoff)$gene, file=file.path(kallisto_dir,condition,specific_dir,paste("genes.p",padj_mod,".lfc",lfc_val,".txt",sep="")), quote = FALSE,row.names =  FALSE,col.names = FALSE)
}
