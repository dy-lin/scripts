#!/usr/bin/env Rscript

# plot the scatter plots

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  transcriptome <- args[1]
  ref <- args[2]
} else {
  stop("USAGE: scatterDE.R <Transcriptome Genotype> <Reads Treatment>")
}

# load METADATA file
# At least two columns, 'sample' and 'treatment'
# test
# transcriptome <- "WS77111"
# ref <- "wpw"
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

outfile<-"scatter.png"

axisfont <- 14
labelfont <- 5
allfont <- 18

xlabel<-"Bark"
ylabel<-"Young_Buds"
title<-paste(ylabel,"vs",xlabel)
subtitle <- reads


# fill=""
# color=""
# shape=""

# Default values
vjust_label<-0
hjust_label<-0.5



# load packages
library(tximport)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(rhdf5)
library(DESeq2)


num_treatments <- length(unique(samples$treatment))
treatments <- unique(samples$treatment)
num_reps <- nrow(samples)/num_treatments

names(files) <- samples$sample
tx2gene <- read.csv(tx2gene_dict)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

df <- as.data.frame(txi.kallisto$abundance)

# Change column ranges and order as needed
# REMEMBER TO ADD/REMOVE VARIABLES HERE
if (num_reps > 1) {
  df_avg <- df
  for (i in 1:num_treatments) {
    varname <- as.character(treatments[i])
    df_avg <- mutate(df_avg, !!varname := rowMeans(df_avg[,((num_reps*i) - (num_reps-1)):(num_reps*i)]) )
  }
  # df_avg <- df %>% mutate(avg1=rowMeans(df[,1:num_reps]),
  #                         avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]),
  #                         avg3=rowMeans(df[,((num_reps*2)+1):(num_reps*3)]),
  #                         # avg4=rowMeans(df[,((num_reps*3)+1):(num_reps*4)]),
  #                         # avg5=rowMeans(df[,((num_reps*4)+1):(num_reps*5)]),
  #                         # avg6=rowMeans(df[,((num_reps*5)+1):(num_reps*6)]),
  #                         # avg7=rowMeans(df[,((num_reps*6)+1):(num_reps*7)]),
  #                         # avg8=rowMeans(df[,((num_reps*7)+1):(num_reps*8)])
  #                         )
} else {
  df_avg <- df
  colnames(df_avg) <- treatments
}

df_avg_drop <- df_avg

# drop 0 or any other desired values
df_avg_drop[df_avg_drop==0] <- NA

# add row names
rownames(df_avg_drop) <- rownames(df)

# find minimum and maximum values for each condition to decide ylim/xlim
# this will zoom into the defensins in the plot (but may leave out important information)

defensin_subset <- subset(df_avg_drop,is.element(row.names(df_avg_drop),defensins))

# Remember to change avg1/avg2 here depending on which condition ...

# Do A LOOP

for (i in 2:num_treatments) {
  
  
  ymax <- max(log10(defensin_subset[[i]]), na.rm = TRUE) + 0.25
  ymin <- min(log10(defensin_subset[[i]]), na.rm = TRUE) - 0.25
  
  xmax <- max(log10(defensin_subset[[treatments[1]]]), na.rm = TRUE) + 0.25
  xmin <- min(log10(defensin_subset[[treatments[1]]]), na.rm = TRUE) - 0.25
  
  min <- min(c(xmin,ymin))
  max <- max(c(xmax,ymax))
  
  ylimit <- c(min,max)
  xlimit <- ylimit
  
  control <- as.character(treatments[1])
  var <- as.character(treatments[i])
  
  xlabel<- paste(toupper(substr(control,1,1)),substr(control,2,nchar(control)),sep="")
  
  ylabel<-paste(toupper(substr(var,1,1)),substr(var,2,nchar(var)),sep="")
  full <- unlist(strsplit(ylabel,"_"))
  if (length(na.omit(full)) != 1) {
    word1 <- full[1]
    word2 <- full[2]
    word2 <- paste(toupper(substr(word2,1,1)),substr(word2,2,nchar(word2)),sep="")
    ylabel <- paste(word1, word2)
    
  } else {
    ylabel <- gsub("_", " ", ylabel)
  }
  
  full <- unlist(strsplit(xlabel,"_"))
  if (length(na.omit(full)) != 1) {
    word1 <- full[1]
    word2 <- full[2]
    word2 <- paste(toupper(substr(word2,1,1)),substr(word2,2,nchar(word2)),sep="")
    xlabel <- paste(word1, word2)
    
  } else {
    xlabel <- gsub("_", " ", xlabel)
  }
  
  title<-paste(ylabel,"vs",xlabel)
  subtitle <- reads
  
  ggplot(df_avg_drop, aes(x=log10(!!ensym(control)), y=log10(!!ensym(var)))) + 
    geom_point(na.rm = TRUE, alpha=0.5, color="grey") + 
    geom_smooth(method=lm, alpha=1, fill="yellow") + 
    geom_point(data = defensin_subset, 
               aes(x=log10(!!ensym(control)), y=log10(!!ensym(var))), 
               color="red",alpha=1, 
               show.legend = TRUE) + 
    geom_text_repel(label=ifelse(is.element(row.names(df_avg_drop),defensins), rownames(df_avg_drop),""),
                    size=labelfont,
                    color="red",
                    #                  vjust=vjust_label, 
                    #                    hjust=hjust_label
    ) + 
    xlab(paste(xlabel, "(log(TPM))")) + 
    ylab(paste(ylabel, "(log(TPM))")) + 
    theme(text = element_text(size=allfont),
          axis.text = element_text(size=axisfont),
          axis.title = element_text(size=axisfont)) + 
    ylim(ylimit) + 
    xlim(xlimit) + 
    ggtitle(title, subtitle=reads)
  
  ggsave(paste(var,"vs",control,outfile,sep="_"), 
         dpi=300, 
         device="png", 
         width=16, 
         height=9, 
         # units="cm"
  )
}
