#!/usr/bin/env Rscript

# plot the scatter plots

# load METADATA file
# At least two columns, 'sample' and 'treatment'

## PG29
# metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples.csv"
# samples <- read.csv(metadata, header = TRUE, sep=",")
# kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
# specific_dir<-"annotated_transcripts"
# files <- file.path(kallisto_dir,samples$treatment,specific_dir, samples$sample, "abundance.h5")
# tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
# defensins<-c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
# reads <- ""

## Q903
# metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_wpw.csv"
# metadata<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/samples_stonecell.csv"
# samples <- read.csv(metadata, header = TRUE, sep=",")
# kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
# specific_dir <- "annotated_transcripts/replicates"
# files <- file.path(kallisto_dir,samples$treatment,specific_dir, samples$sample, "abundance.h5")
# tx2gene_dict<- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
# defensins<-c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
# reads <- ""

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
               "DB47_00080436")
reads <- "(using Interior Spruce RNAseq reads)"

outfile<-"scatter.png"

axisfont <- 14
labelfont <- 10
allfont <- 18

xlabel<-"Bark"
ylabel<-"Young_Buds"
title<-paste(ylabel,"vs",xlabel)
subtitle <- reads

legendtitle<-""

# xvar=""
# yvar=""


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
num_reps <- nrow(samples)/num_treatments

names(files) <- samples$sample
tx2gene <- read.csv(tx2gene_dict)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

df <- as.data.frame(txi.kallisto$abundance)

# Change column ranges and order as needed
# REMEMBER TO ADD/REMOVE VARIABLES HERE
if (num_reps > 1) {
  df_avg <- df %>% mutate(avg1=rowMeans(df[,1:num_reps]),
                          avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]), 
                          avg3=rowMeans(df[,((num_reps*2)+1):(num_reps*3)]),
                          # avg4=rowMeans(df[,((num_reps*3)+1):(num_reps*4)]),
                          # avg5=rowMeans(df[,((num_reps*4)+1):(num_reps*5)]),
                          # avg6=rowMeans(df[,((num_reps*5)+1):(num_reps*6)]),
                          # avg7=rowMeans(df[,((num_reps*6)+1):(num_reps*7)]),
                          # avg8=rowMeans(df[,((num_reps*7)+1):(num_reps*8)])
                          )
} else {
  df_avg <- df
  colnames(df_avg) <- paste("avg",1:num_treatments, sep="")
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
ymax <- max(log10(defensin_subset$avg8), na.rm = TRUE) + 0.25
ymin <- min(log10(defensin_subset$avg8), na.rm = TRUE) - 0.25

xmax <- max(log10(defensin_subset$avg1), na.rm = TRUE) + 0.25
xmin <- min(log10(defensin_subset$avg1), na.rm = TRUE) - 0.25

min <- min(c(xmin,ymin))
max <- max(c(xmax,ymax))

ylimit <- c(min,max)
xlimit <- ylimit

# plot
# remember to change avg1/avg2 TWICE (2X GEOM_POINTS)
ggplot(df_avg_drop, aes(x=log10(avg1), y=log10(avg8))) + 
    geom_point(na.rm = TRUE, alpha=0.5, color="grey") + 
    geom_smooth(method=lm, alpha=1, fill="yellow") + 
    geom_point(data = defensin_subset, 
               aes(x=log10(avg1), y=log10(avg8)), 
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

ggsave(paste(ylabel,"vs",xlabel,outfile,sep="_"), 
       dpi=300, 
       device="png", 
       width=16, 
       height=9, 
       # units="cm"
         )

