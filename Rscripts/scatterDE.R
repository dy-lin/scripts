#!/usr/bin/env Rscript

# plot the scatter plots

metadata<-"/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/samples_stonecell.csv"
kallisto_dir<-"/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
specific_dir<-"annotated_transcripts/replicates"
tx2gene_dict<-"/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
outfile<-"test.png"

defensins<-c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")

axisfont<-10
labelfont<-3
allfont<-12

xlabel<-"Cortical Parenchyma"
ylabel<-"Developing Stone Cell"
title<-paste(ylabel,"vs",xlabel)

legendtitle<-""

# xvar=""
# yvar=""


# fill=""
# color=""
# shape=""

# Default values
vjust_label<-0
hjust_label<-0.5

# number of replicates
num_reps<-3


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

# rownames(samples) <- colnames(txi.kallisto$counts)
# dds <- DESeqDataSetFromTximport(txi.kallisto,samples ~ treatment)

# dds <- DESeq(dds)

# comparisons <- resultNames(dds)

# for (i in comparisons) {
#	res <- results(dds, name=i)
	# plot
#}



df <- as.data.frame(txi.kallisto$abundance)

# Change column ranges and order as needed
# REMEMBER TO ADD/REMOVE VARIABLES HERE
df_avg <- df %>% mutate(avg1=rowMeans(df[,1:num_reps]),
                        avg2=rowMeans(df[,(num_reps+1):(num_reps*2)]), 
                      #  avg3=rowMeans(df[,((num_reps*2)+1):(num_reps*3)])
                        )

df_avg_drop <- df_avg

# drop 0 or any other desired values
df_avg_drop[df_avg_drop==0] <- NA

# add row names
rownames(df_avg_drop) <- rownames(df)

# find minimum and maximum values for each condition to decide ylim/xlim
# this will zoom into the defensins in the plot (but may leave out important information)

defensin_subset <- subset(df_avg_drop,is.element(row.names(df_avg_drop),defensins))

# Remember to change avg1/avg2 here depending on which condition ...
ymax <-  max(defensin_subset$avg2) + 0.25
ymin <- min(defensin_subset$avg2) - 0.25

xmax <- max(defensin_subset$avg1) + 0.25
xmin <- min(defensin_subset$avg1) - 0.25

# plot
# remember to change avg1/avg2 TWICE (2X GEOM_POINTS)

ggplot(df_avg_drop, aes(x=log10(avg1), y=log10(avg2))) + 
    geom_point(na.rm = TRUE, alpha=0.5, color="grey") + 
    geom_smooth(method=lm, alpha=1, fill="yellow") + 
    geom_point(data = defensin_subset, 
               aes(x=log10(avg1), y=log10(avg2)), 
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
          axis.text = element_text(size=axisfont)) + 
    ylim(ylimit) + 
    xlim(xlimit) + 
    ggtitle(title)

ggsave(outfile, 
       dpi=300, 
       device="png", 
       width=16, 
       height=9, 
       # units="cm"
       )

