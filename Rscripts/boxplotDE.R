#!/usr/bin/env Rscript

# plot the count boxplots

# load packages
library(tximport)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(rhdf5)
library(DESeq2)

# load METADATA file
# At least two columns, 'sample' and 'treatment'

# Variables to be modified based on conditions

## Q903
metadata<- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/samples_wpw.csv"
# metadata<- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto/samples_stonecell.csv"
samples <- read.csv(metadata, header = TRUE, sep=",")
kallisto_dir <- "/projects/spruceup/scratch/psitchensis/Q903/annotation/amp/kallisto"
specific_dir <- "annotated_transcripts/replicates"
files <- file.path(kallisto_dir,samples$treatment, specific_dir, samples$sample, "abundance.h5")
tx2gene_dict <- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/deseq/tx2gene.csv"
defensins <- c("E0M31_00027086", "E0M31_00027087", "E0M31_00055415", "E0M31_00093276", "E0M31_00093277")
reads <- ""

## PG29
# metadata<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/samples.csv"
# samples <- read.csv(metadata, header = TRUE, sep=",")
# kallisto_dir<-"/projects/spruceup/scratch/interior_spruce/PG29/annotation/amp/kallisto"
# specific_dir<-"annotated_transcripts"
# files <- file.path(kallisto_dir,samples$treatment, specific_dir, "abundance.h5")
# tx2gene_dict<-"/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/deseq/tx2gene.csv"
# defensins <- c("ABT39_00024884", "ABT39_00024885", "ABT39_00024887", "ABT39_00102286", "ABT39_00108568", "ABT39_00122613")
# reads <- ""

## WS77111
# metadata <- "/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/samples.csv"
# samples <- read.csv(metadata, header = TRUE, sep=",")
# kallisto_dir<-"/projects/spruceup/scratch/pglauca/WS77111/annotation/amp/kallisto"
# specific_dir<-"annotated_transcripts"
# files <- file.path(kallisto_dir,samples$treatment, specific_dir, "abundance.h5")
# tx2gene_dict <- "/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/deseq/tx2gene.csv"
# defensins <- c("DB47_00018419",
#                "DB47_00018420",
#                "DB47_00018421",
#                "DB47_00018422",
#                "DB47_00018423",
#                "DB47_00018424",
#                "DB47_00028544",
#                "DB47_00044066",
#                "DB47_00073581",
#                "DB47_00073614",
#                "DB47_00080438")
# reads <- "(using PG29 RNAseq reads)"

outfile<- "boxplot.png"

# Default font sizes
axisfont <- 14
labelfont <- 10
allfont <- 18

# Modify based on conditions e.g. Cells or Treatment
xlabel<-"Treatments"

# Default for this plot
ylabel<-"Normalized Counts"
# title<-""

# number of replicates default, can be modified
num_reps<-4

# legendtitle=""

# If any commented out variables are used, remember to add them to the ggplot() line

# ylimit=c()

# fill=""
# color=""
# shape=""


names(files) <- samples$sample
tx2gene <- read.csv(tx2gene_dict)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

conditions <- unique(samples$treatment)

rownames(samples) <- colnames(txi.kallisto$counts)
dds <- DESeqDataSetFromTximport(txi.kallisto,samples, ~ treatment)

dds <- DESeq(dds)

esf <- estimateSizeFactors(dds)
counts <- as.data.frame(counts(esf, normalized=TRUE))

gene_names <- NULL
for (i in defensins) {
    temp <- rep(i,num_reps*length(conditions))
    gene_names <- c(gene_names,temp)
}

rm(temp)
conds <- NULL
for (j in conditions) {
    temp <- rep(j,num_reps)
    conds <- c(conds, temp)
}
rm(temp)
gene_df <- as.data.frame(gene_names)
gene_df$treatment <- as.factor(rep(conds,length(defensins)))

normalized_counts <- NULL
for (l in defensins) {
    gene_temp <- as.data.frame(counts[l,] %>% t)
    normalized_counts <- c(normalized_counts,gene_temp[,1])
}

gene_df$counts <- normalized_counts

colnames(gene_df) <- c("Gene", "Treatment", "Normalized_Counts")

ggplot(gene_df, aes(x=Treatment, y=Normalized_Counts)) + 
    geom_boxplot() + 
    ggtitle("Normalized Counts for Defensin Genes",subtitle=reads) + 
    ylab(ylabel) + 
    xlab(xlabel) + 
    theme(text = element_text(size=allfont),
          axis.text = element_text(size=axisfont)) +
    facet_wrap(~ Gene, ncol=3)

ggsave(outfile, 
       device="png", 
       dpi=300, 
       width=16, 
       height=9, 
       #         units="cm"
)
# For individual plots
# for (i in defensins) {
#   gene_df <- as.data.frame(counts[i,] %>% t)
#   v1 <- NULL
#   for (j in conditions) {
#     temp <- rep(j,num_reps)
#     v1 <- c(v1, temp)
#   }
#   gene_df$treatment <- as.factor(v1)
#   colnames(gene_df) <- c("normalized_counts","treatment")
#   ggplot(gene_df, aes(x=treatment, y=normalized_counts)) + 
#       geom_boxplot() + 
#       ggtitle(i) + 
#       ylab(ylabel) + 
#       xlab(xlabel) + 
#       theme(text = element_text(size=allfont),
#             axis.text = element_text(size=axisfont))
#   ggsave(paste(i,outfile,sep="_"), 
#          device="png", 
#          dpi=300, 
#          width=16, 
#          height=9, 
# #         units="cm"
#          )
# }
