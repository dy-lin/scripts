#!/usr/bin/env Rscript

# plot the volcano plots

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

# padj threshold default, can be modified
padj_cutoff<-0.05

# lfc threshold default, can be modified
lfc<-c(-1,1)

# number of replicates default, can be modified
num_reps<-4

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
  res <- results(dds, name=i)
  # for each condition vs control, plot a volcano plot
  # convert to dataframe

  condition <- paste(toupper(substr(strsplit(i,"_")[[1]][2],1,1)),
                     substr(strsplit(i,"_")[[1]][2],2,nchar(strsplit(i,"_")[[1]][1])), 
                     sep = "")
  control <- paste(toupper(substr(strsplit(i,"_")[[1]][4],1,1)),
                   substr(strsplit(i,"_")[[1]][4],2,nchar(strsplit(i,"_")[[1]][1])), 
                   sep = "")
  
  vs <- strsplit(i,"_")[[1]][3]
  title=paste(condition,vs,control)
  
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
  xlimit <- c(xmin, xmax)
  
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
                 color="green4") +
      geom_text_repel(label=ifelse(is.element(rownames(df_log),defensins),rownames(df_log),""), 
                      size=labelfont, 
 #                     vjust=vjust_label, 
 #                     hjust=hjust_label, 
                      color="green4") + 
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
}
