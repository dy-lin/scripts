#!/usr/bin/env Rscript

# plot defensin stacked bar plot

library(ggplot2)
library(dplyr)
library(ggrepel)

axisfont <- 16
allfont <- 18
labelfont <- 8

vjust_label <- 1.5
infile <- "/projects/btl/scratch/dlin/spruce/defensins/defensins.csv"
outfile <- "test.png"

AMP <- "Defensin"
title <- paste(AMP, "Genes in Spruce")
ylabel <- paste("Number of", AMP, "Genes")
xlabel <- "Spruce"
    
defensins <- read.csv(infile, header=TRUE, sep=",")
ggplot(defensins,
       aes(fill=Defensin,
           y=Count,
           x=Spruce)) + 
    geom_bar(stat="identity") + 
    ylim(0,13) + 
    ggtitle(title) + 
    geom_text(aes(label=Count, 
                  vjust=vjust_label),
              size=labelfont) + 
    ylab(ylabel) + 
    theme(text = element_text(size=allfont),
          legend.title = element_blank(), 
          axis.text = element_text(size=axisfont), 
          axis.title = element_text(size=axisfont), 
          legend.text = element_text(size=axisfont)) +
    xlab(xlabel)

ggsave(outfile, 
       width=16, 
       height=9, 
       dpi=300, 
       device="png")
