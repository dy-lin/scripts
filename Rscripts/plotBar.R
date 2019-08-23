#!/usr/bin/env Rscript

# plot distribution of p-values (reg, and adjusted)

## Q903
# infile <- "/projects/spruceup_scratch/psitchensis/Q903/annotation/amp/kallisto/final_summary.csv"
# title <- "Defensin Expression in Sitka Spruce"

## PG29
# infile <- "/projects/spruceup_scratch/interior_spruce/PG29/annotation/amp/kallisto/final_summary.csv"
# title <- "Defensin Expression in Interior Spruce"

## WS77111
infile <- "/projects/spruceup_scratch/pglauca/WS77111/annotation/amp/kallisto/final_summary.csv"
title <- "Defensin Expression in White Spruce (using Interior Spruce RNAseq reads)"

outfile <- "bar.png"

library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(grid)


axisfont <- 14
labelfont <- 10
allfont <- 18

summary <- read.csv(infile, header = T, sep=",")

TPM_cutoff <- 1

# Round all values < 1 TPM to 0, and round to one decimal place

summary <- summary %>% mutate(TPM_ge1 = ifelse(TPM>TPM_cutoff,TPM,0)) %>% mutate_at(vars(TPM_ge1), funs(round(., 1)))

num <- length(unique(summary$Condition))
# if there are more than one category of conditions, use gridExtra to combine them.
if (num > 1 ) {
    
    
    ylabel <- "Transcripts per Million (TPM)"
    subtitle <- paste("TPM >=", TPM_cutoff)
    
    ymax <- max(summary$TPM_ge1) + 5
    ylimit <- c(0,ymax)
    
    angle <- 0
    # vjust 
    plots <- list()
    count <- 1
    for (j in rev(unique(summary$Condition))) {
        subset <- subset(summary,Condition == j)
        
        xlabel <- paste("Treatment:", j)
        
        if (count > 1) {
            plots[[count]] <- ggplot(data=subset, 
                   aes(fill=Gene, 
                       x=Treatment, 
                       y=TPM_ge1)) + 
                geom_bar(width=0.5 ,
                         position=position_dodge(0.5),
                         stat="identity") + 
                xlab(xlabel) + 
                ylab(NULL) + 
                theme(text = element_text(size=allfont),
                      axis.text = element_text(angle = angle, 
                   #                            hjust = hjust, 
                                               size=axisfont
                      ), axis.title = element_text(size = axisfont)
                )  +
                ylim(ylimit) +
                ggtitle("", subtitle= "")
        } else {
            plots[[count]] <- ggplot(data=subset, 
                        aes(fill=Gene, 
                            x=Treatment,
                            y=TPM_ge1)) + 
                geom_bar(width=0.5 ,
                         position=position_dodge(0.5),
                         stat="identity", show.legend = FALSE) + 
                xlab(xlabel) + 
                ylab(ylabel) + 
                ggtitle(title, 
                        subtitle = subtitle)  + 
                theme(text = element_text(size=allfont),
                      axis.text = element_text(angle = angle, 
                         #                      hjust = hjust, 
                                               size=axisfont
                      ), axis.title = element_text(size = axisfont)
                ) +
                ylim(ylimit)
        }
        count <- count + 1
    }
    combined <- do.call(grid.arrange,c(plots,ncol=2))
    ggsave(plot = combined, outfile, 
           width=16, 
           height=9, 
           dpi=300, 
           device="png")
} else {
    
    xlabel <- "Treatments"
    ylabel <- "Transcripts per Million (TPM)"
    subtitle <- paste("TPM >=", TPM_cutoff, reads)
    
    angle <- 20
    hjust <- 1
    # vjust 
    
    ggplot(data=summary, 
           aes(fill=Gene, 
               x=Treatment, 
               y=TPM_ge1)) + 
        geom_bar(width=0.5 ,
                 position=position_dodge(0.5),
                 stat="identity") + 
        xlab(xlabel) + 
        ylab(ylabel) + 
        ggtitle(title, 
                subtitle = subtitle)  + 
        theme(text = element_text(size=allfont),
              axis.text = element_text(angle = angle, 
                                         hjust = hjust, 
                                         size=axisfont
              ), axis.title = element_text(size = axisfont)
        ) 
   #     geom_text_repel(position = position_dodge(0.5), 
   #                     label=ifelse(summary$TPM_ge1>1, summary$TPM_ge1, ""), 
  #                      size=labelfont)
    
    ggsave(outfile, 
           width=16, 
           height=9, 
           dpi=300, 
           device="png")
    
}