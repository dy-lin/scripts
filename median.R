#!/usr/bin/env Rscript

# Given a double text file (each value on a newline), prints the median of that dataset
path <- commandArgs(trailingOnly=TRUE)
coverage <- scan(path, double(), quote = "")
print(median(coverage, na.rm=TRUE))
