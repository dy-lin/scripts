#!/usr/bin/env Rscript

# Given a double text file (each value on a newline), prints the median of that dataset
path <- commandArgs(trailingOnly=TRUE)
coverage <- scan(path, double(), quote = "")
print(median(coverage, na.rm=TRUE))


## Usin AWK:
# sort -n file | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'
