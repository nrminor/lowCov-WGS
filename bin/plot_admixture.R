#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

file_prefix = args[1]
K_value = as.numeric(args[2])

# Load data
data <- read.table("output.q")

# Create plot
ggplot(data, aes(x=V1, y=V2, fill=as.factor(V3))) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("Sample ID") +
  ylab("Ancestry Proportions") +
  ggtitle("ADMIXTURE results") +
  scale_fill_brewer(palette="Dark2")