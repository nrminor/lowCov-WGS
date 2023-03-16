#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

# Read in the PCA output file
pca <- read.table("pca_file.txt", header = TRUE, sep = "\t", row.names = 1)

# Create a data frame with the first two principal components
pca_df <- data.frame(PC1 = pca$PC1, PC2 = pca$PC2)

# Plot the first two principal components
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point()