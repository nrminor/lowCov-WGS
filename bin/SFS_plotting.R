#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# library(tidyverse)

# read in the site frequency spectrum file
SFS <- scan(args[1], skip =1, sep = " ")

# define the species
species <- as.character(args[2])

# produce an SFS bar plot for that species given available data, with and without singletons
filename = paste(species, "_", args[3], "_with_singletons", ".pdf", sep = "")
pdf(file = filename, 
    width = 8, height = 6)
barplot(SFS[1:length(SFS)], names.arg = 1:length(SFS), col = "blue",
        main = "Multi-Population Site Frequency Spectrum")
dev.off()

filename = paste(species, "_", args[3], "_without_singletons", ".pdf", sep = "")
pdf(file = filename, 
    width = 8, height = 6)
barplot(SFS[2:length(SFS)], names.arg = 2:length(SFS), col = "blue",
        main = "Multi-Population SFS (without singletons)")
dev.off()