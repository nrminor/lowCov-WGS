#!/usr/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(data.table)

# Path to the mosdepth output file
mosdepth_file <- "path/to/mosdepth.bed.gz"

# Read in the mosdepth output file
mosdepth <- fread(mosdepth_file, sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "depth"))

# Group the data by scaffold/chromosome
by_scaffold <- split(mosdepth, mosdepth$chrom)

# Generate a plot for each scaffold/chromosome
for (i in 1:length(by_scaffold)) {
  scaffold <- names(by_scaffold[i])
  data <- by_scaffold[[i]]
  p <- ggplot(data, aes(x=start, y=depth)) +
    geom_line() +
    labs(title=scaffold, x="Position", y="Depth") +
    theme_bw()
  ggsave(paste0(scaffold, ".png"), p, width=8, height=6, dpi=300)
}
