#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# defining the relevant input parameters based on what was specified in
# nextflow.config, or in the command line of the workflow run
sfs <- as.integer(scan(args[1], skip =1, sep = " "))
species <- as.character(args[2])
pop <- as.character(args[3])
prep <- as.character(args[4])
sample_size <- as.numeric(args[5])
sfs <- sfs[1:round(sample_size/2)]
sfs <- paste(sfs, collapse = " ")
genome_length <- as.numeric(args[6])
year_per_generation <- as.numeric(args[7])
mutation_rate <- as.character(args[8])
random_seed <- as.numeric(args[9])
whether_folded <- as.logical(args[10])
# project_dir <- as.character(args[11])
# path_to_stairway_plot <- as.character(args[12])

# creating a dataframe that will be filled out into the Stairway plot blueprint
# file
blueprint <- data.frame("blueprint_file" = "# blueprint file")

# Adding population ID line
blueprint <- rbind(blueprint, paste("popid: ", species, "_", 
									prep, "_", pop, sep = ""))

# Adding sample size line
blueprint <- rbind(blueprint, paste("nseq: ", sample_size, sep = ""))

# Adding number of sites line
blueprint <- rbind(blueprint, paste("L: ", genome_length, sep = ""))

# Adding whether folded SFS line
blueprint <- rbind(blueprint, paste("whether_folded: ", whether_folded, sep = ""))

# Adding SFS line
blueprint <- rbind(blueprint, paste("SFS: ", sfs, sep = ""))

# Adding smallest SFS bin to be used
blueprint <- rbind(blueprint, 
				   paste("smallest_size_of_SFS_bin_used_for_estimation: ", 
						 "1", sep = ""))

# Adding largest SFS bin to be used
blueprint <- rbind(blueprint, 
				   paste("largest_size_of_SFS_bin_used_for_estimation: ", 
						 as.character(as.integer(sample_size/2)), sep = ""))

# Adding model training proportion
blueprint <- rbind(blueprint, paste("pct_training: ", 0.67, sep = ""))

# Adding breakpoints
blueprint <- rbind(blueprint,
				   paste("nrand:",
				   round(sample_size -(2*7)), round(sample_size - (2*5)),
				   round(sample_size - (2*3)), round(sample_size - 2),
						 sep = " "))

# Adding input directory
blueprint <- rbind(blueprint, paste("project_dir: ", getwd(), sep = ""))

# Adding Stairway plot software directory
blueprint <- rbind(blueprint, paste("stairway_plot_dir: ", 
				   getwd(), "/stairway_plot_es", sep = ""))

# Adding number of input files for estimation line
blueprint <- rbind(blueprint, paste("ninput: 200", sep = ""))

# Adding random seed line
blueprint <- rbind(blueprint, paste("random_seed: ", random_seed, sep = ""))

# Adding mutation rate line
blueprint <- rbind(blueprint, paste("mu: ", mutation_rate, sep = ""))

# Adding year(s) per generation line
blueprint <- rbind(blueprint, paste("year_per_generation: ", 
				   year_per_generation, sep = ""))

# Adding plot title line
blueprint <- rbind(blueprint, paste("plot_title: ", 
				   species, "_", prep, "_", pop, sep = ""))

# Adding default settings for x- and y-axis limits
blueprint <- rbind(blueprint, "xrange: 0,0")
blueprint <- rbind(blueprint, "yrange: 0,0")
blueprint <- rbind(blueprint, "xspacing: 2")
blueprint <- rbind(blueprint, "yspacing: 2")

# Adding fontsize (12 pt) line
blueprint <- rbind(blueprint, "fontsize: 12")

# naming Stairway plot blueprint file and exporting it
write.table(blueprint, 
			paste(species, "_", pop, "_", prep, ".blueprint", sep = ""),
			append = FALSE, quote = FALSE, sep = "\t",
			col.names = FALSE, row.names = FALSE)

# Command to construct stairway plot BASH script
# command <- paste("java -cp", path_to_stairway_plot, ,
# 				 sep = " ")
# system(command)