print('R script Loaded')

args <- commandArgs(trailingOnly = TRUE)

coords_file <- args[1]
nodes_file <- args[2]
phenotype_file <- args[3]
expression_file <- args[4]

coords <- read.csv(coords_file) 
nodes <- read.csv(nodes_file)
phenotype <- read.csv(phenotype_file)
expression <- read.csv(expression_file)
