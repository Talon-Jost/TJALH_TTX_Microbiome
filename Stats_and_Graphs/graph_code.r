#____Graph_Creation_Code____

#Packages
rm(list=ls())
packages <- c('dplyr', 'MASS', 'ggplot2', 'lme4', 'nlme', 'vegan', 'readr', 'plotrix', 'effects', 'remotes', 'MuMIn', 'ggtext', 'car', 'minpack.lm', 'mgcv', 'tidymv')

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mia")

library(dplyr)
library(mia)
library(MASS)
library(ggplot2)
library(lme4) 
library(nlme)
library(vegan)
library(readr)
library(plotrix)
library(effects)
library(remotes)
library(MuMIn)
library(ggtext)
library(car)
library(minpack.lm)
#remotes::install_github('ProcessMiner/nlcor')
library(nlcor)
#remotes::install_github('femiguez/nlraa')
library(nlraa)
library(mgcv)
library(tidymv)
#remotes::install_github("gavinsimpson/gratia")
library(gratia)
library(tidyverse)
library(phyloseq)
library(plyr)
#remotes::install_github('jbisanz/qiime2R')
library(qiime2R)

#qiime2R functions - helpful

#read_qza() - Function for reading artifacts (.qza).
#qza_to_phyloseq() - Imports multiple artifacts to produce a phyloseq object.
#read_q2metadata() - Reads qiime2 metadata file (containing q2-types definition line)
#write_q2manifest() - Writes a read manifest file to import data into qiime2
#theme_q2r() - A ggplot2 theme for for clean figures.
#print_provenance() - A function to display provenance information.
#is_q2metadata() - A function to check if a file is a qiime2 metadata file.
#parse_taxonomy() - A function to parse taxonomy strings and return a table where each column is a taxonomic class.
#parse_ordination() - A function to parse the internal ordination format.
#read_q2biom() - A function for reading QIIME2 biom files in format v2.1
#make_clr() - Transform feature table using centered log2 ratio.
#make_proportion() - Transform feature table to proportion (sum to 1).
#make_percent() - Transform feature to percent (sum to 100).
#interactive_table() - Create an interactive table in Rstudio viewer or rmarkdown html.
#summarize_taxa()- Create a list of tables with abundances sumed to each taxonomic level.
#taxa_barplot() - Create a stacked barplot using ggplot2.
#taxa_heatmap() - Create a heatmap of taxonomic abundances using gplot2.
#corner() - Show top corner of a large table-like obejct.
#min_nonzero() - Find the smallest non-zero, non-NA in a numeric vector.
#mean_sd() - Return mean and standard deviation for plotting.
#subsample_table() - Subsample a table with or without replacement.
#filter_features() - Remove low abundance features by number of counts and number of samples they appear in.

setwd('/home/tally/Desktop/final_qiime/Stats_and_Graphs/')
list.files()
all_removed_table <- read_qza('all_removed_table_101623.qza')
names(all_removed_table)
all_removed_table$uuid
all_removed_table$format
all_removed_table$contents
all_removed_table$data[1:5, 1:5]
all_removed_table$type
all_removed_table$version
all_removed_table$provenance

options(na.action(na.exclude))
metadata <- read_q2metadata('metadata_bc.tsv')
head(metadata)[1:5, 1:9]

#Phyloseq workflow trial
theme_set(theme_bw())
GP = GlobalPatterns

#Loading dataset
newt_newt_metadata <- read.csv('new_newt_metadata.csv', header = TRUE)
head(newt_newt_metadata, 4)
dim(newt_newt_metadata)
site = newt_newt_metadata$Site
setwd('/home/tally/Desktop/final_qiime/Stats_and_Graphs/Bray-curtis_Site_PERMANOVA')
list.files()
bc_permanova_site_overall <- read.csv('bray-curtis_PERMANOVA_raw_data.csv', header = TRUE)
head(bc_permanova_site_overall, 4)

bc_distance_matrix <- read_tsv('distance-matrix.tsv', col_names = TRUE)
bc_distance_matrix[1] = NULL
head(bc_distance_matrix)
imported_distance_matrix <- as.dist(bc_distance_matrix)
dim(bc_distance_matrix)
dist <- read_qza('bray_curtis_distance_matrix.qza')$data




site_adonis <- adonis2(bc_distance_matrix~Site, data = metadata, permutations = 999)
summary()

