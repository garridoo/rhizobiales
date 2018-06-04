
# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("../plotting_functions.R")

# load plotting functions

library("ggplot2")
library("scales")
library("grid")
library("MASS")

# directories

results.dir <- "../../meta_16S/results/"
figures.dir <- "../../meta_16S/figures/"

# files

design.file <- paste(results.dir, "design.txt", sep="")
mapping.file <- paste(results.dir, "mapping.txt", sep="")
otu_table.file <- paste(results.dir, "otu_table.txt", sep="")
otu_table_rhizobia.file <- paste(results.dir, "otu_table_rhizobia_97.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
mapping <- read.table(mapping.file, header=T, sep="\t")

# read de-novo OTU table to get normalizing factors

otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
rownames(otu_table) <- otu_table[, 1]
otu_table <- otu_table[, -1]

# read reference-based OTU table of rhizobial sequences

otu_table_rhizobia <- read.table(otu_table_rhizobia.file, sep="\t", header=T, check.names=F)
rownames(otu_table_rhizobia) <- otu_table_rhizobia[, 1]
otu_table_rhizobia <- otu_table_rhizobia[, -1]

# match OTU tables

idx <- match(colnames(otu_table_rhizobia), colnames(otu_table))
otu_table <- otu_table[, idx]

# remove shallow samples

min_depth <- 500
norm_factors <- colSums(otu_table)
idx <- norm_factors >= min_depth
norm_factors <- norm_factors[idx]
otu_table <- otu_table[, idx]
otu_table_rhizobia <- otu_table_rhizobia[, idx]

# normalize rhizobial OTU table w.r.t. all taxa

otu_table_norm_all <- otu_table_rhizobia
for (i in 1:ncol(otu_table)) {otu_table_norm_all[, i] <- otu_table_rhizobia[, i] / colSums(otu_table)[i]}

# normalize rhizobial OTU table

otu_table_norm <- otu_table_rhizobia
for (i in 1:ncol(otu_table)) {otu_table_norm[, i] <- otu_table_rhizobia[, i] / colSums(otu_table_rhizobia)[i]}

# subset samples of interest

idx <- design$host=="arabidopsis_thaliana" & design$compartment=="leaf"
arabidopsis_leaf_design <- design[idx, ]
arabidopsis_leaf_design$category <- "arabidopsis_thaliana_leaf"

idx <- design$host=="arabidopsis_thaliana" & design$compartment=="root"
arabidopsis_root_design <- design[idx, ]
arabidopsis_root_design$category <- "arabidopsis_thaliana_root"

idx <- design$host=="arabidopsis_thaliana" & design$compartment=="rhizosphere"
arabidopsis_rhizosphere_design <- design[idx, ]
arabidopsis_rhizosphere_design$category <- "arabidopsis_thaliana_rhizosphere"

idx <- design$host=="lotus_japonicus" & design$compartment=="root" & design$genotype %in% c("gifu")
lotus_root_wt_design <- design[idx, ]
lotus_root_wt_design$category <- "lotus_japonicus_wt_root"

idx <- design$host=="lotus_japonicus" & design$compartment=="rhizosphere" & design$genotype %in% c("gifu")
lotus_rhizosphere_wt_design <- design[idx, ]
lotus_rhizosphere_wt_design$category <- "lotus_japonicus_wt_rhizosphere"

idx <- design$host=="lotus_japonicus" & design$compartment=="pooled_nodules" & design$genotype %in% c("gifu")
lotus_nodules_wt_design <- design[idx, ]
lotus_nodules_wt_design$category <- "lotus_japonicus_wt_nodules"

idx <- design$host=="lotus_japonicus" & design$compartment=="root" & design$genotype %in% c("nfr5_2", "nfr5_3")
lotus_root_mutant_design <- design[idx, ]
lotus_root_mutant_design$category <- "lotus_japonicus_mutant_root"

idx <- design$host=="lotus_japonicus" & design$compartment=="rhizosphere" & design$genotype %in% c("nfr5_2", "nfr5_3")
lotus_rhizosphere_mutant_design <- design[idx, ]
lotus_rhizosphere_mutant_design$category <- "lotus_japonicus_mutant_rhizosphere"

idx <- design$host=="barley" & design$compartment=="root"
barley_root_design <- design[idx, ]
barley_root_design$category <- "barley_root"

idx <- design$host=="barley" & design$compartment=="rhizosphere"
barley_rhizosphere_design <- design[idx, ]
barley_rhizosphere_design$category <- "barley_rhizosphere"

idx <- design$compartment=="soil"
soil_design <- design[idx, ]
soil_design$category <- "soil"

design <- rbind(
                arabidopsis_leaf_design,
                arabidopsis_root_design, arabidopsis_rhizosphere_design, 
                lotus_root_wt_design, lotus_rhizosphere_wt_design, lotus_nodules_wt_design,
                barley_root_design, barley_rhizosphere_design,
                soil_design)
                    
# match design table

design <- design[design$SampleID %in% colnames(otu_table_norm), ]
idx <- match(design$SampleID, colnames(otu_table_norm))
otu_table_norm <- otu_table_norm[, idx]
otu_table_norm_all <- otu_table_norm_all[, idx]

# match mapping table

mapping <- mapping[mapping$ID %in% rownames(otu_table_norm), ]
idx <- match(rownames(otu_table_norm), mapping$ID)
mapping <- mapping[idx, ]

# colors for plotting

cat <- c("arabidopsis_thaliana_leaf",
         "arabidopsis_thaliana_root", "arabidopsis_thaliana_rhizosphere",
         "lotus_japonicus_wt_root", "lotus_japonicus_wt_rhizosphere",
         "lotus_japonicus_wt_nodules", 
         "barley_root", "barley_rhizosphere",
         "soil")

col <- c("#64b200ff",
         "#ef67ebff", "#ef67ebff",
         "#00c1a7ff", "#00c1a7ff",
         "#00a6ffff",
         "#d4aa00ff", "#d4aa00ff",
         "#aa4400ff")

sha <- c(1,
         19, 3,
         19, 3,
         17,
         19, 3,
         18)

graphics <- data.frame(category=cat, color=col, shape=sha)

