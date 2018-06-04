
# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list = ls())

# load libraries

library(utils, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(MASS, quietly=T, warn.conflicts=F)
library(gridExtra, quietly=T, warn.conflicts=F)
library(scales, quietly=T, warn.conflicts=F)

options(warn=-1)

# plotting functions, etc.

source("plotting_functions.R")

# load paths to project directories

source("paths.R")

mapping.file <- paste(data.dir, "mapping.txt", sep="")
taxonomy.file <- paste(data.dir, "taxonomy.txt", sep="")

# load data

mapping <- read.table(mapping.file, sep="\t", header=T, colClasses="character")
taxonomy <- read.table(taxonomy.file, sep="\t", header=T)

### functional profiles

message("generating matrix of functional profiles...")

ko.all <- data.frame(genome=NULL, ko=NULL) 
sizes.all <- data.frame(genome=NULL, size=NULL) 

pb <- txtProgressBar(min=1, max=length(mapping$ID), style=3)
i <- 1

for (g in mapping$ID) {

    setTxtProgressBar(pb, i)
    i <- i + 1
   
    ko <- read.table(paste(annotation.dir, g, ".ko", sep=""),
                     fill=T, header=F, sep="\t",
                     col.names=c("peg", "ko"))[, 2]
    ko.genome <- data.frame(genome=g, ko=ko)
    ko.all <- rbind(ko.all, ko.genome)

    size.genome <- data.frame(genome=g, size=dim(ko.genome)[1])
    sizes.all <- rbind(sizes.all, size.genome)

}

close(pb)

ko.table <- table(ko.all)
ko.table <- t(ko.table[, -1])

sizes.all$perc_annotated <- colSums(ko.table) / sizes.all$size

func <- (ko.table > 0) * 1
#~ func <- ko.table

write.table(func, file=paste(data.dir, "/functional_profiles.txt", sep=""),
            sep="\t", quote=F, col.names=T, row.names=T)


func <- read.table(paste(data.dir, "/functional_profiles.txt", sep=""), sep="\t", header=T, check.names=F)


# calculate pairwise functional distances

message("calculating pairwise functional distances...")

d <- 1 - cor(func)
diag(d) <- 0

### PCoA of functional distances

message("calculating functional PCoA...")

k <- 2

pcoa <- cmdscale(d, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig

points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points$compartment <- mapping$compartment[match(rownames(points), mapping$ID)] 
points$taxonomy <- mapping$taxonomy[match(rownames(points), mapping$ID)] 
points$host <- mapping$host[match(rownames(points), mapping$ID)] 
points$nifh <- mapping$nifh[match(rownames(points), mapping$ID)] 

source("colors.R")

points$taxonomy <- factor(points$taxonomy, levels=colors$group)

host.colors <- data.frame(host=c("Arabidopsis", "Legume"       , "Maple",
                                 "Corn"       , "Soybean"      , "Nematode",
                                 "Insect"     , "Other"        , "Rice",
                                 "Wheat"      , "Oats"         , "Cotton",
                                 "Cucumber"   , "Soil"))

          host.colors$color <- c("red"       , "darkgreen"    , "pink",
                                 "orange"    , "green"        , "darkred",
                                 "black"     , "grey"         , "darkblue",
                                 "yellow"    , "purple"       , "blue",
                                 "green"     , "brown")


points$host <- factor(points$host, levels=host.colors$host)


p1 <- ggplot(points, aes(x=x, y=y, color=taxonomy, shape=nifh)) +
      geom_point(alpha=.7, size=1) +
      scale_shape_manual(values=c(16, 3)) +
      scale_colour_manual(values=as.character(colors$color)) +
      scale_colour_manual(values=as.character(host.colors$color)) +
      labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
      main_theme +
      theme(legend.title=element_blank(),
            legend.position="top")

ggsave(file=paste(figures.dir, "functional_MDS.pdf", sep=""), p1, height=8, width=8)

