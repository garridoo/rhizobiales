# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# PCoA Bray-Curtis

k <- 2
bray_curtis <- vegdist(t(otu_table_norm_all), method="bray")
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID), ])


# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=category, shape=compartment)) +
     geom_point(alpha=.7, size=2) +
     scale_colour_manual(values=c(c_red, c_dark_brown, c_green, c_very_dark_green, c_grey)) +
     scale_shape_manual(values=c(18, 16, 1, 3)) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
     main_theme +
     theme(legend.position="right")

ggsave(paste(figures.dir, "PCoA_BC.pdf", sep=""), p, height=7, width=10)


# NMDS Bray-Curtis

k <- 2
bray_curtis <- vegdist(t(otu_table_norm_all), method="bray")
pcoa <- isoMDS(bray_curtis, k=k)
points <- pcoa$points
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID), ])


# plot NMDS 1 and 2

points$category <- factor(points$category, levels=cat)
points$host <- as.character(points$host)
points$host[points$category=="lotus_japonicus_wt_nodules"] <- "lotus_japonicus_wt_nodules"
points$host[points$category=="arabidopsis_thaliana_leaf"] <- "arabidopsis_thaliana_leaf"


p <- ggplot(points, aes(x=x, y=y, color=host, shape=host)) +
     geom_point(alpha=.7, size=1.5) +
     stat_ellipse(type="norm") +
     scale_colour_manual(values=as.character(graphics$color)) +
     scale_shape_manual(values=graphics$shape) +
     labs(x="NMDS 1", y="NMDS 2") +
     main_theme +
     theme(legend.position="right")

ggsave(paste(figures.dir, "NMDS_BC.pdf", sep=""), p, height=7, width=10)

sha <- c(25,
         19, 3,
         19, 3,
         17,
         19, 3,
         18)

graphics <- data.frame(category=cat, color=col, shape=sha)

idx <- rownames(otu_table_norm_all) %in% mapping$ID[mapping$tax=="Methylobacterium"]
methylo <- colSums(otu_table_norm_all[idx, ])
points$methylo <- methylo
points$methylo[points$methylo > 0.05] <- 0.05

idx <- rownames(otu_table_norm_all) %in% mapping$ID[mapping$tax=="Mesorhizobium"]
meso <- colSums(otu_table_norm_all[idx, ])
points$meso <- meso

p <- ggplot(points, aes(x=x, y=y, color=methylo, fill=methylo, shape=category)) +
     geom_point(alpha=.7, size=1.5) +
     scale_colour_gradient(high="#2e815a", low="grey") +
     scale_fill_gradient(high="#2e815a", low="grey") +
     scale_shape_manual(values=graphics$shape) +
     labs(x="NMDS 1", y="NMDS 2") +
     main_theme +
     theme(legend.position="right")

ggsave(paste(figures.dir, "NMDS_BC_methylo.pdf", sep=""), p, height=7, width=10)

p <- ggplot(points, aes(x=x, y=y, color=meso, fill=meso, shape=category)) +
     geom_point(alpha=.7, size=1.5) +
     scale_colour_gradient(high="#32dc32", low="grey") +
     scale_fill_gradient(high="#32dc32", low="grey") +
     scale_shape_manual(values=graphics$shape) +
     labs(x="NMDS 1", y="NMDS 2") +
     main_theme +
     theme(legend.position="right")

ggsave(paste(figures.dir, "NMDS_BC_meso.pdf", sep=""), p, height=7, width=10)

