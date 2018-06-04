
# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# aggregate rhizobial relative abundances per sample category

df <- design

df$ra <- colSums(otu_table_norm_all)
print(aggregate(df$ra, by=list(design$category), FUN=median))

# plotting

df$category <- factor(df$category, levels=cat)

p <- ggplot(df, aes(x=category, y=ra, color=category)) +
            geom_boxplot(alpha=1, outlier.shape=NA, size=0.7, width=1) +
            stat_summary(fun.y="mean", fun.ymin="mean", fun.ymax="mean", size=0.3, width=0.7,
                         geom="crossbar", position=position_dodge(width=0.7)) +
            geom_point(position=position_jitter(width=0.5), size=1, alpha=0.7) +
            scale_colour_manual(values=as.character(graphics$color)) +
            scale_y_continuous(labels=percent, limits=c(0, 1)) +
            coord_flip() +
            labs(x="", y="Aggregated Relative Abundance of Rhizobia") +
            main_theme

ggsave(paste(figures.dir, "aggregated_RA_boxplots.pdf", sep=""), p, width=9, height=5)

