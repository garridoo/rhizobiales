## Scripts of Garrido-Oter *et. al*, Modular traits of the Rhizobiales root microbiota and their evolutionary relationship with symbiotic rhizobia

originally by Ruben Garrido-Oter

garridoo@mpipz.mpg.de

These scripts are made available to facilitate reproducibility of our research. If you re-use any or part of this code, please reference with comments and cite our paper!

### Scripts used for processing data, creating the figures and performing the statistical analysis reported in the manuscript.

#### Whole-genome assembly, quality control, gene calling, annotation and phylogeny inference:

[config.sh](path)

Configuration script containing paths to scripts and data as well as parameters for custom scripts and third-party tools.

[assembly.functions.sh](path)

Bash script containing ancilliary functions for whole-genome assembly.

[assembly.sh](path)

Script used to assemble the genomes using SOAPdenovo and A5. It can be run in parallel using either the custom script bellow or the gnu parallel suit.

[parallel.sh](path)

Custom script to run bash functions in parallel in a multi-core machine.

[assembly_stats.R](path)

R script used to generate assembly statistics as well as GC and k-mer spectral projections. The output of this script contains clean assemblies (all contigs <1,000 bp are removed) as well as a PDF file containing a report which was used to manually inspect for likely contaminated assemblies.

[datagen.functions.sh](path/datagen.functions.sh), [compgen.functions.sh](path/compgen.functions.sh), and [phylo.functions.sh](path/phylo.functions.sh)

Bash scripts with ancilliary functions for gene calling and parsing of assemblies and ORF fasta files.

[datagen.sh](path/datagen.sh)

Bash script to perform parallel (opt. serial) processing of the whole-genome assemblies, including gene calling, quality filtering, parsing of the FASTA files, annotation, and phylogeny inference.

[paths.R](path)

R script containing paths to data, intermediate results and output folders for the figures and statistical analyses.

[plotting_functions.R](path) and [colors.R](path)

R scripts containing functions and variables using for plotting (e.g. colors, ggplot2 themes, etc.)

[functional_MDS_KEGG.R](path)

R script to perform dimmensionality reduction of genome functional profiles and plotting.

#### Reference-based 16S rRNA amplicon profiling meta-analysis of Rhizobiales community structure and diversity:

[meta.sh](path)

Script to process 16S amplicon data from previous studies from raw data using a standarized pipeline. The output consists of joint OTU tables, FASTA files containing reference sequences, taxonomy classification, as well as multiple alpha- and beta-diversity estimates.

[16s.functions.sh](path)

Auxiliary functions to perform meta-analysis of 16S data.

[config.sh](path)

Configuration script containing paths to scripts and data as well as parameters for custom scripts and third-party tools.

[bulgarelli.sh](path)

Script to pre-process data from Schlaeppi et al., 2012 (Arabidopsis root and rhizosphere).

[schlaeppi_2014.sh](path)

Script to pre-process data from Schlaeppi et al., 2014 (Arabidopsis and relatives root and rhizosphere.

[bulgarelli_2015.sh](path)

Script to pre-process data from Schlaeppi et al., 2015 (Barley root and rhizosphere).

[bai_2015_root.sh](path)

Script to pre-process data from Bai et al., 2015 (Arabidopsis root).

[bai_2015_leaf.sh](path)

Script to pre-process data from Bai et al., 2014 (Arabidopsis leaf).

[zgadzaj_2016.sh](path)

Script to pre-process data from Zgadzaj et al., 2016 (Lotus root, rhizosphere and nodules).

[normalize_otu_table.R](path)

R script to perform CSS normalization of the joint OTU table.

[16S_meta_analysis.R](path)

Script to process data from joint analysis of 16S published studies.

[boxplots.R](path)

Script to plot cumulative relative abundances of rhizobiales across hosts and micro-habitats.

[PCoA.R](path)

Script to plot beta-diversity of rhizobiales communities across hosts and micro-habitats.

---------------------------

For any questions regarding these scripts, please contact

Ruben Garrido-Oter

garridoo@mpipz.mpg.de

