## Scripts of Garrido-Oter *et. al*, Modular traits of the Rhizobiales root microbiota and their evolutionary relationship with symbiotic rhizobia

originally by Ruben Garrido-Oter

garridoo@mpipz.mpg.de

These scripts are made available to facilitate reproducibility of our research. If you re-use any or part of this code, please reference with comments and cite our paper! Raw data and intermediate results necesary to run these scripts can be downloaded [here](http://www.mpipz.mpg.de/R_scripts).

### Scripts used for processing data, creating the figures and performing the statistical analysis reported in the manuscript.

#### Whole-genome assembly, quality control, gene calling, annotation and phylogeny inference:

[config.sh](https://github.com/garridoo/rhizobiales)

Configuration script containing https://github.com/garridoo/rhizobialess to scripts and data as well as parameters for custom scripts and third-party tools.

[assembly.functions.sh](https://github.com/garridoo/rhizobiales)

Bash script containing ancilliary functions for whole-genome assembly.

[assembly.sh](https://github.com/garridoo/rhizobiales)

Script used to assemble the genomes using SOAPdenovo and A5. It can be run in parallel using either the custom script bellow or the gnu parallel suit.

[parallel.sh](https://github.com/garridoo/rhizobiales)

Custom script to run bash functions in parallel in a multi-core machine.

[assembly_stats.R](https://github.com/garridoo/rhizobiales)

R script used to generate assembly statistics as well as GC and k-mer spectral projections. The output of this script contains clean assemblies (all contigs <1,000 bp are removed) as well as a PDF file containing a report which was used to manually inspect for likely contaminated assemblies.

[datagen.functions.sh](https://github.com/garridoo/rhizobiales/datagen.functions.sh), [compgen.functions.sh](https://github.com/garridoo/rhizobiales/compgen.functions.sh), and [phylo.functions.sh](https://github.com/garridoo/rhizobiales/phylo.functions.sh)

Bash scripts with ancilliary functions for gene calling and parsing of assemblies and ORF fasta files.

[datagen.sh](https://github.com/garridoo/rhizobiales/datagen.sh)

Bash script to perform parallel (opt. serial) processing of the whole-genome assemblies, including gene calling, quality filtering, parsing of the FASTA files, annotation, and phylogeny inference.

[https://github.com/garridoo/rhizobialess.R](https://github.com/garridoo/rhizobiales)

R script containing https://github.com/garridoo/rhizobialess to data, intermediate results and output folders for the figures and statistical analyses.

[plotting_functions.R](https://github.com/garridoo/rhizobiales) and [colors.R](https://github.com/garridoo/rhizobiales)

R scripts containing functions and variables using for plotting (e.g. colors, ggplot2 themes, etc.)

[functional_MDS_KEGG.R](https://github.com/garridoo/rhizobiales)

R script to perform dimmensionality reduction of genome functional profiles and plotting.

#### Reference-based 16S rRNA amplicon profiling meta-analysis of Rhizobiales community structure and diversity:

[meta.sh](https://github.com/garridoo/rhizobiales)

Script to process 16S amplicon data from previous studies from raw data using a standarized pipeline. The output consists of joint OTU tables, FASTA files containing reference sequences, taxonomy classification, as well as multiple alpha- and beta-diversity estimates.

[16s.functions.sh](https://github.com/garridoo/rhizobiales)

Auxiliary functions to perform meta-analysis of 16S data.

[config.sh](https://github.com/garridoo/rhizobiales)

Configuration script containing https://github.com/garridoo/rhizobialess to scripts and data as well as parameters for custom scripts and third-party tools.

[bulgarelli.sh](https://github.com/garridoo/rhizobiales)

Script to pre-process data from Schlaeppi et al., 2012 (Arabidopsis root and rhizosphere).

[schlaeppi_2014.sh](https://github.com/garridoo/rhizobiales)

Script to pre-process data from Schlaeppi et al., 2014 (Arabidopsis and relatives root and rhizosphere.

[bulgarelli_2015.sh](https://github.com/garridoo/rhizobiales)

Script to pre-process data from Schlaeppi et al., 2015 (Barley root and rhizosphere).

[bai_2015_root.sh](https://github.com/garridoo/rhizobiales)

Script to pre-process data from Bai et al., 2015 (Arabidopsis root).

[bai_2015_leaf.sh](https://github.com/garridoo/rhizobiales)

Script to pre-process data from Bai et al., 2014 (Arabidopsis leaf).

[zgadzaj_2016.sh](https://github.com/garridoo/rhizobiales)

Script to pre-process data from Zgadzaj et al., 2016 (Lotus root, rhizosphere and nodules).

[normalize_otu_table.R](https://github.com/garridoo/rhizobiales)

R script to perform CSS normalization of the joint OTU table.

[16S_meta_analysis.R](https://github.com/garridoo/rhizobiales)

Script to process data from joint analysis of 16S published studies.

[boxplots.R](https://github.com/garridoo/rhizobiales)

Script to plot cumulative relative abundances of rhizobiales across hosts and micro-habitats.

[PCoA.R](https://github.com/garridoo/rhizobiales)

Script to plot beta-diversity of rhizobiales communities across hosts and micro-habitats.

---------------------------

For any questions regarding these scripts, please contact

Ruben Garrido-Oter

garridoo@mpipz.mpg.de

