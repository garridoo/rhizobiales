#!/bin/bash

# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# exits whenever a function returns 1
set -e
# exits if unset variables are used
set -o nounset

# parse arguments
config_file=$1

# get path to scripts
scripts_dir=$(dirname $0)

# check paths file
if [ ! -f $config_file ]
then
    echo "invalid config file"
    exit 1
fi

# load paths
source $config_file

# load functions
source $scripts_dir/compgen.functions.sh
source $scripts_dir/datagen.functions.sh
source $scripts_dir/phylo.functions.sh
source $scripts_dir/parallel.sh

genomes=$(cut -f 1 $mapping_file | grep -v "ID")

### gene calling

for genome_id in $genomes

    do echo "geneCalling" $genome_id

done > $working_dir/datagen_parallel_list.txt

parallel $working_dir/datagen_parallel_list.txt $n_cores
rm -f $working_dir/datagen_parallel_list.txt

### parse proteomes

for genome_id in $genomes

    do echo "prepareProteome" $genome_id

done > $working_dir/datagen_parallel_list.txt

parallel $working_dir/datagen_parallel_list.txt $n_cores
rm -f $working_dir/datagen_parallel_list.txt

### annotate proteomes

for genome_id in $genomes

    do echo "annotateGenomeKEGG" $genome_id

done > $working_dir/datagen_parallel_list.txt

parallel $working_dir/datagen_parallel_list.txt $n_cores
rm -f $working_dir/datagen_parallel_list.txt

### extract 16S rRNA sequences

rm -f $working_dir/16S_all.fasta

for genome_id in $genomes

    do echo "get16S" $genome_id

done > $working_dir/datagen_parallel_list.txt

parallel $working_dir/datagen_parallel_list.txt $n_cores
rm -f $working_dir/datagen_parallel_list.txt

cat $working_dir/16S/* >> $working_dir/16S_all.fasta
mv $working_dir/16S_all.fasta $working_dir/16S/

### generate species tree

# extract amphora genes

for genome_id in $genomes

    do echo "getAmphoraGenes" $genome_id

done > $working_dir/datagen_parallel_list.txt

parallel $working_dir/datagen_parallel_list.txt $n_cores
rm -f $working_dir/datagen_parallel_list.txt

# align amphora genes independently
alignAmphoraSeqs

# concatenate MSA
concatenateAmphoraSeqs

# build Maximum-Likelihood tree
buildMLTree

# parse tree and rename tip labels
$scripts_dir/phylo.R $mapping_file \
                     $working_dir/amphora_tree.newick \
                     $working_dir

log "DONE!"

