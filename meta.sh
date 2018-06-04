#!/bin/bash

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# exits whenever a function returns 1
set -e

# get path to scripts
scripts_dir=$(dirname $0)

# load config file
source $scripts_dir/config.sh

# load functions
source $scripts_dir/16s.functions.sh

# activate QIIME, etc.
source $scripts_dir/activate.sh

# # process data from the different studies
# 
# log "processing data from Schlaeppi et al., 2014..."
# $scripts_dir/schlaeppi_2014.sh
# 
# log "processing data from Bulgarelli et al., 2015..."
# $scripts_dir/bulgarelli_2015.sh
# 
# log "processing data from Bai et al., 2015..."
# $scripts_dir/bai_2015_root.sh
# $scripts_dir/bai_2015_leaf.sh
# 
# log "processing data from Zgadzaj et al., 2016..."
# $scripts_dir/zgadzaj_2016.sh
# 
log "combining pre-processed reads from all studies..."

# load config file
source $scripts_dir/config.sh

rm -f $working_dir/seqs.fasta
cat $working_dir/schlaeppi_2014/seqs.fasta >> $working_dir/seqs.fasta 
cat $working_dir/bulgarelli_2015/seqs.fasta >> $working_dir/seqs.fasta 
cat $working_dir/bai_2015_root/seqs.fasta >> $working_dir/seqs.fasta 
cat $working_dir/bai_2015_leaf/seqs.fasta >> $working_dir/seqs.fasta 
cat $working_dir/zgadzaj_2016/seqs.fasta >> $working_dir/seqs.fasta 

# dereplication
log "dereplicating..."
usearch -derep_fulllength $working_dir/seqs.fasta \
        -fastaout $working_dir/seqs_unique.fasta \
        -sizeout \
        &>> $output

# abundance sort and discard singletons
log "sorting by abundance and discarding singletons..."
usearch -sortbysize $working_dir/seqs_unique.fasta \
        -fastaout $working_dir/seqs_unique_sorted.fasta \
        -minsize $min_size \
        &>> $output

# OTU clustering
log "OTU clustering using UPARSE..."
usearch -cluster_otus $working_dir/seqs_unique_sorted.fasta \
        -otus $working_dir/otus.fasta \
        -id $id_threshold \
        &>> $output

# chimera detection
log "removing chimeras..."
usearch -uchime_ref $working_dir/otus.fasta \
        -db $gold_db \
        -strand plus \
        -nonchimeras $working_dir/otus_nc.fasta \
        -threads $n_cores \
        &>> $output

# align sequences to database using PyNAST and remove remaining
log "aligning OTU representative sequences to database..."
align_seqs.py -i $working_dir/otus_nc.fasta \
              -t $gg_core_aligned_db \
              -p $min_id_aln \
              -o $working_dir

sed -i 's/-//g' $working_dir/otus_nc_aligned.fasta

# rename OTUs and remove alignment gaps
log "renaming OTUs..."

cat $working_dir/otus_nc_aligned.fasta | \
    awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' \
    >> $working_dir/rep_seqs.fasta

# generate OTU table
log "generating OTU table..."
usearch -usearch_global $working_dir/seqs.fasta \
        -db $working_dir/rep_seqs.fasta \
        -strand plus \
        -id $id_threshold \
        -uc $working_dir/read_mapping.uc \
        &>> $output

# convert UC file to txt
log "converting UC OTU table file into text format..."
python $usearch_dir/uc2otutab.py $working_dir/read_mapping.uc \
    1> $working_dir/otu_table.txt \
    2>> $output

# reference-based OTU clustering using the UPARSE-REF algorithm
# using the rhizobia representative sequences as a reference
log "reference-based OTU clustering using UPARSE-REF..."
usearch -uparse_ref $working_dir/seqs_unique_sorted.fasta \
        -db $data_dir/rhizobia_100_rep_seqs.fasta \
        -strand plus \
        -uparseout $working_dir/otu_clustering_rhizobia_97.up \
        -threads $n_cores \
        &>> $output

# generate UC-compatible mapping file for USEARCH compatibility at 97%
rm -f $working_dir/otu_clustering_rhizobia_97.uc
cat $working_dir/otu_clustering_rhizobia_97.up | \
    awk -v i=0.97 '{if ($3<i*100 || $2=="chimera" || $2=="other") {$5=$2}
                    print $1, $5}' | \
    sed 's/ /\t/g;s/^/H\t\t\t\t\t\t\t\t/g' \
    1>> $working_dir/otu_clustering_rhizoia_97.uc \
    2>> $output

# generate OTU table and filter sequences
log "generating OTU table..."
usearch -usearch_global $working_dir/seqs.fasta \
        -db $data_dir/rhizobia_100_rep_seqs.fasta \
        -strand plus \
        -id $id_threshold \
        -uc $working_dir/otu_clustering_rhizobia_97.uc \
        -matched $working_dir/rhizobia_seqs.fasta \
        &>> $output

# convert UC file to txt
log "converting UC OTU table file into text format..."
python $usearch_dir/uc2otutab.py $working_dir/otu_clustering_rhizobia_97.uc \
    1> $working_dir/otu_table_rhizobia_97.txt \
    2>> $output

log "DONE!"

