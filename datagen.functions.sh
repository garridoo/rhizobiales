#!/bin/bash

# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

geneCalling() {

    local genome_id=$1
    local assembly=$data_dir/assemblies/"$genome_id".fna

    log "["$genome_id"] gene calling..."

    # cleanup
    
    mkdir -p $working_dir/ORFs_NT \
             $working_dir/ORFs_AA \
             $working_dir/ORFs_coords

    rm -f $working_dir/ORFs_NT/"genome_id".ffn \
          $working_dir/ORFs_AA/"genome_id".faa

    # use PRODIGAL for gene prediction
    prodigal -i $assembly \
             -a $working_dir/ORFs_AA/"$genome_id".faa \
             -d $working_dir/ORFs_NT/"$genome_id".ffn \
             -o $working_dir/ORFs_coords/"$genome_id".gff \
             -f gff \
             1>> $output \
             2>> /dev/null \

    # remove stop codons
    sed -i 's/\*//g' $working_dir/ORFs_AA/"$genome_id".faa

}

prepareProteome() {

    local genome_id=$1
    local proteome_nt=$working_dir/ORFs_NT/"$genome_id".ffn
    local proteome_aa=$working_dir/ORFs_AA/"$genome_id".faa
 
    log "["$genome_id"] parsing proteome fasta files..."
 
    # cleanup
    rm -f $working_dir/"$genome_id"_compliant.ffn \
          $working_dir/"$genome_id"_compliant.faa

    # parse NT proteome
    awk -v org=$genome_id 'BEGIN {n=1}; />/ {print ">"org"|"org".peg."n; n++} !/>/ {print}' $proteome_nt \
        >> $working_dir/"$genome_id"_compliant.ffn
     mv $working_dir/"$genome_id"_compliant.ffn $proteome_nt
    
    # parse AA proteome
    awk -v org=$genome_id 'BEGIN {n=1}; />/ {print ">"org"|"org".peg."n; n++} !/>/ {print}' $proteome_aa \
        >> $working_dir/"$genome_id"_compliant.faa
     mv $working_dir/"$genome_id"_compliant.faa $proteome_aa

}

