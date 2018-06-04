#!/bin/bash

# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

getAmphoraGenes() {

    local genome_list=$(cut -f 1 $mapping_file | grep -v "ID")
  
    for genome_id in $genome_list 
    do
 
        log "["$genome_id"] retrieving AMPHORA genes..."

        # cleanup
        rm -f $working_dir/"$genome_id"_amphora.txt

        local amphora_list=$(cat $path_to_db/amphora/amphora_list.txt)

        for gene in $amphora_list
        do

            hmmsearch -E 0.001 \
                      -o $working_dir/"$genome_id"_"$gene".hmr \
                      $path_to_db/amphora/"$gene".hmm \
                      $working_dir/ORFs_AA/"$genome_id".faa \
                      1> /dev/null \
                      2> $output

            grep -o '....\.peg\.[0-9]*' $working_dir/"$genome_id"_"$gene".hmr | head -n 1 \
                 >> $working_dir/"$genome_id"_"$gene".txt

            rm -f $working_dir/"$genome_id"_"$gene".hmr

            touch $working_dir/"$gene"_"$genome_id".ffn

            n_seqs=$(wc -l $working_dir/"$genome_id"_"$gene".txt | sed 's/ .*//g')

            if [[ "$n_seqs" -gt 0 ]]
            then

                peg=$(cat $working_dir/"$genome_id"_"$gene".txt)
                
                echo ">"$genome_id"" > $working_dir/"$gene"_"$genome_id".ffn
          	    awk "/$peg$/ {flag=1;next} />/{flag=0} flag {print}" \
                    $working_dir/ORFs_NT/"$genome_id".ffn \
                    >> $working_dir/"$gene"_"$genome_id".ffn

            fi

            # cleanup
            rm -f $working_dir/"$genome_id"_"$gene".txt

        done

    done

}

alignAmphoraSeqs() {

    local amphora_list=$(cat $path_to_db/amphora/amphora_list.txt)

    for gene in $amphora_list
    do
    
        log "["$gene"] alignning gene..."

        cat $working_dir/"$gene"_*.ffn >> $working_dir/"$gene".ffn
        rm -f $working_dir/"$gene"_*.ffn $working_dir/"$gene"*.msa

        clustalo --seqtype=DNA \
                 --threads=$(nproc) \
                 -i $working_dir/"$gene".ffn \
                 -o $working_dir/"$gene".msa \
                 --full \
                 &>> $output
        
    done

}

concatenateAmphoraSeqs() {

    rm -f $working_dir/amphora.msa

    local genome_list=$(cut -f 1 $mapping_file | grep -v "ID")
    local amphora_list=$(cat $path_to_db/amphora/amphora_list.txt)

    for genome_id in $genome_list 
    do

    echo ">"$genome_id"" >> $working_dir/amphora.msa

        for gene in $amphora_list
        do 
            
            local n=$(grep -c ">" $working_dir/"$gene".msa)
            local n_tot=$(echo $genome_list | awk '{print NF}')

            # use only genes found in every genome
            if [[ "$n" -eq $n_tot ]]
            then 

                awk "/$genome_id$/ {flag=1;next} />/{flag=0} flag {print}" \
                    $working_dir/"$gene".msa \
                    >> $working_dir/amphora.msa

            fi

        done  
       
    done

    for gene in $amphora_list
    do 
 
        rm -f $working_dir/"$gene".msa $working_dir/"$gene".ffn

    done

}

buildMLTree() {

    # cleanup
    rm -f $working_dir/amphora.tree

    # generate ML ssp. tree
    FastTree -nt -gtr $working_dir/amphora.msa \
             1>> $working_dir/amphora_tree.newick \
             2>> $output

}

get16S() {

    local genome_id=$1

    # cleanup

    mkdir -p $working_dir/16S

    rm -f $working_dir/16S/"$genome_id".fasta

    rnammer -S bac \
            -m ssu \
            -f $working_dir/16S/"$genome_id".fasta \
            $data_dir/assemblies/"$genome_id".fna \
            -h /dev/null \
            &>> $output

    cat $working_dir/16S/"$genome_id".fasta | \
        awk 'BEGIN{n=0} />/ {n++; print ">" n $3} !/>/ {print $0}' | \
        sed 's/>/>'$genome_id'_/g;s/.score=/ /g' \
        >> $working_dir/16S/16S_all.fasta

}

getSequences() {

    local N=$1
    
    log "[gfam_"$N"] retrieving sequences..."

    mkdir -p $working_dir/gfams_AA \
             $working_dir/gfams_NT

    rm -f $working_dir/gfams_AA/gfam_"$N".faa \
          $working_dir/gfams_NT/gfam_"$N".ffn

    touch $working_dir/gfams_AA/gfam_"$N".faa \
          $working_dir/gfams_NT/gfam_"$N".ffn

    pegs=$(sed -n "$N"p $data_dir/orthogroups.txt | sed 's/gfam_[0-9]*.//g;s/[0-9]*|//g')

    if [ ! -z "$pegs" ]; then
  
        for i in $pegs
        do
  
        	org=(`echo $i | grep -o '[^\.]*\.' `)
  
        	echo ">"$i >> $working_dir/gfams_AA/gfam_"$N".faa
        	awk "/$i$/ {flag=1;next} /peg/{flag=0} flag {print}" \
                $data_dir/ORFs_AA/$org"faa" \
                1>> $working_dir/gfams_AA/gfam_"$N".faa \
                2>> $output
  
        	echo ">"$i >> $working_dir/gfams_NT/gfam_"$N".ffn
        	awk "/$i$/ {flag=1;next} /peg/{flag=0} flag {print}" \
                $data_dir/ORFs_NT/$org"ffn" \
                1>> $working_dir/gfams_NT/gfam_"$N".ffn \
                2>> $output
  
        done
  
    fi

}

alignSequences() {

    local N=$1

    log "[gfam_"$N"] aligning sequences..."
    
    mkdir -p $working_dir/MSAs

    rm -f $data_dir/MSAs/gfam_"$N".msa

    clustalo --seqtype=DNA \
             --threads=$n_cores \
             -i $working_dir/gfams_NT/gfam_"$N".ffn \
             -o $data_dir/MSAs/gfam_"$N".msa \
             --full \
             &>> $output
    
}

generateTree() {

    local N=$1

    log "[gfam_"$N"] generating tree..."

    mkdir -p $working_dir/trees

    rm -f $working_dir/trees/gfam_"$N".tree \
          $working_dir/trees/gfam_"$N"_adjusted.tree

    $scripts_dir/reduceMSA.R $working_dir/MSAs/gfam_"$N".msa \
                             $working_dir/MSAs/gfam_"$N"_reduced.msa \
                             &>> $output

    $scripts_dir/reduceTree.R $working_dir/MSAs/gfam_"$N"_reduced.msa \
                              $data_dir/spp_rooted.tree \
                              $working_dir/trees/gfam_"$N".tree \
                              &>> $output

    FastTree -nt -gtr -nome -mllen -intree \
             $working_dir/trees/gfam_"$N".tree \
             $working_dir/MSAs/gfam_"$N"_reduced.msa \
             1>> $working_dir/trees/gfam_"$N"_adjusted.tree \
             2>> $output

    mv $working_dir/trees/gfam_"$N"_adjusted.tree $working_dir/trees/gfam_"$N".tree

    #~ sed -i 's/)[0-9.]*/)/g' $data_dir/trees/gfam_"$N".tree

}

gfamPhylo() {

    local N=$1

    getSequences $N

    if [[ $(grep -c ">" $working_dir/gfams_AA/gfam_"$N".faa) -eq 1 ]] 
    then

        log "[gfam_"$N"] only one sequence, skipping..."
        mv $working_dir/ORFs_NT/gfam_"$N".ffn $data_dir/MSAs/gfam_"$N".msa

    else if [[ $(grep -c ">" $working_dir/gfams_AA/gfam_"$N".faa) -eq 0 ]] 
    then

        log "[gfam_"$N"] no valid sequences, skipping..."

    else
        
        alignSequences $N

        generateTree $N

        log "[gfam_"$N"] done "

    fi
    fi

}


gfamAnnotation() {

    local N=$1

    log "[gfam_"$N"] extracting annotations..."

    mkdir -p $working_dir/gfams_annotations

    rm -f $working_dir/gfams_annotations/gfam_"$N"_KEGG.txt \
          $working_dir/gfams_annotations/gfam_"$N"_KEGG.txt_temp \
          $working_dir/gfams_annotations/gfam_"$N".ko

    local pegs=$(sed -n "$N"p $data_dir/orthogroups.txt | cut -d ' ' -f 1 --complement)

    for peg in $pegs
    do 
        
        local genome=$(echo $peg | sed 's/|.*//g')
        
        cat $data_dir/annotations/"$genome"_KEGG.txt | \
            grep ''$peg'[^0-9]' $data_dir/annotations/"$genome"_KEGG.txt | \
            cut -f 1 --complement | \
            sed '/^\t*$/d' \
            >> $working_dir/gfams_annotations/gfam_"$N"_KEGG.txt_temp

    done

    cat $working_dir/gfams_annotations/gfam_"$N"_KEGG.txt_temp | \
        sort | uniq | sed 's/^/gfam_'$N'\t/g' \
        >> $working_dir/gfams_annotations/gfam_"$N"_KEGG.txt

    cat $working_dir/gfams_annotations/gfam_"$N"_KEGG.txt | \
        cut -f 1,5 | sort | uniq \
        >> $working_dir/gfams_annotations/gfam_"$N".ko

    rm -f $working_dir/gfams_annotations/gfam_"$N"_KEGG.txt_temp

}

