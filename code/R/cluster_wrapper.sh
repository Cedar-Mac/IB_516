#!/bin/sh    

# This is a wrapper for the --cluster_size clustering algorithm of vsearch. 
# The only argument is the --id parameter, the clustering distance between strands in a cluster


VSEARCH=$(which vsearch)
THREADS=4
INPUT_FILE=/Users/cedarmackaness/gc/F_Cluster_otus/_data/2_OTU_clustering/C_all_derep_min2_OTUs+chimeras.fasta
ALPHA=$1


function cluster_wrapper {
    $VSEARCH --cluster_unoise $INPUT_FILE \
    --threads $THREADS \
    --ALPHA $ALPHA \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc id_${ID}_min2_OTUs_no_chimeras.uc \
    --relabel OTU_ \
    --centroids id_${ID}_min2_OTUs_no_chimeras.fasta \
    --otutabout id_${ID}_min2_OTUs_no_chimeras.otutab.txt
}

cluster_wrapper $ALPHA