# IB_516 Course Materials and Project

This repo is a baby project. Most of the material contained here is simply course materials.
However, for the little nugget of a project check the following subdirectories:
- code
- results
- data

The goal of the project is to implement a pipeline for filtering and denoising fastq sequences by merging paired end reads, filtering by length and quality scores, filtering by entropy scores and by "codon integrity" (no stop codons, and 4 conserved amino acids), removing chimeric sequences, removing singletons, clustering into OTU's based on observed sequence similarity, and assigning taxon information from an uncurated database (BOLD database).

The "novel" portion will be optimizing the entropy ratio of 2nd and 3rd codon positions (measure of sequence integrity) and number of sequences removed.
If that goes well, then I may also include some network analysis by implementing a network flavored version of the Specificity and Specificity Diversity metric.

