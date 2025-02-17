# Overview of Code

There are two subfolders in this directory:

* [Python](/Python)
* [R](/R)

Most of the R code is legacy from previous work done for J. Freedman (mostly recreating his own work to gain better personal familiarity with the pipeline). This code is not actively managed and is included simply for my personal reference.

The [Python](/Python) subfolder contains all the novel code for this project. There are many scripts with a single function per script I've also consolidated like functions into their own source scripts. the main ones are:

1. [pipeline.py](/Python/pipeline.py): All functions that are part of the basic metabarcoding sequence process including:

    * Merging sequences from paired end reads.
    * Trimming primers.
    * Filtering for exact amplicon length.
    * Filtering out singletons (or otherwise rare sequences).
    * Filtering by quality scores.
    * Filter chimeric sequences.
    * Denoise (various mehtods implemented).

2. [alignments.py](/Python/alignments.py) Functions necessary to find conserved base positions by aligning amplicon to a reference set of sequences (see Pentasari et al.) Functions include:

    * Aligning reference sequences of amino acids.  
    * Finding the conserved amino acids in metazoan reference sequences.
    * Align the amplicon of study to the reference alignment. (convert to amino acid sequence first)
    * Get the conserved nucleotide base positions from the aligned amplicon amino acid sequence.

3. [denoising_stats.py](/Python/denoising_stats.py): All functions related to getting statistics for each denoising method, including:

    * Total true errors
    * Kept sequences
    * Error rate (total true errors / kept sequences)
    * Average entropy ratio (2nd base pos. / 3rd base pos.)

Most of the individual functions have also been broken out into their own scripts (may delete later because it seems redundant, was originally for unit testing).

There may be more to come. Currently graphing functions are included in the denoising_stats.py script. Additionally, I may add new stats scripts for various stages of the pipeline, or add that functionality to current functions.
