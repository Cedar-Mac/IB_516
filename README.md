# IB_516 Course Materials and Project

This repo is a baby project. Most of the material contained here is simply course materials.
However, for the little nugget of a project check the following subdirectories:
- code
- results
- data

The goal of the project is to implement a pipeline for filtering and denoising fasta sequences by removing erroneous sequences, chimeric sequences, removing singletons, clustering into OTU's based on observed sequence similarity, etc.

The "novel" portion will be optimizing the entropy ratio of 2nd and 3rd codon positions (measure of sequence integrity) and number of sequences removed.
If that goes well, then I may also include some network analysis by implementing a network flavored version of the Speceficity and Speceficity Diversity metric.

