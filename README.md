# IB_516 Project

This repo is a baby project. The goal of the project is two-fold:

1. Implement entropy based denoising of a COI amplicon.
2. Benchmark various denoising methods (unoise3 algorithm, DADA2 and entropy penalized)
3. Bonus: If time, explore other changes to the sequence processing pipeline (besides denoising)

Here is the structure of the project:

* [code](/code/)
  * [Python](/code/Python/)
  * [R](/code/R/)
* [tmp_plots](/tmp_plots/)
* [results](/results/)
  * [final_figs/](results/final_figs)
* [Readings](/Readings)
* data (currently not published)

The "novel" portion of this project will be optimizing denoising parameters based on "true errors." Or, in otherwords, errors which disrupt the reading frame (stop codons) or alter amino acids that are coserved across all metazoans.

