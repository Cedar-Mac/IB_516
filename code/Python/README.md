# General Description of Functions

These are the python scripts for merging, filtering, denoising, and clustering short read sequences for metabarcoding.
Generally each script is a function, and functions are called in order in the ./main.py executable.

Here is a general overview of each step (refer to the function definitions for further info):

- merge_pairs.py: Merge paired end reads from illumina sequencing.
- trim_primers.py: Remove forward and reverse primer sequences.
- quality_filter.py: Filter reads based on average expected error.
- length_filter.py: Filter reads to be exact amplicon length. Reformat as fasta file
- chimera_filter.py: Remove chimeric sequences de novo (no ref sequences)
- frequency_filter.py: Filter sequences if it occurs more than once at a site, or occurs at more than one site.
- denoise.py:
- cluster.py:
