import subprocess
import os
import sys


sys.argv[1] = "COI"
sys.argv[2] = 1

barcode_type = sys.argv[1]

barcode_option = int(sys.argv[2])

if barcode_type == "COI":
    if barcode_option == 1:
        fw_primer = ""
        rv_primer = ""

        merge_pairs(data_dir="")

        trim_primers()

        length_filter()

        quality_filter()

        chimera_filter()

        frequency_filter()

        translation_filter()

        denoise_seqs()

        cluster_otus()

        get_taxon_assignments()


if barcode_type == "18S":
    pass

if barcode_type == "16S":
    pass

