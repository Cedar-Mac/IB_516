import os, shutil
from collections import Counter


def frequency_filter(data_dir:str):
    """
    Some choices are made here.

    Currently only singleton reads present at a single site are removed.
    This is potentially not great.
    """

    # Remove old directory and files
    if "freq_filtered" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/freq_filtered/")

    # make output directory for frequency filtered files
    os.makedirs(f"{data_dir}/freq_filtered/")

    # make temp directory for sorted counts
    os.makedirs(f"{data_dir}/freq_filtered/temp")

    ones_list = []

    for file in os.listdir(f"{data_dir}/length_filtered/"):

        if "fasta" in file:
            # open files
            in_file = open(f"{data_dir}/length_filtered/{file}", "r")
            temp_file = open(f"{data_dir}/freq_filtered/temp/{file}", "a")

            # Read lines in file into list
            lines = in_file.readlines()

            # Count occurance of sequences
            seq_counts = Counter(lines[1::2])

            for key in seq_counts:
                if seq_counts[key] == 1:
                    ones_list.append(key)

            seq_counts_sorted = seq_counts.most_common()

            for seq in seq_counts_sorted:
                temp_file.write(f"{seq[1]}\n")
                temp_file.write(f"{seq[0]}")

            temp_file.close()
            in_file.close()
    
    # Ones list now contains all sequences only seen once.
    # Find sequences that only occur once and at a single site.
    multiple_sites = set()
    one_site_only = []
    for one in ones_list:
        if one not in multiple_sites:
            one_site_only.append(one)
            multiple_sites.add(one) 

    for tmp_file in os.listdir(f"{data_dir}/freq_filtered/temp/"):
        
        if "fasta" in file:

            with open(f"{data_dir}/freq_filtered/temp/{tmp_file}") as f:
                out_file = open(f"{data_dir}/freq_filtered/{tmp_file}", "a")

                for header, seq in zip(f, f):
                    if (seq not in one_site_only) or (header != "1"):
                        out_file.write(header)
                        out_file.write(seq)

            out_file.close()


    # Remove temp directory and files
    if "temp" in os.listdir(f"{data_dir}/freq_filtered"):
        shutil.rmtree(f"{data_dir}/freq_filtered/temp/")


frequency_filter("../../data/test_data")