import os, shutil
from itertools import pairwise

fastq_suffix = ".fastq"

def length_filter(data_dir:str, amplicon_length:int):
    """
    Filter sequences to match amplicon length.

    Inputs:
        - data_dir: path to data directory as a string
        - amplicon_length: fixed length of amplicon sequence.

    Outputs:
        - Trimmed sequences as fasta files in subfolder of data directory
    """

    # Remove old files and make output directory for length-filtered files
    if "length_filtered" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/length_filtered/")

    if "length_filtered" not in os.listdir(data_dir):
        os.makedirs(f"{data_dir}/length_filtered/")

    for data_file in os.listdir(f"{data_dir}/quality_filtered/"):
        if "fastq" in data_file:
            rm_counts = 0 # Keep track of removed lines (includes metadata lines at the moment)
            kepper_counts = 0 # Keep track of kept sequences

            log_file = open(f"{data_dir}/length_filtered/length_filter.log", "a")
            in_file = open(f"{data_dir}/quality_filtered/{data_file}", "r")
            out_file = open(f"{data_dir}/length_filtered/{data_file[0:-len(fastq_suffix)]}.fasta", "a")

            lines = in_file.readlines()

            # Iterate over pairs of lines.
            for i, pair in enumerate(pairwise(lines)):
                if (i % 4 == 0) and (len(pair[1].rstrip()) == amplicon_length): # if current line is a header and next line is a keeper sequence
                    out_file.write(f">{pair[0]}") # write header line
                if i % 4 == 1:
                    if len(pair[0].rstrip()) == amplicon_length: # if current line is a keeper sequence
                        out_file.write(f"{pair[0]}") # write keeper sequence
                        kepper_counts += 1
                    if len(pair[0]) != 142: # sequence was not correct length, keep track of skipped sequences.
                        rm_counts += 1
            out_file.close()
            in_file.close()
            log_file.write(f"{data_file} had {rm_counts} sequences removed and {kepper_counts} sequences were kept.\n\n")
    log_file.close()


if __name__ == "__main__":
    length_filter(data_dir = "../../data/test_data", amplicon_length = 142)