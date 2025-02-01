import os, shutil
from collections import Counter


def frequency_filter(data_dir:str, min_seq_count:int, min_site_occurance:int):
    """
    Filter sequences based on their frequency of occurances within and between sites.

    Input:
        - data_dir: string with path to data directory
        - min_seq_count: Integer, minimum sequence occurance at a site
        - min_site_occurance: minimum number of sites a low abundance sequence
                            must occur at to be retained.

    Output: 
        - fasta files with size (seq count) in header in freq_filtered subdirectory.
        - Formated for input into DnoisE filtering algorithm.

    Details:
        A sequence occurance at a site is only dropped if two conditions are met:

        - The sequence occurs less than a minimum number of times (min_seq_count)
        - And the sequence occurs at less than a minimum number of sites (min_site_occurance).

        A temporary file for each site is created to group duplicate sequences and order them by abundance.
        Currently the only information stored in the sequence header is the sequence name (site_name_sequence_id)
        and count of the sequence ("size")
    """

    fasta_suffix = ".fasta"

    # Remove old directory and files
    if "freq_filtered" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/freq_filtered/")

    # make output directory for frequency filtered files
    os.makedirs(f"{data_dir}/freq_filtered/")

    # make temp directory for sorted counts
    os.makedirs(f"{data_dir}/freq_filtered/temp/")

    too_few_list = []

    for file in os.listdir(f"{data_dir}/chimera_filtered/"):

        if "fasta" in file:
            # open files
            in_file = open(f"{data_dir}/chimera_filtered/{file}", "r")
            temp_file = open(f"{data_dir}/freq_filtered/temp/{file}", "a")

            # Read lines in file into list
            lines = in_file.readlines()

            # Count occurance of sequences
            seq_counts = Counter(lines[1::2])

            for key in seq_counts:
                if seq_counts[key] < min_seq_count:
                    too_few_list.append(key)

            seq_counts_sorted = seq_counts.most_common()

            for seq in seq_counts_sorted:
                temp_file.write(f"{seq[1]}\n")
                temp_file.write(f"{seq[0]}")

            temp_file.close()
            in_file.close()
    
    # Ones list now contains all sequences seen a minimum number of times.
    # From these limited sequences, find ones that only occur at a minimum number of sites.
    site_occurances = Counter(too_few_list)

    i = 1

    for tmp_file in os.listdir(f"{data_dir}/freq_filtered/temp/"):
        
        if "fasta" in tmp_file:

            with open(f"{data_dir}/freq_filtered/temp/{tmp_file}") as f:
                out_file = open(f"{data_dir}/freq_filtered/{tmp_file}", "a")
                log_file = open(f"{data_dir}/freq_filtered/freq_filter.log", "a")

                for count, seq in zip(f, f):
                    if (site_occurances[seq] >= min_site_occurance) or (int(count) >= min_seq_count):
                        out_file.write(f">id={tmp_file[0:-len(fasta_suffix)]}_{i}; size={count.rstrip()};\n")
                        out_file.write(seq)
                        i += 1
                    else:
                        log_file.write(f"{seq} occured less than {min_seq_count} times at {tmp_file[0:-len(fasta_suffix)]} and was only present at {site_occurances[seq] - 1} other sites.\n")

                out_file.close()
    log_file.close()


    # Remove temp directory and files
    if "temp" in os.listdir(f"{data_dir}/freq_filtered"):
        shutil.rmtree(f"{data_dir}/freq_filtered/temp/")


if __name__ == "__main__":
    frequency_filter("../../data/test_data", min_seq_count=3, min_site_occurance=3)