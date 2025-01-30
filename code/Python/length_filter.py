import os

trimmed_suffix = "_trimmed.fastq"

def length_filter(data_dir:str, amplicon_length:int):
    """
    Filter sequences to match amplicon length.
    Inputs:
        data_dir: path to data directory as a string
        amplicon_length: fixed length of amplicon sequence.
    Outputs:
        Trimmed sequences as fasta files in subfolder of data directory
    """

    # make output directory for length-filtered files
    if "length_filtered" not in os.listdir(data_dir):
        os.makedirs(f"{data_dir}/length_filtered/")
    

    for data_file in os.listdir(f"{data_dir}/quality_filtered/"):
        if "fastq" in data_file:
            rm_counts = 0 # Keep track of removed lines (includes metadata lines at the moment)
            kepper_counts = 0 # Keep track of kept sequences

            log_file = open(f"{data_dir}/length_filtered/length_filter.log", "a")
            in_file = open(f"{data_dir}/trimmed/{data_file}", "r")
            out_file = open(f"{data_dir}/length_filtered/{data_file[0:-len(trimmed_suffix)]}.fasta", "a")

            for line in in_file:
                if len(line) == amplicon_length:
                    out_file.write(f">{line}")
                    kepper_counts += 1
                if len(line) != 142:
                    rm_counts += 1
            out_file.close()
            in_file.close()
            log_file.write(f"{data_file} had {rm_counts} lines removed and {kepper_counts} sequences were kept")
    log_file.close()


length_filter(data_dir = "../../data/test_data", amplicon_length = 142)