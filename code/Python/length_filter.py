import os

trimmed_suffix = "_trimmed.fastq"

def length_filter(data_dir:str):
    # make output directory for length-filtered files
    if "length_filtered" not in os.listdir(data_dir):
        os.makedirs(f"{data_dir}/length_filtered/")
    

    for data_file in os.listdir(f"{data_dir}/quality_filtered/"):
        if "fastq" in data_file:
            rm_counts = 0
            kepper_counts = 0
            in_file = open(f"{data_dir}/trimmed/{data_file}", "r")
            out_file = open(f"{data_dir}/length_filtered/{data_file[0:-len(trimmed_suffix)]}.fasta", "a")
            for line in in_file:
                if len(line) == 142:
                    out_file.write(f">{line}")
                    kepper_counts += 1
                if len(line) != 142:
                    rm_counts += 1
            out_file.close()
            print(f"{in_file} had {rm_counts} removed sequences and {kepper_counts} retained sequences")


length_filter("../../data/test_data")

